function [results, params] = prevalenceCore(a, P2, alpha)

% permutation-based prevalence inference using the minimum statistic, core
% 
% [results, params] = prevalenceCore(a, P2 = 1e6, alpha = 0.05)
%
% a:            three-dimensional array of test statistic values
%               (voxels x subjects x first-level permutations)
%               a(:, :, 1) must contain actual values
% P2:           number of second-level permutations to generate
% alpha:        significance level
% results:      per-voxel analysis results
%   .puGN         uncorrected p-values for global null
%   .pcGN         corrected p-values for global null
%   .puMN         uncorrected p-values for majority null
%   .pcMN         corrected p-values for majority null
%   .gamma0       prevalence lower bounds
%   .aTypical     typical values of test statistic where pcMN <= alpha
% params:        analysis parameters and properties
%   .V            number of voxels
%   .N            number of subjects
%   .P1           number of first-level permutations
%   .P2           number of second-level permutations actually generated
%   .alpha        significance level
%   .pcMNMin      smallest possible corrected p-value for majority null
%   .gamma0Max    largest possible prevalence lower bound
%
% See also prevalence.
%
%
% Copyright (C) 2016 Carsten Allefeld
%
% This program is free software: you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation, either version 3 of the License, or (at your
% option) any later version. This program is distributed in the hope that
% it will be useful, but without any warranty; without even the implied
% warranty of merchantability or fitness for a particular purpose. See the
% GNU General Public License <http://www.gnu.org/licenses/> for more details.

% default parameter values
if (nargin < 2) || isempty(P2)
    P2 = 1e6;
end
if (nargin < 3) || isempty(alpha)
    alpha = 0.05;
end

% init
[V, N, P1] = size(a); % V: voxels, N: subjects, P1: first-level permutations
fprintf('generating %d of %d second-level permutations\n', P2, P1 ^ N)
if P2 > P1 ^ N
    error('Monte Carlo implementation is inadequate!')
    % implement full enumeration of permutations?
end
fprintf('the computation can be stopped at any time by closing the output window\n\n')

% prepare plot window
fh = figure('Name', 'permutation-based prevalence inference using the minimum statistic');
text(0.5, 0.5, {'please wait for results', '', ...
    'close window to stop computation at any time'}, 'HorizontalAlignment', 'center')
axis off
drawnow

% generate second-level permutations
uRank = zeros(V, 1);
cRank = zeros(V, 1);
nPermsReport = 1;
tic
for j = 1 : P2
    % select first-level permutations
    if j == 1
        % select neutral permutations
        sp = ones(N, 1);
    else
        % select permutations randomly (Monte Carlo)
        sp = randi(P1, N, 1);
    end
    % indices of permutation values for each subject
    ind = sub2ind([N, P1], (1 : N)', sp);
    
    % test statistic: minimum across subjects
    m = min(a(:, ind), [], 2);
    % store result of neutral permutation (actual value) for each voxel
    if j == 1
        m1 = m;
    end
    
    % compare actual value with permutation value for each voxel separately,
    % determines uncorrected p-values for global null (see below)
    uRank = uRank + (m >= m1);          % part of Eq. 24
    % compare actual value at each voxel with maximum across voxels,
    % determines corrected p-values for global null (see below)
    cRank = cRank + (max(m) >= m1);     % Eq. 25 & part of Eq. 26

    % calibrate reporting interval to be between 2 and 5 seconds
    if (nPermsReport == 1) && (toc >= 2)
        nPermsReport = 10 ^ floor(log10(j)) * [1 2 5 10];
        nPermsReport = min(nPermsReport(nPermsReport >= j));
        nPermsReport = max(nPermsReport, 2);
    end
    
    % at intervals
    if ((nPermsReport > 1) && (mod(j, nPermsReport) == 0)) || (j == P2)
        drawnow
        stop = ~ishandle(fh);
        
        % compute results,
        % based on permutations performed so far: j plays the role of P2 
        % uncorrected p-values for global null hypothesis
        puGN = uRank / j;                               % part of Eq. 24
        % corrected p-values for global null hypothesis
        pcGN = cRank / j;                               % part of Eq. 26
        % significant voxels for global null hypothesis
        sigGN = (pcGN <= alpha);
        % * Step 5a: compute p-values for given prevalence bound
        % specifically 0.5, i.e for the majority null hypothesis
        % uncorrected p-values for majority null hypothesis
        puMN = ((1 - 0.5) * puGN .^ (1/N) + 0.5) .^ N;  % Eq. 19
        % corrected p-values for majority null hypothesis
        pcMN = pcGN + (1 - pcGN) .* puMN;
        % significant voxels for majority null hypothesis
        sigMN = (pcMN <= alpha);
        % lower bound on corrected p-values for majority null hypothesis
        puMNMin = ((1 - 0.5) * 1/j .^ (1/N) + 0.5) .^ N;
        pcMNMin = 1/j + (1 - 1/j) .* puMNMin;
        % * Step 5b: compute lower bounds for prevalence
        % lower bounds for prevalence
        alphac = (alpha - pcGN) ./ (1 - pcGN);          % Eq. 22
        gamma0 = (alphac .^ (1/N) - puGN .^ (1/N)) ./ (1 - puGN .^ (1/N)); % Eq. 23
        gamma0(alphac < puGN) = nan;                    % undefined       
        % upper bound for lower bounds
        alphacMax = (alpha - 1/j) / (1 - 1/j);          % Eq. 27
        gamma0Max = (alphacMax .^ (1/N) - 1/j .^ (1/N)) ./ (1 - 1/j .^ (1/N)); % Eq. 27
        
        % print summary
        fprintf('  %d of %d permutations,  %.1f of %.1f min\n', ...
            j, P2, toc / 60, toc / 60 * P2 / j)
        fprintf('    minimal rank\n')
        fprintf('      uncorrected:  %d,  reached at %d voxels\n', ...
            min(uRank), sum(uRank == min(uRank)))
        fprintf('      corrected:    %d,  reached at %d voxels\n', ...
            min(cRank), sum(cRank == min(cRank)))
        fprintf('    minimal p-value for global null\n')
        fprintf('      uncorrected:  %g\n', ...
            min(puGN))
        fprintf('      corrected:    %g\n', ...
            min(pcGN))
        fprintf('    minimal p-value for majority null\n')
        fprintf('      uncorrected:  %g\n', ...
            min(puMN))
        fprintf('      corrected:    %g\n', ...
            min(pcMN))
        fprintf('    number of voxels (of %d) at which\n', V)
        fprintf('      global null is rejected:      %d\n', ...
            sum(sigGN))
        fprintf('      majority null is rejected:    %d\n', ...
            sum(sigMN))
        fprintf('      prevalence bound is defined:  %d\n', ...
            sum(~isnan(gamma0)))
        fprintf('    largest prevalence bound:  %g\n', max(gamma0))
        fprintf('\n')
        
        % graphical display
        if stop
            fh = figure('Name', 'permutation-based prevalence inference using the minimum statistic');
        else
            % make figure current without getting in the way
            set(groot, 'CurrentFigure', fh)
            clf
        end
        % prevalence bounds
        subplot(2, 1, 1)
        plot([0.5, V + 0.5], 0.5 * [1 1], 'k')
        hold all
        plot([0.5, V + 0.5], gamma0Max * [1 1], 'r')
        plot(gamma0, 'b.')
        axis([0.5, V + 0.5, 0, 1])
        title('prevalence lower bounds')
        ylabel('\gamma_0')
        if V > 200
            set(gca, 'XTick', [])
        else
            set(gca, 'XTick', 1 : V)
            if V > 20
                set(gca, 'XTickLabel', [])
            end
        end
        if j < P2
            text(V + 0.5, 1, sprintf('%.0f %% ', j / P2 * 100), ...
                'Color', [0.8 0.8 0.8], 'HorizontalAlignment', 'right', ...
                'VerticalAlignment', 'top')
        end
        % majority null p-values
        subplot(2, 1, 2)
        semilogy([0.5, V + 0.5], alpha * [1 1], 'k')
        hold all
        plot([0.5, V + 0.5], pcMNMin * [1 1], 'r')
        plot(pcMN, '.b')
        xlim([0.5, V + 0.5])
        title('p-values majority null')
        xlabel('voxels')
        ylabel('p^*_N(m | \gamma \leq 0.5)')
        if V > 200
            set(gca, 'XTick', [])
        else
            set(gca, 'XTick', 1 : V)
            if V > 20
                set(gca, 'XTickLabel', [])
            end
        end
        drawnow

        if stop
            s = warning('query', 'backtrace');
            warning off backtrace
            warning('computation has been stopped')
            warning(s)
            P2 = j;
            break
        end
    end
    
end

% where majority null can be rejected, typical value of test statistic
aTypical = nan(V, 1);
aTypical(sigMN) = median(a(sigMN, :, 1), 2);

% collect return values
params = struct;
params.V = V;
params.N = N;
params.P1 = P1;
params.P2 = P2;
params.alpha = alpha;
params.pcMNMin = pcMNMin;
params.gamma0Max = gamma0Max;
results = struct;
results.puGN = puGN;
results.pcGN = pcGN;
results.puMN = puMN;
results.pcMN = pcMN;
results.gamma0 = gamma0;
results.aTypical = aTypical;

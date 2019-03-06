function [results, params] = prevalenceCore(a, P2, alpha)

% permutation-based prevalence inference using the minimum statistic, core
% implementation of the method proposed by Allefeld, Goergen and Haynes (2016)
%
% [results, params] = prevalenceCore(a, P2 = 1e6, alpha = 0.05)
%
% a:            three-dimensional array of test statistic values
%               (test units x subjects x first-level permutations)
%               a(:, :, 1) must contain actual values
% P2:           number of second-level permutations to generate
% alpha:        significance level
% results:      per-test unit analysis results
%   .puGN         uncorrected p-values for global null hypothesis         (Eq. 24)
%   .pcGN         corrected p-values for global null hypothesis           (Eq. 26)
%   .puMN         uncorrected p-values for majority null hypothesis       (Eq. 19)
%   .pcMN         corrected p-values for majority null hypothesis         (Eq. 21)
%   .gamma0u      uncorrected prevalence lower bounds                     (Eq. 20)
%   .gamma0c      corrected prevalence lower bounds                       (Eq. 23)
%   .aTypical     median values of test statistic where pcMN <= alpha     (Fig. 4b)
% params:        analysis parameters and properties
%   .V            number of test units
%   .N            number of subjects
%   .P1           number of first-level permutations
%   .P2           number of second-level permutations actually generated
%   .alpha        significance level
%   .puMNMin      smallest possible uncorrected p-value for majority H0
%   .pcMNMin      smallest possible corrected p-value for majority H0
%   .gamma0uMax   largest possible uncorrected prevalence lower bound
%   .gamma0cMax   largest possible corrected prevalence lower bound       (Eq. 27)
%
% 'Test units' may correspond to different locations (voxels, sensors)
% different time points, or both. The user must take care that the
% 'spatially extended version of the prevalence null' (section 'Information
% maps' of Allefeld, Goergen and Haynes 2016) does make sense for the given
% set of test units, where "in a small area" translates to "in a small
% subset of test units". Otherwise the *corrected* p-values for the
% majority null and the *corrected* prevalence lower bounds are wrong.
%
% The 'majority null hypothesis' referenced here is a special case of the
% prevalence null hypothesis (Eq. 17), where the critical value is gamma0 =
% 0.5. It describes the situation where there is no effect in the majority
% of subjects in the population. Rejecting it allows to infer that there is
% an effect in the majority of subjects in the population. aTypical is only
% defined where the (spatially extended) majority null hypothesis can be
% rejected. Compare Fig. 4b and 'What does it mean for an effect to be
% "typical" in the population?' in the Discussion of Allefeld, Goergen and
% Haynes (2016).
%
% The function opens a figure window that shows results based on the
% permutations generated so far and is regularly updated. If it is closed,
% the computation is stopped (not aborted) and results are returned
% based on the permutations so far.
%
% The window shows in the upper panels
%  -- p-values for the majority null hypothesis for all test units (blue dots),
%  -- the currently smallest possible p-value (red line),
%  -- and the significance threshold alpha (black line).
% In the left side uncorrected values, puMN and puMNMin,
% and on the right side corrected values, pcMN and pcMNMin.
%
% In the lower panels it shows
%  -- prevalence lower bounds for all test units (blue dots),
%  -- the currently largest possible prevalence lower bound (red line),
%  -- and the majority prevalence 0.5 (black line).
% On the left side uncorrected values, gamma0u and gamma0uMax,
% and on the right side corrected values, gamma0c and gamma0cMax.
%
% Upper and lower panels show complementary representations of the same
% result, corresponding to algorithm Step 5a and 5b, respectively (see
% paper).
%
% If many blue dots lie on the red line in either panel, this indicates
% that a more differentiated result might be obtained from more
% permutations. However, the amount of additional permutations needed for
% that might be computationally unfeasible.
%
% If P2 is larger than the number of possible second-level permutations,
% P1 ^ N, an error message recommends to set P2 = P1 ^ N. If this
% recommendation is followed, Monte Carlo estimation is replaced by a
% complete enumeration of all possible second-level permutations. In this
% case, stopping the computation may bias the result.
%
% See also prevalence.
%
%
% Copyright (C) 2016–2019 Carsten Allefeld
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

% prepare plot window
fh = figure('Name', 'permutation-based prevalence inference using the minimum statistic');
set(gcf, 'Position', get(gcf, 'Position') .* [1 1 0 0] + [0 0 800 580])
text(0.5, 0.5, {'please wait for results', '', ...
    'close window to stop computation at any time'}, ...
    'HorizontalAlignment', 'center')
axis off
drawnow

% init
[V, N, P1] = size(a);
% V: test units, N: subjects, P1: first-level permutations
if P2 > P1 ^ N
    fprintf('%d is larger than the number of possible second-level permutations\n', P2)
    fprintf('please restart with parameter P2 = %d\n', P1 ^ N)
    error('^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^')
end
enum = (P2 == P1 ^ N);
if enum
    fprintf('enumerating all %d second-level permutations\n\n', P2)
else
    fprintf('randomly generating %d of %d possible second-level permutations\n\n', ...
        P2, P1 ^ N)
end

% generate second-level permutations
uRank = zeros(V, 1);
cRank = zeros(V, 1);
nPermsReport = 1;
fprintf('the computation can be stopped by closing the output window\n\n')
tstart = tic;
for j = 1 : P2
    % select first-level permutations
    if enum     % complete enumeration
        % translate index of second-level permutation (j)
        % into indices of first-level permutations (is)
        jc = j - 1; % de2bi does the same
        is = nan(N, 1);
        for k = 1 : N
            is(k) = rem(jc, P1) + 1;
            jc = floor(jc / P1);
        end
    else        % Monte Carlo
        % randomly select permutations, except for first
        if j == 1
            is = ones(N, 1);
        else
            is = randi(P1, N, 1);
        end
    end
    
    % translate indices of first-level permutations to indices into a
    ind = sub2ind([N, P1], (1 : N)', is);
    
    % test statistic: minimum across subjects
    m = min(a(:, ind), [], 2);
    % store result of neutral permutation (actual value) for each test unit
    if j == 1
        m1 = m;
    end
    
    % compare actual value with permutation value for each test unit separately,
    % determines uncorrected p-values for global null hypothesis (see below)
    uRank = uRank + (m >= m1);          % part of Eq. 24
    % compare actual value at each test unit with maximum across units,
    % determines corrected p-values for global null hypothesis (see below)
    cRank = cRank + (max(m) >= m1);     % Eq. 25 & part of Eq. 26
    
    % calibrate reporting interval to be between 2 and 5 seconds
    if (nPermsReport == 1) && (toc(tstart) >= 2)
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
        % significant test units for global null hypothesis
        sigGN = (pcGN <= alpha);
        % * Step 5a: compute p-values for given prevalence bound
        % (here specifically gamma0 = 0.5, i.e the majority null hypothesis)
        % uncorrected p-values for majority null hypothesis
        puMN = ((1 - 0.5) * puGN .^ (1/N) + 0.5) .^ N;  % Eq. 19
        % corrected p-values for majority null hypothesis
        pcMN = pcGN + (1 - pcGN) .* puMN;               % Eq. 21
        % significant test units for majority null hypothesis
        sigMN = (pcMN <= alpha);
        % lower bound on corrected p-values for majority null hypothesis
        puMNMin = ((1 - 0.5) * 1/j .^ (1/N) + 0.5) .^ N;
        pcMNMin = 1/j + (1 - 1/j) .* puMNMin;
        % * Step 5b: compute prevalence lower bounds for given alpha
        % uncorrected prevalence lower bounds
        gamma0u = (alpha .^ (1/N) - puGN .^ (1/N)) ./ (1 - puGN .^ (1/N)); % Eq. 20
        gamma0u(alpha < puGN) = nan;                    % undefined
        % corrected prevalence lower bounds
        alphac = (alpha - pcGN) ./ (1 - pcGN);          % Eq. 22
        gamma0c = (alphac .^ (1/N) - puGN .^ (1/N)) ./ (1 - puGN .^ (1/N)); % Eq. 23
        gamma0c(puGN > alphac) = nan;                   % undefined
        % upper bound for lower bounds
        gamma0uMax = (alpha .^ (1/N) - 1/j .^ (1/N)) ./ (1 - 1/j .^ (1/N));
        alphacMax = (alpha - 1/j) / (1 - 1/j);          % Eq. 27
        gamma0cMax = (alphacMax .^ (1/N) - 1/j .^ (1/N)) ./ (1 - 1/j .^ (1/N)); % Eq. 27
        % The criterion for the corrected prevalence lower bound to be
        % defined, `puGN <= alphac`, is not equivalent to the significance
        % criterion for the GN, `pcGN <= alpha`, but is slightly more
        % conservative. A possible disagreement may be detected by
        % comparing the lines
        %   global null hypothesis is rejected
        % and
        %   prevalence bound is defined
        % in the diagnostic output. The two numbers should normally be
        % identical, but the second one can be smaller.
        
        % print summary
        fprintf('  %d of %d permutations\n', j, P2)
        fprintf('    minimal rank\n')
        fprintf('      uncorrected:  %d,  reached in %d test units\n', ...
            min(uRank), sum(uRank == min(uRank)))
        fprintf('      corrected:    %d,  reached in %d test units\n', ...
            min(cRank), sum(cRank == min(cRank)))
        fprintf('    minimal p-value for global null hypothesis\n')
        fprintf('      uncorrected:  %g\n', ...
            min(puGN))
        fprintf('      corrected:    %g\n', ...
            min(pcGN))
        fprintf('    minimal p-value for majority null hypothesis\n')
        fprintf('      uncorrected:  %g\n', ...
            min(puMN))
        fprintf('      corrected:    %g\n', ...
            min(pcMN))
        fprintf('    number of test units (of %d) in which (corrected)\n', V)
        fprintf('      global null hypothesis is rejected:    %d\n', ...
            sum(sigGN))
        fprintf('      majority null hypothesis is rejected:  %d\n', ...
            sum(sigMN))
        fprintf('      prevalence bound is defined:           %d\n', ...
            sum(~isnan(gamma0c)))
        fprintf('    largest corrected prevalence bound:  %g\n', max(gamma0c))
        fprintf('\n')
        
        % graphical display
        if stop
            fh = figure('Name', 'permutation-based prevalence inference using the minimum statistic');
            set(gcf, 'Position', get(gcf, 'Position') .* [1 1 0 0] + [0 0 800 580])
        else
            % make figure current without getting in the way
            set(groot, 'CurrentFigure', fh)
            clf
        end
        % majority null hypothesis p-values, uncorrected
        subplot(2, 2, 1)
        semilogy([0.5, V + 0.5], alpha * [1 1], 'k')
        hold all
        plot([0.5, V + 0.5], puMNMin * [1 1], 'r')
        plot(puMN, '.b')
        title({'p-values majority null hypothesis', 'uncorrected'})
        ylabel('p_N(m | \gamma \leq 0.5)')
        xdeco
        % majority null hypothesis p-values, corrected
        subplot(2, 2, 2)
        semilogy([0.5, V + 0.5], alpha * [1 1], 'k')
        hold all
        plot([0.5, V + 0.5], pcMNMin * [1 1], 'r')
        plot(pcMN, '.b')
        title({'p-values majority null hypothesis', 'corrected'})
        ylabel('p^*_N(m | \gamma \leq 0.5)')
        xdeco
        % prevalence bounds, uncorrected
        subplot(2, 2, 3)
        plot([0.5, V + 0.5], 0.5 * [1 1], 'k')
        hold all
        plot([0.5, V + 0.5], gamma0uMax * [1 1], 'r')
        plot(gamma0u, 'b.')
        title({'prevalence lower bounds', 'uncorrected'})
        ylabel('\gamma_0')
        ylim([0, 1])
        xdeco
        % prevalence bounds, corrected
        subplot(2, 2, 4)
        plot([0.5, V + 0.5], 0.5 * [1 1], 'k')
        hold all
        plot([0.5, V + 0.5], gamma0cMax * [1 1], 'r')
        plot(gamma0c, 'b.')
        title({'prevalence lower bounds', 'corrected'})
        ylabel('\gamma_0^*')
        ylim([0, 1])
        xdeco
        % progress
        if j < P2
            annotation(gcf, 'textbox', [0 0 1 1], 'Color', [1 0.65 0], ...
                'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
                'String', sprintf(' %.0f %%          %.1f of %.1f min', ...
                j / P2 * 100, toc(tstart) / 60, toc(tstart) / 60 * P2 / j))
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

% where majority null hypothesis can be rejected, typical value of test statistic
aTypical = nan(V, 1);
aTypical(sigMN) = median(a(sigMN, :, 1), 2);

% collect return values
params = struct;
params.V = V;
params.N = N;
params.P1 = P1;
params.P2 = P2;
params.alpha = alpha;
params.puMNMin = puMNMin;
params.pcMNMin = pcMNMin;
params.gamma0uMax = gamma0uMax;
params.gamma0cMax = gamma0cMax;
results = struct;
results.puGN = puGN;
results.pcGN = pcGN;
results.puMN = puMN;
results.pcMN = pcMN;
results.gamma0u = gamma0u;
results.gamma0c = gamma0c;
results.aTypical = aTypical;


    function xdeco
        % common plot x-axis decorations
        xlabel('test units (voxels)')
        xlim([0.5, V + 0.5])
        if V > 200
            set(gca, 'XTick', [])
        else
            set(gca, 'XTick', 1 : V)
            if V > 20
                set(gca, 'XTickLabel', [])
            end
        end
    end
end

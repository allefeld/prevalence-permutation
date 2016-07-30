function prevalence(ifn, P2, ofn, alpha)
if nargin == 0, test_prevalence, return, end    % HACK

% permutation-based prevalence inference
%
% gamma0max = prevalence(inputfilenames, P2 = 1e6, outputfilename = 'prevalence', alpha = 0.05)
%
% inputfilenames:   cell array of input image filenames, subjects x permutations
% P2:               number of second-level permutations to perform
% outputfilename:   output image filename
% alpha:            significance level
% gamma0max:        theoretical upper bound on gamma0
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

if (nargin < 2) || isempty(P2)
    P2 = 1e6;
end
if (nargin < 3) || isempty(ofn)
    ofn = 'prevalence';
end
if nargin < 4
    alpha = 0.05;
end

fprintf('\n*** permutation-based prevalence inference ***\n\n')


%% load and prepare accuracies

[N, P1] = size(ifn);
fprintf('loading data: %d subjects, %d first-level permutations\n', N, P1)

% load test statistic images
a = cell(N, P1);    %  how to initialize vol?
for k = 1 : N
    fprintf('  subject #%d ', k)
    for i = 1 : P1
        gz = ~isempty(regexp(ifn{k, i}, '.gz$', 'once')); % gzipped image?
        if gz
            tfn = gunzip(ifn{k, i}, tempdir);   % gunzip to temporary file
            ifn{k, i} = tfn{1};
        end
        vol(k, i) = spm_vol(ifn{k, i});                                     %#ok<AGROW>
        Y = spm_read_vols(vol(k, i));
        a{k, i} = Y(:);
        if gz
            delete(tfn{1})
        end
        fprintf('.')
    end
    fprintf('\n')
end
spm_check_orientations(vol(:));
% a is now a cell array of size N x P1, where each cell contains voxel
% values in one column vector
a = cell2mat(reshape(a, [1, N, P1]));
% a is now a matrix of size (number of voxels) x N x P1
fprintf('%d voxels, ', size(a, 1));

% determine mask from data; out-of-mask voxels may be NaN or 0
mask = all(all(~isnan(a), 2), 3) & ~any(all(a == 0, 3), 2);
% truncate data to in-mask voxels
a = a(mask, :, :);
% a is now a matrix of size (number of in-mask voxels) x N x P1
fprintf('%d in-mask\n', size(a, 1));
fprintf('\n')


%% perform prevalence inference

prevalence_compute(a, P2, alpha)


%% save results

data = nan(size(mask));
data(mask) = gamma0;
saveMRImage(data, [ofn '_gamma0.nii'], vol.mat, 'prevalence map')
data = nan(size(mask));
data(mask) = at;
saveMRImage(data, [ofn '_typical.nii'], vol.mat, 'typical map')
saveMRImage(uint8(mask), [ofn '_mask.nii'], vol.mat, 'prevalence map mask')
% stored mask values end up to be 1.00000005913898 instead of 1 – why?


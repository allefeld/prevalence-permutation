function prevalence(ifn, P2, alpha, prefix)

% permutation-based prevalence inference, for MR image files
%
% prevalence(ifn, P2 = 1e6, alpha = 0.05, prefix = 'prevalence_')
%
% ifn:     cell array of input image filenames of test statistic values
%          (subjects x permutations)
%          images in ifn(:, 1) must contain actual values
% P2:      number of second-level permutations to perform
% prefix:  prefix for output files
% alpha:   significance level
%
% Results are written to files that begin with the string given in
% `prefix`. By including a path in `prefix`, the output file directory can
% be specified; otherwise files are written to the current directory.
%
% With the default prefix, analysis results are written to the files
%   prevalence_puGN.nii, prevalence_pcGN.nii, prevalence_puMN.nii,
%   prevalence_pcMN.nii, prevalence_gamma0.nii, & prevalence_aTypical.nii.
% Additionally, the brain mask is written to prevalence_mask.nii, and
% analysis parameters and properties to prevalence_param.mat.
% For a detailed explanation, see output parameters `results` and `params`
% of prevalence_compute.
%
% See also prevalence_compute.
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

if nargin == 0, test_prevalence, return, end    % HACK

if (nargin < 2) || isempty(P2)
    P2 = 1e6;
end
if (nargin < 3) || isempty(prefix)
    prefix = 'prevalence_';
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

[results, param] = prevalence_compute(a, P2, alpha);


%% save results

data = nan(size(mask));
data(mask) = gamma0;
saveMRImage(data, [prefix 'gamma0.nii'], vol.mat, 'prevalence map')
data = nan(size(mask));
data(mask) = at;
saveMRImage(data, [prefix 'typical.nii'], vol.mat, 'typical map')
saveMRImage(uint8(mask), [prefix '_mask.nii'], vol.mat, 'prevalence map mask')
% stored mask values end up to be 1.00000005913898 instead of 1 – why?


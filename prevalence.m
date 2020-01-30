function prevalence(ifn, P2, alpha, prefix, mask)

% permutation-based prevalence inference using the minimum statistic
%
% prevalence(ifn, P2 = 1e6, alpha = 0.05, prefix = 'prevalence_', mask)
%
% ifn:     cell array of input image filenames of test statistic values
%          (subjects x permutations)
%          images in ifn(:, 1) must contain actual values
% P2:      number of second-level permutations to perform
% alpha:   significance level
% prefix:  prefix for output files
% mask:    brain mask
%
% Results are written to files that begin with the string given in
% `prefix`. By including a path in `prefix`, the output file directory can
% be specified; otherwise files are written to the current directory.
%
% With the default prefix, analysis results are written to the files
%   prevalence_puGN.nii, prevalence_pcGN.nii, prevalence_puMN.nii,
%   prevalence_pcMN.nii, prevalence_gamma0c.nii, & prevalence_aTypical.nii.
% Additionally, the brain mask is written to prevalence_mask.nii, and
% analysis parameters and properties to prevalence_params.mat.
% For a detailed explanation, see output parameters `results` and `params`
% of prevalenceCore.
%
% See also prevalenceCore
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
if (nargin < 3) || isempty(prefix)
    prefix = 'prevalence_';
end
if nargin < 4
    alpha = 0.05;
end

fprintf('\n*** permutation-based prevalence inference using the minimum statistic ***\n\n')

% load and prepare test statistic data
[N, P1] = size(ifn);
fprintf('loading data: %d subjects, %d first-level permutations\n', N, P1)
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
% a is now a cell array of size N x P1, containing voxel value columns vectors
a = cell2mat(reshape(a, [1, N, P1]));
% a is now a matrix of size (number of voxels) x N x P1
fprintf('%d voxels, ', size(a, 1));
if nargin < 5
    % determine mask from data; out-of-mask voxels may be NaN or 0
    mask = all(all(~isnan(a), 2), 3) & ~any(all(a == 0, 3), 2);
end
% truncate data to in-mask voxels
a = a(mask, :, :);
% a is now a matrix of size (number of in-mask voxels) x N x P1
fprintf('%d in-mask\n', size(a, 1));
fprintf('\n')

% perform prevalence inference
[results, params] = prevalenceCore(a, P2, alpha);                           %#ok<ASGLU>

% save
fprintf('saving results and parameters\n')
% per-voxel analysis results
f = fields(results);
for fi = 1 : numel(f)
    data = nan(vol(1).dim);
    data(mask) = getfield(results, f{fi});                                  %#ok<GFLD>
    spmWriteImage(data, [prefix f{fi} '.nii'], vol(1).mat, 'descrip', f{fi})
end
% brain mask
spmWriteImage(reshape(mask, vol(1).dim), [prefix 'mask.nii'], vol(1).mat, ...
    'descrip', 'prevalence brain mask')
% analysis parameters and properties
save([prefix 'params'], '-struct', 'params')

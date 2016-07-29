function gamma0max = prevalence(inputfilenames, P2, outputfilename, alpha)

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
if (nargin < 3) || isempty(outputfilename)
    outputfilename = 'prevalence';
end
if nargin < 4
    alpha = 0.05;
end
[N, P1] = size(inputfilenames);


%% load and prepare accuracies

fprintf('\n*** prevalence ***\n\n')
fprintf('loading data\n')

% load accuracy images
a = cell(N, P1);
for k = 1 : N
    fprintf('  subject #%d: ', k)
    for i = 1 : P1
        vol = spm_vol(inputfilenames{k, i});
        Y = spm_read_vols(vol);
        a{k, i} = Y(:);
        fprintf('.')
    end
    fprintf('\n')
end
% a is now a cell array of size N x P1, where each cell contains voxel
% values in one column vector
a = cell2mat(reshape(a, [1, N, P1]));
% a is now a matrix of size (number of voxels) x N x P1

% determine mask from data; out-of-mask voxels may be NaN or 0
mask = all(all(~isnan(a), 2), 3) & ~any(all(a == 0, 3), 2);
mask = reshape(mask, size(Y));

% truncate data to in-mask voxels
a = a(mask, :, :);
V = sum(mask(:));
% a is now a matrix of size V x N x P1


%% call prevalence_compute.m


%% save results

data = nan(size(mask));
data(mask) = gamma0;
saveMRImage(data, [outputfilename '_gamma0.nii'], vol.mat, 'prevalence map')
data = nan(size(mask));
data(mask) = at;
saveMRImage(data, [outputfilename '_typical.nii'], vol.mat, 'typical map')
saveMRImage(uint8(mask), [outputfilename '_mask.nii'], vol.mat, 'prevalence map mask')
% stored mask values end up to be 1.00000005913898 instead of 1 – why?


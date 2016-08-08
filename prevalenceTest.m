% test the implementation of 'permutation-based prevalence inference
% using the minimum statistic'
%
% This script uses the data set 'cichy-2011-category-smoothedaccuracy' to
% test the present implementation of prevalence inference, and to check
% whether changes of the code change the result. The data are not included
% with the implementation, but can be downloaded from
%   https://github.com/allefeld/cichy-2011-category-smoothedaccuracy/releases/tag/v1.0.0
% Unzip the archive within this directory, so that its contents are in the
% subdirectory cichy-2011-category-smoothedaccuracy-1.0.0.
%
% This script can also be taken as an example for how to use the implementation.
%
% With parameter P2 changed to 1e7, it should reproduce the results shown
% in Allefeld, Goergen & Haynes (2016), Fig 4a and b, except for minor
% differences due to Monte Carlo estimation. 
%
% WARNING: The script resets the random number generator to a fixed seed in
% order to obtain reproducible results. Make sure to reinitialize it again
% to a time-based seed afterwards: rng('shuffle')
%
% With Matlab 8.5.0 (R2015a) and SPM8 r4290, this script generates image
% files with the following md5 checksums: 
%   5e9967acd0fb7e63ee9d74839029519f  prevalence_aTypical.nii
%   c1d1b791aa080e79e88105d36681d94b  prevalence_gamma0c.nii
%   c830ddf72c6a30b0f0e228f2d9589232  prevalence_gamma0u.nii
%   afb1c042ee7329832b2b895d576664b3  prevalence_mask.nii
%   1b8a5e4edf0270f23e6bed54b9ae4b1f  prevalence_pcGN.nii
%   d23d16395900d13c2e7b8976afc0edc5  prevalence_pcMN.nii
%   8ac8c45ee7314adfc16172c451442a22  prevalence_puGN.nii
%   d98ef36d1815e2632d0cd38358f95827  prevalence_puMN.nii
%
%
% This file is part of v1.1.0 of prevalence-permutation, see
% https://github.com/allefeld/prevalence-permutation/releases
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

clear

% check if data is present
datadir = 'cichy-2011-category-smoothedaccuracy-1.0.0';
if ~exist(datadir, 'dir')
    fprintf('please download test data, see `help prevalenceTest`\n')
    return
end

% collect input image filenames
N = 12;
P1 = 16;
input = cell(N, P1);
for k = 1 : N
    for i = 1 : P1
        input{k, i} = sprintf([datadir '/%02d/sa_C0002_P%04d.nii.gz'], k, i);
    end
end

% make temporary directory for results
tn = tempname;
mkdir(tn);

% initialize the random number generator for reproducibilty
% In normal use, the rng should be initialized by `rng('shuffle')` on
% Matlab startup and not changed during the session.
rng(42)

% perform prevalence inference
P2 = 200000;
alpha = 0.05;
prefix = [tn '/prevalence_'];
prevalence(input, P2, alpha, prefix)
rng('shuffle')

% compute checksums of results
fprintf('\nchecksums of results:\n')
if isunix
    system(['md5sum ' tn '/*.nii']);
else
    system(['certutil -hashfile ' tn '/*.nii MD5']); % DOES THIS WORK?
end

fprintf('\nconsider deleting the directory %s and its contents\n', tn)

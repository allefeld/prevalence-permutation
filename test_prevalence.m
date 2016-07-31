% test the implementation of permutation-based prevalence inference
%
% This script uses the data contained in the subdirectory
% cichy-2011-category-smoothedaccuracy to test the present implementation
% of prevalence inference, and to check whether changes of the code change
% the result.
%
% It can also be taken as an example for how to use the implementation.
%
% With parameter P2 changed to 1e7, this script should reproduce the
% results shown in Allefeld, Görgen & Haynes (2016), Fig 4a and b, except
% for minor differences due to Monte Carlo estimation.
%
% On 2016-7-31 with Matlab 8.5.0 (R2015a) and SPM8 r4290 under Linux 3.16,
% this script generated image files with the following md5 checksums:
%   5e9967acd0fb7e63ee9d74839029519f  prevalence_aTypical.nii
%   19a5df882deeb5c857ce3345d194e017  prevalence_gamma0.nii
%   afb1c042ee7329832b2b895d576664b3  prevalence_mask.nii
%   1b8a5e4edf0270f23e6bed54b9ae4b1f  prevalence_pcGN.nii
%   d23d16395900d13c2e7b8976afc0edc5  prevalence_pcMN.nii
%   8ac8c45ee7314adfc16172c451442a22  prevalence_puGN.nii
%   d98ef36d1815e2632d0cd38358f95827  prevalence_puMN.nii

clear, close all

% collect input image filenames
N = 12;
P1 = 16;
ifnPat = 'cichy-2011-category-smoothedaccuracy/%02d/sa_C0002_P%04d.nii.gz';
ifn = cell(N, P1);
for k = 1 : N
    for i = 1 : P1
        ifn{k, i} = sprintf(ifnPat, k, i);
    end
end

% make temporary directory for results
tn = tempname;
mkdir(tn);

% initialize the random number generator for reproducibilty, to always
% obtain exactly the same Monte Carlo result
% in normal use, the rng should be initialized by `rng('shuffle')` on
% startup.
rng(42)

% perform prevalence inference
P2 = 200000;
alpha = 0.05;
prefix = [tn '/prevalence_'];
prevalence(ifn, P2, alpha, prefix)

% compute checksums of results, needs md5sum to be installed on the system
fprintf('\nchecksums of results:\n')
system(['md5sum ' tn '/*.nii']);

fprintf('\nconsider deleting the directory %s and its contents\n', tn)

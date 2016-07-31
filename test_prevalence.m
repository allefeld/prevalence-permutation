clear, close all

N = 12;
P1 = 16;
ifn = cell(N, P1);

for k = 1 : N
    for i = 1 : P1
        ifn{k, i} = sprintf('cichy-2011-category-smoothedaccuracy/%02d/sa_C0002_P%04d.nii.gz', k, i);
    end
end

% initialize rng for reproducibility

tn = tempname;
mkdir(tn);
prefix = fullfile(tn, 'prevalence_');

prevalence(ifn, 2000, 0.05, prefix)


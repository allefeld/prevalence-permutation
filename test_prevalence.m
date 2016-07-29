clear, close all

N = 12;
P1 = 16;
ifn = cell(N, P1);

for i = 1 : N
    for j = 1 : P1
        ifn{i, j} = sprintf('cichy-2011-smoothedaccuracy/%02d/sa_C0002_P%04d.nii.gz', i, j);
    end
end

for i = 1 : N
    for j = 1 : P1
        fn = gunzip(ifn{i, j}, tempname);
        ifn{i, j} = fn{1};
    end
end

ifn

% v = spm_vol(fn);

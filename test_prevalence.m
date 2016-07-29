clear, close all

n = 20;
N = 12;
P1 = 5;
a = randn(n, N, P1);
a(:, :, 1) = a(:, :, 1) + 2;

prevalence_compute(a)

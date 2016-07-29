clear, close all

V = 1000;
N = 12;
P1 = 5;
a = randn(V, N, P1);
a(:, :, 1) = bsxfun(@plus, a(:, :, 1), linspace(0, 2, V)');

prevalence_compute(a)

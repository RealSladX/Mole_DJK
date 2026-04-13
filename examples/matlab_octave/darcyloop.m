close all; clc;

alphas = [1, 10, 100, 1000];
ms = [16, 32, 64, 128, 256];

maxErr = zeros(length(ms), length(alphas));
relL2  = zeros(length(ms), length(alphas));

for ia = 1:length(alphas)
    alpha = alphas(ia);

    for im = 1:length(ms)
        m = ms(im);
        n = m;

        [errMax, errL2] = darcy(m, n, alpha);

        maxErr(im, ia) = errMax;
        relL2(im, ia)  = errL2;

        fprintf('done: alpha=%g, m=%d, max=%e, relL2=%e\n', ...
            alpha, m, errMax, errL2);
    end
end

fprintf('\n================ MAX ABS ERROR ================\n');
fprintf('%8s', 'm');
for ia = 1:length(alphas)
    fprintf('%18s', sprintf('alpha=%g', alphas(ia)));
end
fprintf('\n');

for im = 1:length(ms)
    fprintf('%8d', ms(im));
    for ia = 1:length(alphas)
        fprintf('%18.6e', maxErr(im, ia));
    end
    fprintf('\n');
end

fprintf('\n================ RELATIVE L2 ERROR ================\n');
fprintf('%8s', 'm');
for ia = 1:length(alphas)
    fprintf('%18s', sprintf('alpha=%g', alphas(ia)));
end
fprintf('\n');

for im = 1:length(ms)
    fprintf('%8d', ms(im));
    for ia = 1:length(alphas)
        fprintf('%18.6e', relL2(im, ia));
    end
    fprintf('\n');
end

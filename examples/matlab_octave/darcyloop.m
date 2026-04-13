close all; clc;

alphas = [1, 10, 100, 1000];
ms = [16, 32, 64, 128, 256];

% maxErr = zeros(length(ms), length(alphas));
% relL2  = zeros(length(ms), length(alphas));
% qxErr = zeros(length(ms), length(alphas));
% qyErr = zeros(length(ms), length(alphas));
% uMin  = zeros(length(ms), length(alphas));
% uMax  = zeros(length(ms), length(alphas));
% for ia = 1:length(alphas)
%     alpha = alphas(ia);
%
%     for im = 1:length(ms)
%         m = ms(im);
%         n = m;
%
%         [errMax, errL2,  qxErr(im, ia), qyErr(im, ia), ...
%          uMin(im, ia), uMax(im, ia)] = darcy(m, n, alpha);
%
%         maxErr(im, ia) = errMax;
%         relL2(im, ia)  = errL2;
%
%         fprintf('done: alpha=%g, m=%d, max=%e, relL2=%e\n', ...
%             alpha, m, errMax, errL2);
%     end
% end

% fprintf('\n================ MAX ABS ERROR ================\n');
% fprintf('%8s', 'm');
% for ia = 1:length(alphas)
%     fprintf('%18s', sprintf('alpha=%g', alphas(ia)));
% end
% fprintf('\n');
%
% for im = 1:length(ms)
%     fprintf('%8d', ms(im));
%     for ia = 1:length(alphas)
%         fprintf('%18.6e', maxErr(im, ia));
%     end
%     fprintf('\n');
% end
%
% fprintf('\n================ RELATIVE L2 ERROR ================\n');
% fprintf('%8s', 'm');
% for ia = 1:length(alphas)
%     fprintf('%18s', sprintf('alpha=%g', alphas(ia)));
% end
% fprintf('\n');
%
% for im = 1:length(ms)
%     fprintf('%8d', ms(im));
%     for ia = 1:length(alphas)
%         fprintf('%18.6e', relL2(im, ia));
%     end
%     fprintf('\n');
% end
uErrMax = zeros(length(ms), length(alphas)); 
uErr  = zeros(length(ms), length(alphas));
qxErr = zeros(length(ms), length(alphas));
qyErr = zeros(length(ms), length(alphas));
uMin  = zeros(length(ms), length(alphas));
uMax  = zeros(length(ms), length(alphas));

for ia = 1:length(alphas)
    alpha = alphas(ia);

    for im = 1:length(ms)
        m = ms(im);
        n = m;

        [uErrMax(im, ia), uErr(im, ia), qxErr(im, ia), qyErr(im, ia), ...
         uMin(im, ia), uMax(im, ia)] = darcy(m, n, alpha);

        fprintf('done: alpha=%g, m=%d\n', alpha, m);
    end
end
for ia = 1:length(alphas)
    fprintf('\nTable for alpha = %g\n', alphas(ia));
    fprintf('%8s %16s %8s %16s %8s %16s %8s %12s %12s %12s\n', ...
        'h', '||u-u^h||', 'Rate', '||qx-qx^h||', 'Rate', ...
        '||qy-qy^h||', 'Rate', 'u_min', 'u_max', "u_err_max");

    for im = 1:length(ms)
        h = 1 / ms(im);

        if im == 1
            uRate  = NaN;
            qxRate = NaN;
            qyRate = NaN;
            fprintf('%8s %16.3e %8s %16.3e %8s %16.3e %8s %12.3e %12.3e %12.3e\n', ...
                sprintf('1/%d', ms(im)), ...
                uErr(im, ia), '-', ...
                qxErr(im, ia), '-', ...
                qyErr(im, ia), '-', ...
                uMin(im, ia), uMax(im, ia), uErrMax(im, ia));
        else
            uRate  = log(uErr(im-1, ia)  / uErr(im, ia))  / log(2);
            qxRate = log(qxErr(im-1, ia) / qxErr(im, ia)) / log(2);
            qyRate = log(qyErr(im-1, ia) / qyErr(im, ia)) / log(2);

            fprintf('%8s %16.3e %8.2f %16.3e %8.2f %16.3e %8.2f %12.3e %12.3e %12.3e\n', ...
                sprintf('1/%d', ms(im)), ...
                uErr(im, ia),  uRate, ...
                qxErr(im, ia), qxRate, ...
                qyErr(im, ia), qyRate, ...
                uMin(im, ia), uMax(im, ia), uErrMax(im, ia));
        end
    end
end

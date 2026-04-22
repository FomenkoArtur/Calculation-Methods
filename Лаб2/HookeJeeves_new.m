function [x_min, f_min, t, n, history] = HookeJeeves_new(x0, E, f)
    tic;
    lambda = 1.0;
    gamma = 0.5;
    beta = 2.0;
    lambda_max = 10.0;
    x_base = x0;
    f_prev = f(x_base);
    n = 0;
    dim = length(x0);
    kmax = 1e5 * dim;
 
    history = zeros(kmax + 1, dim);
    history(1, :) = x0(:)';
    hist_idx = 1;
 
    while (lambda >= E) && (n <= kmax)
        n = n + 1;
        x_cur = x_base;
        f_cur = f_prev;
 
        for i = 1:dim
            e    = zeros(1, dim);
            e(i) = 1;
            x_plus  = x_cur + lambda * e;
            x_minus = x_cur - lambda * e;
            f_plus  = f(x_plus);
            f_minus = f(x_minus);
 
            if f_plus < f_cur
                x_cur = x_plus;
                f_cur = f_plus;
            elseif f_minus < f_cur
                x_cur = x_minus;
                f_cur = f_minus;
            end
        end
 
        if f_cur < f_prev
            x_prev = x_base;
            x_base = x_cur;
            f_prev = f_cur;
 
            hist_idx = hist_idx + 1;
            history(hist_idx, :) = x_base;
 
            while n <= kmax
                n = n + 1;
                x_pattern = x_base + (x_base - x_prev);
                f_pattern = f(x_pattern);
 
                if f_pattern < f_prev
                    x_prev = x_base;
                    x_base = x_pattern;
                    f_prev = f_pattern;
                    lambda = min(lambda * beta, lambda_max);
 
                    hist_idx = hist_idx + 1;
                    history(hist_idx, :) = x_base;
                else
                    lambda = gamma * lambda;
                    break;
                end
            end
        else
            lambda = gamma * lambda;
        end
    end
 
    history = history(1:hist_idx, :);
 
    x_min = x_base;
    f_min = f(x_min);
    t = toc;
end
function [x_min, f_min, t, n, history] = HookeJeeves(x0, E, f)
    tic;
    lambda = 1.0;
    gamma = 0.5;
    x_base = x0;
    f_prev = f(x_base);
    n = 0;
    dim = length(x0);
    kmax = 1e5*dim;
 
    history = zeros(kmax + 1, dim);
    history(1, :) = x0(:)';
    hist_idx = 1;
 
    while (lambda >= E) && (n <= kmax)
        n = n + 1;
        f_best = f_prev;
        x_best = x_base;
 
        for i = 1:dim
            e = zeros(1, dim);
            e(i) = 1;
            x_minus = x_base - lambda * e;
            x_plus  = x_base + lambda * e;
            f_minus = f(x_minus);
            f_plus  = f(x_plus);
 
            if f_minus < f_best
                f_best = f_minus;
                x_best = x_minus;
            end
            if f_plus < f_best
                f_best = f_plus;
                x_best = x_plus;
            end
        end
 
        if f_best >= f_prev
            lambda = gamma * lambda;
        else
            f_prev = f_best;
            x_base = x_best;
 
            hist_idx = hist_idx + 1;
            history(hist_idx, :) = x_base;
        end
    end
 
    history = history(1:hist_idx, :);
 
    x_min = x_base;
    f_min = f(x_min);
    t = toc;
end
function [x_min, f_min, t, n, history] = NewtonRaphsonMethod(x0, E, f)
    tic;
    n = 0;
    x = x0(:);
    dim = length(x0);
    kmax = 1e4 * dim;
    h = 5;


    history = zeros(kmax + 1, dim);
    history(1, :) = x';
    hist_idx = 1;

    g = gradient(f, x, E/10);
    g = g(:);
    
    H_val = hessian(f, x, E*10);
    H_inv = H_val^(-1);

    while n <= kmax
        n = n + 1;
        
        p = -H_inv * g;
        
        func = @(alpha) f(x + alpha * p);
        alpha = GoldenSection(0, h, E, func);
        
        x = x + alpha * p;
        
        hist_idx = hist_idx + 1;
        history(hist_idx, :) = x';

        g = gradient(f, x, E/10);
        g = g(:);
        H_val = hessian(f, x, E*10);
        H_inv = H_val^(-1);

        grad_norm = norm(g);

        if grad_norm < E
            break;
        end
    end

    x_min = x;
    f_min = f(x_min);
    t = toc;
    history = history(1:hist_idx, :);
end
function [x_min, f_min, t, n, history] = BroydenMethod(x0, E, f)
    tic;
    n = 0;
    x = x0(:);
    dim = length(x0);
    kmax = 1e4 * dim;
    h = 0.025;


    history = zeros(kmax + 1, dim);
    history(1, :) = x';
    hist_idx = 1;

    A = eye(dim);
    g = gradient(f, x, E/10);
    g = g(:);
    
    while n <= kmax
        n = n + 1;
        
        p = -A * g;

        if mod(n, dim + 1) == 0
            A = eye(dim);
        end
        
        func = @(alpha) f(x + alpha * p);
        alpha = GoldenSection(0, h, E, func);
        
        x_old = x;
        g_old = g;
        
        x = x + alpha * p;
        
        hist_idx = hist_idx + 1;
        history(hist_idx, :) = x';
        
        g = gradient(f, x, E/10);
        g = g(:);
        
        dx = x - x_old;
        dg = g - g_old;
        
        u = dx - A * dg;
        denom = u' * dg;
        
        if abs(denom) > 1e-8 * norm(dg) * norm(u)
            A = A + (u * u') / denom;
        end

        if norm(g) < E
            break
        end
    end
    
    x_min = x;
    f_min = f(x_min);
    t = toc;
    history = history(1:hist_idx, :);
end
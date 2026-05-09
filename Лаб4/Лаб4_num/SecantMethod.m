function [x_min, f_min, t, n, history] = SecantMethod(x0, E, f)
    tic;
    x = x0(:);
    dim = length(x0);
    kmax = 1e4 * dim;


    history = zeros(kmax + 1, dim);
    history(1, :) = x';
    hist_idx = 1;

    g = gradient(f, x, E/10);
    g = g(:);
    
    x_prev = x;
    g_prev = g;
    x = x - 0.01 * g;
    g = gradient(f, x, E/10);
    g = g(:);
    n = 1;
    
    hist_idx = hist_idx + 1;
    history(hist_idx, :) = x';

    while n <= kmax
        n = n + 1;
        
        x_new = zeros(dim, 1);
        
        for j = 1:dim
            if abs(g(j) - g_prev(j)) > 1e-10
                H = (x(j) - x_prev(j)) / (g(j) - g_prev(j));
                x_new(j) = x(j) - H * g(j);
            else
                x_new(j) = x(j) - 0.01 * g(j);
            end
        end

        x_prev = x;
        g_prev = g;
        x = x_new;
        
        hist_idx = hist_idx + 1;
        history(hist_idx, :) = x';
        
        g = gradient(f, x, E/10);
        g = g(:);

        if norm(g) < E
            break;
        end
    end

    x_min = x;
    f_min = f(x_min);
    t = toc;
    history = history(1:hist_idx, :);
end
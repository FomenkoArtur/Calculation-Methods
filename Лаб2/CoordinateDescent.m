function [x_min, f_min, t, n, history] = CoordinateDescent(x0, E, f)
    tic;
    h = 1.0;
    x_current = x0;
    n = 0;
    dim = length(x0);
    kmax = 1e5*dim;
 
    history = zeros(kmax + 1, dim);
    history(1, :) = x0(:)';
    hist_idx = 1;
 
    while n <= kmax
        n = n + 1;
        x_start = x_current;
 
        for j = 1:dim
            a = x_current(j) - h;
            b = x_current(j) + h;
            func = @(x) f_single_coord(f, x_current, j, x);
            x_opt = GoldenSection(a, b, E, func);
            x_current(j) = x_opt;
        end
 
        hist_idx = hist_idx + 1;
        history(hist_idx, :) = x_current;
 
        if norm(x_current - x_start) > E
            continue;
        else
            break;
        end
    end
 
    history = history(1:hist_idx, :);
 
    x_min = x_current;
    f_min = f(x_min);
    t = toc;
end
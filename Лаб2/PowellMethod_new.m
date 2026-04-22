function [x_min, f_min, t, n, history] = PowellMethod_new(x0, E, f)
    tic;

    x_base = x0;
    n = 0;
    dim = length(x0);

    S = eye(dim);
    h = 1.0;
    %iter = 0;
    kmax = 1e5*dim;

    history = zeros(kmax + 1, dim);
    history(1, :) = x0(:)';
    hist_idx = 1;

    while n <= kmax
        n = n + 1;
        x_start = x_base;
        x_current = x_base;
        

        % --- Сброс базиса (закомментировать для работы без сброса) ---
        %iter = iter + 1;
        %if iter >= dim + 1
        %    S = eye(dim);
        %    iter = 1;
        %end
        % -------------------------------------------------------------

        
        f_vals = zeros(1, dim + 1);
        f_vals(1) = f(x_current);


        for i = 1:dim
            s = S(:, i)';
            func = @(lambda) f(x_current + lambda * s);
            lambda_i = GoldenSection(-h, h, E, func);
            x_current = x_current + lambda_i * s;

            f_vals(i + 1) = f(x_current);
        end

        improvements = f_vals(1:dim) - f_vals(2:dim+1);
        [~, k] = max(improvements);

        if norm(x_current - x_start) <= E
            break;
        else

            for i = k:dim-1
                S(:, i) = S(:, i + 1);
            end

            new_dir = (x_current - x_start)';
            S(:, dim) = new_dir / norm(new_dir);

            x_base = x_current;

            hist_idx = hist_idx + 1;
            history(hist_idx, :) = x_base;

        end
    end
    
    history = history(1:hist_idx, :);
    x_min = x_current;
    f_min = f(x_min);
    t = toc;
end
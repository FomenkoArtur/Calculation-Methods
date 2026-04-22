function [x_min, f_min, t, n, history] = RosenbrockMethod(x0, E, f)
    tic;
    x_base = x0;
    n = 0;
    dim = length(x0);
    kmax = 1e5*dim;
    
    S = eye(dim);
    
    lambda = zeros(dim, 1);
    
    h = 1.0;

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
    
    
        for i = 1:dim
            s = S(:, i)';
    
            func = @(lam) f(x_current + lam * s);
            
            lambda(i) = GoldenSection(-h, h, E, func);
           
            x_current = x_current + lambda(i) * s;
        end

        hist_idx = hist_idx + 1;
        history(hist_idx, :) = x_current;

        if norm(x_current - x_start) <= E
            break;
        else
            A = zeros(dim, dim);
    
            for i = 1:dim
                if abs(lambda(i)) < 1e-10
                    A(:, i) = S(:, i); 
                else
                    a_vec = zeros(dim, 1);
                    for j = i:dim
                        a_vec = a_vec + lambda(j) * S(:, j);
                    end
                    A(:, i) = a_vec;
                end
            end
            
            S_new = zeros(dim, dim);
            for i = 1:dim
                v = A(:, i)';
                
                for j = 1:i-1
                    proj = dot(v, S_new(:, j)');
                    v = v - proj * S_new(:, j)';
                end
                
                norm_v = norm(v);
                if norm_v > 1e-10
                    S_new(:, i) = v' / norm_v;
                else
                    S_new(:, i) = S(:, i); 
                end
            end
            
            S = S_new;
            x_base = x_current;
        end
    end

    
    history = history(1:hist_idx, :);
    x_min = x_current;
    f_min = f(x_min);
    t = toc;
end
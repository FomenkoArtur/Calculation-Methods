function [x_min, f_min, t, n, history] = NewtonMethod(x0, E, f)
    tic;
    n = 0;
    x = x0(:);
    dim = length(x0);
    kmax = 1e4 * dim;

    vars = symvar(f); 
    grad = gradient(f, vars);
    hess = hessian(f, vars);

    f_num = matlabFunction(f, 'Vars', {vars});
    grad_num = matlabFunction(grad, 'Vars', {vars});
    H_num = matlabFunction(hess, 'Vars', {vars});

    history = zeros(kmax + 1, dim);
    history(1, :) = x';
    hist_idx = 1;

    g = grad_num(x'); 
    g = g(:);
    
    H_val = H_num(x');
    H_inv = H_val^(-1);

    while n <= kmax
        n = n + 1;
        
        x = x - H_inv * g;
        
        hist_idx = hist_idx + 1;
        history(hist_idx, :) = x';

        g = grad_num(x');
        g = g(:);
        H_val = H_num(x');

        H_inv = H_val^(-1);

        grad_norm = norm(g);

        if grad_norm < E
            break;
        end
    end

    x_min = x;
    f_min = f_num(x_min');
    t = toc;
    history = history(1:hist_idx, :);
end
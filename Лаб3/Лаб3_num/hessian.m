function H = hessian(func, x, h)
    x = x(:);
    n = length(x);
    H = zeros(n, n);

    for i=1:n
        dx = zeros(n, 1);
        dx(i) = h;

        g_plus = gradient(func, x + dx, h);
        g_minus = gradient(func, x - dx, h);

        H(:, i) = (g_plus - g_minus) / (2 * h);
    end

    H = (H + H') / 2;
end

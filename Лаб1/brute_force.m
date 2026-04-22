function [x0, y0,t] = brute_force(a, b, lambda, E, f)
tic;
    % ---- Этап 1: грубый перебор на [a, b] с шагом lambda ----

    x0 = a;
    y0 = f(x0);
    x = a;
    while x <= b
        y = f(x);
        if y < y0
            x0 = x;
            y0 = y;
        end
        x = x + lambda;
    end

    % ---- Этап 2: уточнение около найденной точки ----

    while lambda > delta
        left = x0 - lambda;
        if left < a, left = a; end
        right = x0 + lambda;
        if right > b, right = b; end

        lambda = lambda * 0.5;       % уменьшаем шаг

        % Перебор на уточнённом интервале с новым шагом
        x_best = x0;
        y_best = y0;
        x = left;
        while x <= right
            y = f(x);
            if y < y_best
                x_best = x;
                y_best = y;
            end
            x = x + step;
        end

        % Обновляем результаты
        x0 = x_best;
        y0 = y_best;
    end
t=toc;    
end
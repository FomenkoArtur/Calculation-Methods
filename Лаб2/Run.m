clear; clc;
format long

a = 1;
b = 100;
f_2d  = @(x1, x2) Rozenbrock([x1, x2], a, b);
f_vec = @(x)      Rozenbrock(x,        a, b);

tests = {
    [-1.5, -1.5], 1e-4, 'x0=[-1.5,-1.5]';
    [ 0.0,  0.0], 1e-4, 'x0=[ 0.0, 0.0]';
    [ 2.0,  2.0], 1e-4, 'x0=[ 2.0, 2.0]';
};

fprintf('\n=====================================================================\n');
fprintf('  СРАВНЕНИЕ МЕТОДОВ МИНИМИЗАЦИИ ФНП  —  Функция Розенброка\n');
fprintf('  f(x) = (%.0f - x1)^2 + %.0f*(x2 - x1^2)^2\n', a, b);
fprintf('=====================================================================\n');

for i = 1:size(tests, 1)
    x0    = tests{i, 1};
    E     = tests{i, 2};
    label = tests{i, 3};

    [x_cd,  f_cd,  t_cd,  n_cd,  hist_cd]  = CoordinateDescent(x0, E, f_vec);
    [x_ps,  f_ps,  t_ps,  n_ps,  hist_ps]  = PatternSearch    (x0, E, f_vec);
    [x_hj,  f_hj,  t_hj,  n_hj,  hist_hj]  = HookeJeeves      (x0, E, f_vec);
    [x_hjn, f_hjn, t_hjn, n_hjn, hist_hjn] = HookeJeeves_new  (x0, E, f_vec);
    [x_rb,  f_rb,  t_rb,  n_rb,  hist_rb]  = RosenbrockMethod (x0, E, f_vec);
    [x_pw,  f_pw,  t_pw,  n_pw,  hist_pw]  = PowellMethod     (x0, E, f_vec);
    [x_pwn, f_pwn, t_pwn, n_pwn, hist_pwn] = PowellMethod_new (x0, E, f_vec);

    fprintf('\nТЕСТ: %s | Точность: %.0e\n', label, E);
    fprintf('---------------------------------------------------------------------\n');
    fprintf('%-20s | x1=%-10.6f x2=%-10.6f | F=%-12.6e | T=%-10.7f | N=%-5d\n', ...
        'CoordDescent',    x_cd(1),  x_cd(2),  f_cd,  t_cd,  n_cd);
    fprintf('%-20s | x1=%-10.6f x2=%-10.6f | F=%-12.6e | T=%-10.7f | N=%-5d\n', ...
        'PatternSearch',   x_ps(1),  x_ps(2),  f_ps,  t_ps,  n_ps);
    fprintf('%-20s | x1=%-10.6f x2=%-10.6f | F=%-12.6e | T=%-10.7f | N=%-5d\n', ...
        'HookeJeeves',     x_hj(1),  x_hj(2),  f_hj,  t_hj,  n_hj);
    fprintf('%-20s | x1=%-10.6f x2=%-10.6f | F=%-12.6e | T=%-10.7f | N=%-5d\n', ...
        'HookeJeeves_new', x_hjn(1), x_hjn(2), f_hjn, t_hjn, n_hjn);
    fprintf('%-20s | x1=%-10.6f x2=%-10.6f | F=%-12.6e | T=%-10.7f | N=%-5d\n', ...
        'Rosenbrock',      x_rb(1),  x_rb(2),  f_rb,  t_rb,  n_rb);
    fprintf('%-20s | x1=%-10.6f x2=%-10.6f | F=%-12.6e | T=%-10.7f | N=%-5d\n', ...
        'Powell',          x_pw(1),  x_pw(2),  f_pw,  t_pw,  n_pw);
    fprintf('%-20s | x1=%-10.6f x2=%-10.6f | F=%-12.6e | T=%-10.7f | N=%-5d\n', ...
        'Powell_new',      x_pwn(1), x_pwn(2), f_pwn, t_pwn, n_pwn);
    fprintf('---------------------------------------------------------------------\n');

    % ================================================================
    %  ВИЗУАЛИЗАЦИЯ СХОДИМОСТИ  (3D, вращаемый)
    % ================================================================
    figure('Name', sprintf('Сходимость: %s', label), 'NumberTitle', 'off');

    % --- Поверхность Розенброка в логарифмическом масштабе ---
    [X1g, X2g] = meshgrid(linspace(-2.2, 2.2, 300), linspace(-1.2, 3.2, 300));
    Zg = (a - X1g).^2 + b*(X2g - X1g.^2).^2;
    surf(X1g, X2g, log1p(Zg), 'EdgeColor', 'none', 'FaceAlpha', 0.85);
    colormap(gca, jet);    
    hold on;

    % --- Настройки методов: имя, история, цвет ---
    method_names = {
        'CoordDescent'
        'PatternSearch'
        'HookeJeeves'
        'HookeJeeves\_new'
        'Rosenbrock'
        'Powell'
        'Powell\_new'
    };
    histories = {
        hist_cd
        hist_ps
        hist_hj
        hist_hjn
        hist_rb
        hist_pw
        hist_pwn
    };
    colors = {
        [0.85, 0.15, 0.15]   % красный
        [0.10, 0.70, 0.20]   % зелёный
        [0.10, 0.35, 0.90]   % синий
        [0.90, 0.55, 0.00]   % оранжевый
        [0.70, 0.10, 0.80]   % фиолетовый
        [0.00, 0.75, 0.85]   % циановый
        [0.90, 0.20, 0.60]   % малиновый
    };

    % --- Рисуем траекторию каждого метода ---
    for k = 1:numel(histories)
        H = histories{k};
        if size(H, 1) < 2; continue; end

        Zh = arrayfun(@(r) log1p(f_vec(H(r, :))), 1:size(H, 1));

        plot3(H(:, 1), H(:, 2), Zh, ...
            '-o', ...
            'Color',           colors{k}, ...
            'LineWidth',       1.6, ...
            'MarkerSize',      3, ...
            'MarkerFaceColor', colors{k}, ...
            'DisplayName',     method_names{k});
    end

    % --- Стартовая точка и аналитический минимум ---
    z_x0  = log1p(f_vec(x0));
    z_opt = 0;
    plot3(x0(1), x0(2), z_x0, ...
        'ks', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'x^0');
    plot3(a, a^2, z_opt, ...
        'k*', 'MarkerSize', 12, 'LineWidth', 2, 'DisplayName', 'x^*');

    % --- Оформление ---
    xlabel('x_1',      'FontSize', 12);
    ylabel('x_2',      'FontSize', 12);
    zlabel('log(1+f)', 'FontSize', 12);
    title(sprintf('Траектории сходимости  |  %s  |  E=%.0e', label, E), ...
        'FontSize', 13);
    legend('Location', 'northeast', 'FontSize', 9);
    grid on;
    view(-40, 30);
    rotate3d on;
    % ================================================================
end

fprintf('\n=====================================================================\n');
fprintf('  Аналитический минимум: x* = (%.1f, %.1f),  f* = 0\n', a, a^2);
fprintf('=====================================================================\n');
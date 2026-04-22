clear; clc; format long

a_roz = 1; b_roz = 100;
a_ack = 20; b_ack = 0.2; c_ack = 2*pi;

syms z [2 1]

f_sym_roz  = (a_roz - z(1))^2 + b_roz*(z(2) - z(1)^2)^2;
f_sym_ack  = -a_ack*exp(-b_ack*sqrt((z(1)^2+z(2)^2)/2)) ...
             - exp((cos(c_ack*z(1))+cos(c_ack*z(2)))/2) + a_ack + exp(1);
f_sym_himm = (z(1)^2 + z(2) - 11)^2 + (z(1) + z(2)^2 - 7)^2;
f_sym_boot = (z(1) + 2*z(2) - 7)^2 + (2*z(1) + z(2) - 5)^2;
f_sym_rast = 10*2 + (z(1)^2 - 10*cos(2*pi*z(1))) + (z(2)^2 - 10*cos(2*pi*z(2)));

f_vec_roz  = @(x) (a_roz-x(1))^2 + b_roz*(x(2)-x(1)^2)^2;
f_vec_ack  = @(x) -a_ack*exp(-b_ack*sqrt((x(1)^2+x(2)^2)/2)) ...
                  - exp((cos(c_ack*x(1))+cos(c_ack*x(2)))/2) + a_ack + exp(1);
f_vec_himm = @(x) (x(1)^2+x(2)-11)^2 + (x(1)+x(2)^2-7)^2;
f_vec_boot = @(x) (x(1)+2*x(2)-7)^2 + (2*x(1)+x(2)-5)^2;
f_vec_rast = @(x) 20 + (x(1)^2-10*cos(2*pi*x(1))) + (x(2)^2-10*cos(2*pi*x(2)));

grids = {
    [-2.5,  2.5,  -1.5,  3.5,  300];
    [-5.0,  5.0,  -5.0,  5.0,  250];
    [-5.0,  5.0,  -5.0,  5.0,  250];
    [-5.0,  5.0,  -5.0,  5.0,  250];
    [-5.5,  5.5,  -5.5,  5.5,  250];
};


tests = {
    f_vec_roz,  f_sym_roz,  [-1.5,-1.5], 1e-4, ...
        'Rozenbrock',  'f=(1-x1)^2+100*(x2-x1^2)^2',           [1.0, 1.0], 0;
    f_vec_ack,  f_sym_ack,  [ 0.25, 0.25], 1e-4, ...
        'Ackley',      'f=-20*exp(...)-exp(...)+20+e',            [0.0, 0.0], 0;
    f_vec_himm, f_sym_himm, [ 0.0, 0.0], 1e-4, ...
        'Himmelblau',  'f=(x1^2+x2-11)^2+(x1+x2^2-7)^2',        [3.0, 2.0], 0;
    f_vec_boot, f_sym_boot, [ 0.0, 0.0], 1e-4, ...
        'Booth',       'f=(x1+2x2-7)^2+(2x1+x2-5)^2',           [1.0, 3.0], 0;
    f_vec_rast, f_sym_rast, [ 0.15, 0.15], 1e-4, ...
        'Rastrigin',   'f=20+x1^2-10cos(2pi*x1)+x2^2-10cos(2pi*x2)', [0.0, 0.0], 0;
};

method_names = {'SteepestDescent','ConjugateGradient','PolakRibiere'};
colors = {[0.85,0.10,0.10]; [0.05,0.65,0.10]; [0.10,0.30,0.95]};
markers = {'o','s','^'};


fprintf('\n');
fprintf('==========================================================================\n');
fprintf('        СРАВНЕНИЕ МЕТОДОВ МИНИМИЗАЦИИ ФНП (методы первого порядка)\n');
fprintf('==========================================================================\n');

for i = 1:size(tests,1)
    f_vec  = tests{i,1};
    f_sym  = tests{i,2};
    x0     = tests{i,3};
    E      = tests{i,4};
    fname  = tests{i,5};
    fdesc  = tests{i,6};
    x_opt  = tests{i,7};
    f_opt  = tests{i,8};
    g      = grids{i};

    [x_sd, f_sd, t_sd, n_sd, hist_sd] = SteepestDescentMethod (x0, E, f_sym);
    [x_cg, f_cg, t_cg, n_cg, hist_cg] = ConjugateGradientMethod(x0, E, f_sym);
    [x_pr, f_pr, t_pr, n_pr, hist_pr] = PolakRibiereMethod     (x0, E, f_sym);

    fprintf('\n');
    fprintf('--------------------------------------------------------------------------\n');
    fprintf('  ФУНКЦИЯ : %s\n', fname);
    fprintf('  ФОРМУЛА : %s\n', fdesc);
    fprintf('  Начальная точка : x0=[%.2f, %.2f] | Точность: %.0e\n', x0(1), x0(2), E);
    fprintf('--------------------------------------------------------------------------\n');
    fprintf('  %-22s | x1=%-12.6f x2=%-12.6f | F=%-14.6e | T=%-11.7f | N=%-6d\n', ...
        'SteepestDescent',   x_sd(1), x_sd(2), f_sd, t_sd, n_sd);
    fprintf('  %-22s | x1=%-12.6f x2=%-12.6f | F=%-14.6e | T=%-11.7f | N=%-6d\n', ...
        'ConjugateGradient', x_cg(1), x_cg(2), f_cg, t_cg, n_cg);
    fprintf('  %-22s | x1=%-12.6f x2=%-12.6f | F=%-14.6e | T=%-11.7f | N=%-6d\n', ...
        'PolakRibiere',      x_pr(1), x_pr(2), f_pr, t_pr, n_pr);
    fprintf('--------------------------------------------------------------------------\n');
    fprintf('  Аналитический минимум: x*=[%.4f, %.4f], f*=%.4f\n', x_opt(1), x_opt(2), f_opt);


    figure('Name', sprintf('Сходимость: %s', fname), ...
           'NumberTitle', 'off', 'Color', 'white');

    histories    = {hist_sd, hist_cg, hist_pr};


    hist_clipped = cell(3,1);
    for k = 1:3
        H = histories{k};
        if isempty(H) || size(H,1) < 1
            hist_clipped{k} = [];
            continue;
        end
        H(:,1) = max(g(1), min(g(2), H(:,1)));
        H(:,2) = max(g(3), min(g(4), H(:,2)));
        hist_clipped{k} = H;
    end

    x1lin = linspace(g(1), g(2), g(5));
    x2lin = linspace(g(3), g(4), g(5));
    [X1g, X2g] = meshgrid(x1lin, x2lin);
    Zg = zeros(size(X1g));
    for r = 1:size(X1g,1)
        for c = 1:size(X1g,2)
            Zg(r,c) = f_vec([X1g(r,c), X2g(r,c)]);
        end
    end


    surf(X1g, X2g, log1p(Zg), 'EdgeColor', 'none', 'FaceAlpha', 0.85);
    colormap(gca, jet);
    hold on;


    for k = 1:numel(histories)
        H = hist_clipped{k};
        if isempty(H) || size(H, 1) < 2; continue; end


        Zh = arrayfun(@(r) log1p(f_vec(H(r, :))), 1:size(H, 1));

        plot3(H(:, 1), H(:, 2), Zh, ...
            '-o', ...
            'Color',           colors{k}, ...
            'LineWidth',       1.6, ...
            'MarkerSize',      3, ... 
            'MarkerFaceColor', colors{k}, ...
            'DisplayName',     method_names{k});
    end

    
    z_x0 = log1p(f_vec(x0));
    plot3(x0(1), x0(2), z_x0, ...
        'ks', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', 'x^0');

    plot3(x_opt(1), x_opt(2), 0, ...
        'k*', 'MarkerSize', 12, 'LineWidth', 2, 'DisplayName', 'x^*');

    xlabel('x_1',      'FontSize', 12);
    ylabel('x_2',      'FontSize', 12);
    zlabel('log(1+f)', 'FontSize', 12);
    title(sprintf('Траектории сходимости | %s | x0=[%.1f,%.1f] | E=%.0e', ...
                  fname, x0(1), x0(2), E), ...
          'FontSize', 13);
    legend('Location', 'northeast', 'FontSize', 9);
    grid on;
    view(-40, 30);
    rotate3d on;
    
    xlim([g(1), g(2)]);
    ylim([g(3), g(4)]);
end

fprintf('\n');
fprintf('==========================================================================\n');
fprintf('  ВСЕ ТЕСТЫ ЗАВЕРШЕНЫ\n');
fprintf('==========================================================================\n');
function PlotNegativeWCases(viz_data, f, vars, plot_mode)
    if isempty(viz_data), return; end
    
    % Берём последние 3 случая (или меньше)
    N = min(3, length(viz_data));
    start_idx = max(1, length(viz_data) - N + 1);
    
    figure('Name', 'Polak-Ribiere: w < 0 cases', 'Position', [100,100,1600,500]);
    
    % Создаём сетку для контурного графика функции
    x1_range = linspace(-2, 2, 100);
    x2_range = linspace(-1, 3, 100);
    [X1, X2] = meshgrid(x1_range, x2_range);
    
    % Создаем функцию от двух отдельных аргументов
    f_handle = matlabFunction(f, 'Vars', {vars(1), vars(2)});
    Z = arrayfun(f_handle, X1, X2);
    
    len = 0.6;  % Чуть увеличили длину для наглядности
    
    for i = 1:N
        idx = start_idx + i - 1;
        data = viz_data(idx);
        
        subplot(1, N, i);
        contour(X1, X2, Z, 50, 'LineWidth', 0.5); hold on;
        colormap(parula); colorbar;
        
        % 1. Текущий градиент (красный)
        g_norm = norm(data.grad_x);
        if g_norm > 0
            quiver(data.x(1), data.x(2), ...
                   len*data.grad_x(1)/g_norm, len*data.grad_x(2)/g_norm, ...
                   0, 'r', 'LineWidth', 3, 'MaxHeadSize', 0.5, 'AutoScale', 'off');
        end
        
        % 2. Предыдущий градиент (синий)
        gp_norm = norm(data.grad_prev);
        if gp_norm > 0
            quiver(data.x(1), data.x(2), ...
                   len*data.grad_prev(1)/gp_norm, len*data.grad_prev(2)/gp_norm, ...
                   0, 'b', 'LineWidth', 3, 'MaxHeadSize', 0.5, 'AutoScale', 'off');
        end
        
        % 3. Предыдущее направление поиска (зелёный)
        % ИСПРАВЛЕНИЕ: пунктир и подпись
        sp_norm = norm(data.s_prev);
        if sp_norm > 0
            quiver(data.x(1), data.x(2), ...
                   len*data.s_prev(1)/sp_norm, len*data.s_prev(2)/sp_norm, ...
                   0, '--g', 'LineWidth', 2.5, 'MaxHeadSize', 0.5, 'AutoScale', 'off');
            
            % Подпись рядом с концом стрелки
            text(data.x(1) + 1.1*len*data.s_prev(1)/sp_norm, ...
                 data.x(2) + 1.1*len*data.s_prev(2)/sp_norm, ...
                 's_{k-1}', 'Color', 'g', 'FontWeight', 'bold', 'FontSize', 14);
        end
        
        % Точка (рисуем последней, чтобы была поверх векторов)
        plot(data.x(1), data.x(2), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'w', 'LineWidth', 2);
        
        title(sprintf('Iter %d: w = %.3e', data.iter, data.w_raw), ...
              'FontSize', 14, 'FontWeight', 'bold');
        
        if i == 1
             legend({'Уровни', '\nabla f_k', '\nabla f_{k-1}', 's_{k-1} (зелёный)', 'Точка'}, ...
                   'Location', 'northwest', 'FontSize', 10);
        end
        
        xlabel('x_1'); ylabel('x_2');
        axis equal; grid on;
        % Расширяем границы, чтобы стрелки не обрезались
        xlim([-2.5, 2.5]); 
        ylim([-1.5, 3.5]); 
    end
    
    sgtitle('Геометрическая интерпретация w < 0 в методе Полака-Рибьера', 'FontSize', 16, 'FontWeight', 'bold');
end
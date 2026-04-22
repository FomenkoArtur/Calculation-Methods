function [x_min, f_min, t, n, w_log, viz_data] = PolakRibiereViz(x0, E, f, plot_flag)
    % plot_flag: 'none' | 'final' | 'all' — уровень визуализации
    
    [x_min, f_min, t, n] = PolakRibiereMethod(x0, E, f);
    
    vars     = symvar(f);
    grad_sym = gradient(f, vars);
    f_num    = matlabFunction(f,        'Vars', {vars});
    grad_num = matlabFunction(grad_sym, 'Vars', {vars});
    
    x         = x0(:);
    h         = 1.0;
    kmax      = 1e5 * length(x0);
    
    grad_prev = grad_num(x');
    s_prev    = -grad_prev;
    func      = @(alpha) f_num((x + alpha * s_prev)');
    alpha     = GoldenSection(0, h, E, func);
    x         = x + alpha * s_prev;
    
    w_log    = [];      % [iter, grad_x, grad_prev, s_prev, w_raw]
    viz_data = struct('x',{}, 'grad_x',{}, 'grad_prev',{}, 's_prev',{}, 'w_raw',{}, 'iter',{});
    
    for k = 1:kmax
        grad_x      = grad_num(x');
        numerator   = grad_x' * (grad_x - grad_prev);
        denominator = grad_prev' * s_prev;
        
        w_raw = 0;
        if abs(denominator) >= 1e-15
            w_raw = numerator / denominator;
        end
        
        % 🔍 Сохраняем данные для визуализации, если w < 0
        if w_raw < 0
            w_log(end+1, :) = [k, grad_x', grad_prev', s_prev', w_raw]; %#ok
            
            if ~strcmp(plot_flag, 'none')
                viz_data(end+1) = struct(...
                    'x', x, ...
                    'grad_x', grad_x, ...
                    'grad_prev', grad_prev, ...
                    's_prev', s_prev, ...
                    'w_raw', w_raw, ...
                    'iter', k);
            end
        end
        
        w      = max(0, w_raw);  % Restart при w < 0
        s      = -grad_x + w * s_prev;
        func   = @(alpha) f_num((x + alpha * s)');
        alpha  = GoldenSection(0, h, E, func);
        x      = x + alpha * s;
        grad_prev = grad_x;
        s_prev    = s;
        
        if norm(grad_num(x')) < E, break; end
    end
    
    % 🎨 Построение финальной визуализации
    if ~isempty(viz_data) && ~strcmp(plot_flag, 'none')
        PlotNegativeWCases(viz_data, f, vars, plot_flag);
    end
end
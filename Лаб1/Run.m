clear; clc;

format long
tests = {
    -50, 50, 0.01, @f1, 'f1';
    0, 8, 0.01, @f2, 'f2';
    1.5, 2, 0.01, @f3, 'f3';
};

n = 15;
lambda = 0.5;


fprintf('\n======================================================================\n');
fprintf('  СРАВНЕНИЕ МЕТОДОВ ОПТИМИЗАЦИИ\n');
fprintf('======================================================================\n');

for i = 1:size(tests, 1)
    [a, b, E, f, name] = deal(tests{i, 1}, tests{i, 2}, tests{i, 3}, ...
                               tests{i, 4}, tests{i, 5});

    [xb, fb, tb] = brute_force(a, b, lambda, E, f);
    [xd, fd, td, nd] = dihtom(a, b, 0.01, E, f);
    [xg, fg, tg, ng] = GoldenSection(a, b, E, f);
    [xf, ff, tf] = Fibonacci(a, b, n, f);
    
    fprintf('ТЕСТ: %s | Интервал: [%.1f, %.1f] | Точность: %.3f\n', name, a, b, E);
    fprintf('----------------------------------------------------------------------\n');
    fprintf('%-12s | X = %-10.6f | F = %-10.6f | T = %-8.6f | N = %-5s\n', ...
        'BruteForce', xb, fb, tb, '-');
    fprintf('%-12s | X = %-10.6f | F = %-10.6f | T = %-8.6f | N = %-5d\n', ...
        'Dichotomy',  xd, fd, td, nd);
    fprintf('%-12s | X = %-10.6f | F = %-10.6f | T = %-8.6f | N = %-5d\n', ...
        'Golden',     xg, fg, tg, ng);
    fprintf('%-12s | X = %-10.6f | F = %-10.6f | T = %-8.6f | N = %-5d\n', ...
        'Fibonacci',  xf, ff, tf, n);
    fprintf('----------------------------------------------------------------------\n\n');
end

fprintf('======================================================================\n');

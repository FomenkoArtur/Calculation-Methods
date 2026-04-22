function [y] = Rozenbrock(x, a, b)
    y = (a - x(1))^2 + b * (x(2) - x(1)^2)^2;
end
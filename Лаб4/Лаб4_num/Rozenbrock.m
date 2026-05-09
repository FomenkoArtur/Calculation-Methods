function [y] = Rozenbrock(z, a, b)
    y = (a - z(1))^2 + b * (z(2) - z(1)^2)^2;
end
function [x_m, f_m, t] = Fibonacci(a, b, n, f)
tic;
c = a + (fibonacci(n)/fibonacci(n+2))*(b-a);
d = a + (fibonacci(n+1)/fibonacci(n+2))*(b-a);
fc = f(c);
fd = f(d);
for i = (1:n-1)
    if (fc<=fd)
        b=d;
        d=c;
        c=a+b-d;
        fd = fc;
        fc = f(c);
    else
        a=c;
        c=d;
        d=a+b-c;
        fc = fd;
        fd = f(d);
    end
end

x_m=c;
f_m=f(x_m);

t=toc;
end
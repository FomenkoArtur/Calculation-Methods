function [x_m, f_m, t, n] = dihtom(a, b, k, E, f)
tic;
delta = E*k;
n=0;
while abs(b-a)>E
    n=n+1;
    c = (a+b-delta)/2;
    d = (a+b+delta)/2;
    if (f(c)<=f(d))
        b=d;
    else
        a=c;
    end
end

x_m=(b+a)/2;
f_m=f(x_m);
t=toc;
end
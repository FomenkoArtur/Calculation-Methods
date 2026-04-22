function val = f_single_coord(f, x, i, v)
    x(i) = v;
    val = f(x);
end
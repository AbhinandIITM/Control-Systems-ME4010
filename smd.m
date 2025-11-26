function dx = smd(x, m, k, c, u)
   

    dx(1,1) = x(2);
    dx(2,1) = -k/m * x(1) - c/m * x(2) + u/m;
end

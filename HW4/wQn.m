function wQn = wQn(xk)
    xt = xk';
    m = 2;
    n = 10^m;
    d = ones(n,1);
    A = spdiags([d -2*d d],[-1 0 1],n,n);
    wQn = sqrt((xt)*A*xk);
end
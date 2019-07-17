function cfsd = sd(A,xk,b)
rk = A*xk - b;
pk = -rk;
tol = 1e-9;
abar = 1;
hh = 0;
while norm(rk) > tol
    ak = abar;
    xk = xk + ak*pk;
    hh = hh+1;
    musd = eig(A);
    musd = sort(musd); %sorting from smallest to biggest
    if hh<= n
       kappasd = musd(hh)/musd(1); %condition number for A
    end
    cfsd(hh) = abs((sqrt(kappasd)-1)/(sqrt(kappasd)+1)); %convergence factor for conjugate gradient
end
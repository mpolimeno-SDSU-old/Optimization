m = 0;
n = 10^m;
dim = 2;
d = ones(n^dim,1);
A = spdiags([d d -4*d d d],[-n -1 0 1 n],n^dim,n^dim);
xk = zeros(size(d));
b = ones(size(d));
rk = A*xk - b;
pk = -rk;
kk = 0;
tol = 1e-9;
tic
while norm(rk) > tol
    rt = rk';
    pt = pk';
    Apk = A*pk; %vector
    pAp = pt*Apk; %scalar
    rtr = rt*rk; 
    ak = rtr/pAp; %scalar
    xkp1 = xk + ak*pk;
    rkp1 = rk + ak*Apk;
    rkpt = rkp1';
    num = rkpt*rkp1; %scalar
    bkp1 = num/rtr; %scalar
    pkp1 = -rkp1 + bkp1*pk;
    xk = xkp1;
    rk = rkp1;
    pk = pkp1;
    kk = kk+1;
    nrk(kk) = norm(rk);
end
toc
endt = toc;

mu = eig(A);
kappa = mu(n)/mu(1); %condition number for A

B = nnz(A); %number of nonzero elements in A
nvec = 0:length(nrk)-1;

figure(1)
semilogy(nvec,nrk,'b-','linewidth',1.5);
legend('L_{2} norm')
xlabel('Iterations')
ylabel('L_{2} norm')
grid on
title(['Residual 2-norm vs. Number of Iterations for n = ' num2str(n)])
n = 20;
A = zeros(n,n);
for ii=1:n
    for jj=1:n
        A(ii,jj) = 1/(ii+jj-1);
    end
end
xk = zeros(n,1);
b = ones(n,1);
rk = A*xk - b;
pk = -rk;
kk = 0;
tol = 1e-6;

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
mu = sort(mu); %sorting from smallest to biggest
kappa = mu(n)/mu(1); %condition number for A

nvec = 0:length(nrk)-1;
muvec = 0:n-1;

figure(1)
semilogy(nvec,nrk,'b-','linewidth',1.5);
legend('L_{2} norm')
xlabel('Iterations')
ylabel('L_{2} norm')
grid on
title(['Residual 2-norm vs. Number of Iterations for n = ' num2str(n)])

mu = log(abs(mu));
figure(2)
plot(muvec,mu,'k*','linewidth',1.5)
title(['Spread of Eigenvalues for n = ' num2str(n)])
xlabel('Iterations')
ylabel('log_{10}(|\lambda|)')
grid on
m = 3;
n = 10^m;
d = ones(n,1);
A = spdiags([d -2*d d],[-1 0 1],n,n);
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
mu = sort(mu); %sorting from smallest to biggest
kappa = mu(n)/mu(1); %condition number for A
B = nnz(A); %number of nonzero elements in A
nvec = 0:length(nrk)-1;

C = A^-1;
D = nnz(C);

figure(1)
semilogy(nvec,nrk,'b-','linewidth',1.5);
legend('L_{2} norm')
xlabel('Iterations')
ylabel('L_{2} norm')
grid on
title(['Residual 2-norm vs. Number of Iterations for n = ' num2str(n)])

figure(2)
spy(A,'b')
hold on
spy(A==-2*d,'r');
hold off
xlabel(['NNZ = ' num2str(B)])
legend('Sub/Super-diagonal elements','Diagonal Elements')
title('Pattern of non-zero elements for Matrix A')

figure(3)
spy(C,'bo')
xlabel(['NNZ = ' num2str(D)])
title('Pattern of non-zero elements for A^{-1}')
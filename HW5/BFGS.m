%code for the BFGS algorithm
x0 = rosenbrock_2Nd([],-1);
xk = x0;
Hk = eye(size(rosenbrock_2Nd(xk,2)));
Bk = eye(size(rosenbrock_2Nd(xk,2)));
abar = 1;
al = 0;
ah = 1;
amax = 5;
c1 = 1e-4;
c2 = 0.9;
rho = .5;
tol = 1e-11;
xstar = ones(18,1);
QNconv = [];
rkvec = [];
kk = 0;
Tin = [];

tic
while norm(rosenbrock_2Nd(xk,1)) > tol
    pk = -Hk*rosenbrock_2Nd(xk,1);
    pt = pk';
    tic
    a = wolfe_strong2(xk,pk,al,ah,amax,c1,c2);
    tin = toc;
    Tin = [Tin;tin];
    xkp1 = xk+a*pk;
    sk = xkp1 - xk;
    yk = rosenbrock_2Nd(xkp1,1)-rosenbrock_2Nd(xk,1);
    yt = yk';
    st = sk';
    rhok = 1/((yt)*sk);
    Hk = (eye(size(Hk))-rhok*sk*yt)*Hk*(eye(size(Hk))-rhok*yk*st) + rhok*sk*st;
    Bk = Bk - (Bk*sk*st*Bk)/(st*Bk*sk) + (yk*yt)/(yt*sk);
    Qconv = norm((Bk-Hk)*pk)/norm(pk);
    QNconv = [QNconv;Qconv];
    xk = xkp1;
    kk = kk + 1;
    rk = norm(xk-xstar);
    rkvec = [rkvec;rk];
end
toc

gradk = rosenbrock_2Nd(xk,1);
%mu = eig(Hk);
figure(1)
kvec = 1:kk;
plot(kvec,QNconv,'k-','linewidth',2)
title('Quasi-Newton Criteria for Convergence')
xlabel('Iterations','interpreter','latex','fontsize',15)
ylabel('$\frac{||B_{k}-\nabla^{2}f({\bf x_{k}}){\bf p_{k}}||}{||{\bf p_{k}}||}$','interpreter','latex','fontsize',14)
grid on

% % Get fitted values
% coeffs = polyfit(kvec', QNconv, 1);
% fittedY = polyval(coeffs, kvec);

figure(2)
kvec = 1:kk;
semilogy(kvec,QNconv,'k-','linewidth',2)
% hold on
% % Plot the fitted line
% plot(kvec, fittedY, '--', 'LineWidth', 2);
% hold off
title('Quasi-Newton Criteria for Convergence')
xlabel('Iterations','interpreter','latex','fontsize',15)
ylabel('$\frac{||B_{k}-\nabla^{2}f({\bf x_{k}}){\bf p_{k}}||}{||{\bf p_{k}}||}$','interpreter','latex','fontsize',14)
grid on

figure(3)
kvec = 1:kk;
semilogy(kvec,rkvec,'k-','linewidth',2)
title('Residual 2-Norm')
xlabel('Iterations','interpreter','latex','fontsize',15)
ylabel('$||{\bf x_{k}} - {\bf x_{\ast}}||_{2}$','interpreter','latex','fontsize',15)
ylim([10^-10 10^2])
grid on
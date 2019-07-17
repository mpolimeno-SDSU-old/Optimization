%trust region
clear all
clc
nn = 100000;
mu = (min(diag(zeros(2))));
x = [0,-1];
Dk = 2;
tol = 1e-8;
xrange = linspace(-Dk,Dk,101);
tic
t = false;
for ii = 1:nn
    nabla = grad_f(x);
    if ii==1
        A = hess_f(x);
        mu(ii) = min(abs(diag(A)));
        if A(ii,ii) < 0
            mu(ii) = min(abs(diag(A))) + 0.1;
        end
    end
    A = hess_f(x) + mu(ii)*eye(length(hess_f(x)));
    if ii==2 && (A(ii,ii) < 0)
        mu(ii) = min(abs(diag(A))) + 0.1;
    end
    if mu(ii) < 0
       mu(ii) = min(abs(diag(A))) + 0.1;
    end
    L = cholesky(x,A);
    LT = (L)';
    LI = (L)^-1;
    q = ((LT)^-1)*(-nabla)';
    p = (LI)*q;
    mu(ii+1) = mu(ii) + (((norm(p))/(norm(q)))^2 )*((norm(p))-Dk)/Dk;
    if abs(mu(ii+1)-mu(ii)) < tol
        t = true;
        break
    end
    mu(ii) = mu(ii+1);
    xkp1 = x + p';
end
toc
mu(ii)
xcontour = linspace(-4,4,101);
[X,Y] = meshgrid(xcontour,xcontour);
pt = p';
m = zeros(length(xcontour),length(xcontour));

for jj=1:length(xcontour)
    x2 = xcontour;
    for kk=1:length(xcontour)
        f = h([x2(jj),x2(kk)]);
        g = grad_f([x2(jj),x2(kk)]);
        Bk = hess_f([x2(jj),x2(kk)])+mu(kk)*eye(length(hess_f([x2(jj),x2(kk)])));
        m(jj,kk) = f + pt*(g')+.5*pt*(Bk)*p;
    end
end

y1 = sqrt(Dk^2-(xrange-x(1)).^2)+x(2);
y2 = -sqrt(Dk^2-(xrange-x(1)).^2)+x(2);
figure(1)
mp = contour(X,Y,m,50);
hold on
circ1 = plot(xrange,y1,'r--','linewidth',1.5);
hold on
circ2 = plot(xrange,y2,'r--','linewidth',1.5);
hold on
x0 = plot(x(1),x(2),'b*','linewidth',2);
hold on
x1 = plot(xkp1(1),xkp1(2),'k*','linewidth',2);
hold on
p = plot([x(1) xkp1(1)],[x(2) xkp1(2)],'g-','linewidth',1.5);
hold off
set(gca, 'Xlim',[-4,4])
set(gca,'Xtick',(-4:0.5:4))
set(gca, 'Ylim',[-4,4])
set(gca,'Ytick',(-4:0.5:4))
title('Contour Plot of m({\bf p}) for {\bf x} = [0,-1]','interpreter','latex')
legend([circ1,p, x0, x1],{'Trust Region','{\bf p}','x_{0}','x_{k+1}'})
colorbar
axis equal

%cholesky
function L = cholesky(~,A)
n = length(A);
L = zeros(n,n);
for i = 1:n
    L(i,i) = sqrt(A(i,i));
    if A(i,i) < 0
       break
    end
    for j = i+1:n
        L(j,i) = A(j,i)/(L(i,i));
        for k = i+1:j
            A(j,k) = A(j,k) - L(j,i)*L(k,i);
        end
    end
end
end

function f = h(x)
f = 10*(x(2)-x(1)^2)^2+(1-x(1))^2;
end

function nabla = grad_f(x)
nabla = [40*x(1)*(x(1)^2-x(2))+2*(x(1)-1), 20*(x(2)-x(1)^2)];
end

function A = hess_f(x)
A = [120*x(1)^2-40*x(2)+2, -40*x(1);...
        -40*x(1), 20];
end
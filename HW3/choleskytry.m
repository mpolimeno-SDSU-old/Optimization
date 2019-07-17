%function L = cholesky(A)
x = [0,-1];
A = hess_f(x);
n = length(A);
L = zeros(n,n);
for i = 1:n
    %n = length(A);
    L(i,i) = sqrt(A(i,i));
    if A(i,i) < 0
       L(i,i) = min(abs(diag(A)));
    else
        for j = i+1:n  %i+1:n, but i = 1
            L(j,i) = A(j,i)/(L(i,i));
            for k = i+1:j  %i+1:j, but j = 2
                A(j,k) = A(j,k) - L(j,i)*L(k,i);
            end
        end
    end
end
%end

function A = hess_f(x)
A = [120*x(1)^2-40*x(2)+2, -40*x(1);...
        -40*x(1), 20];
end
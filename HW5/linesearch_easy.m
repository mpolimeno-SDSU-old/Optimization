%function x = linesearch_easy(x)
x = rosenbrock_2Nd([],-1);
rho = .5;
abar = 1;
tol = 1e-8;
c1 = 1e-4;
ii = 0;

tic
while norm(rosenbrock_2Nd(x,1)) > tol
   twhile = tic;
   pk = dir_newt2(x);
   a = abar;
   pt = pk';
   while rosenbrock_2Nd(x+a*pk,0) > rosenbrock_2Nd(x,0) + c1*a*dot(pt,rosenbrock_2Nd(x,1))
       a = rho*a;
   end
   x = x+a*pk;
   ii = ii+1;
end
toc
disp(ii)
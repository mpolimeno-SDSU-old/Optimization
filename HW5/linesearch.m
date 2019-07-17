function x = linesearch(x,al,ah,amax,tol)
    alpha_bar = 1;
    c1 = 1e-4;
    c2 = 0.9;
    %Hk = eye(size(rosenbrock_2Nd(x,2)));
    pk = dir_newt2(x);%-Hk*rosenbrock_2Nd(x,1);
    ii = 0;
    tic
    while (norm(rosenbrock_2Nd(x,1))) > tol
        a = wolfe_strong2(x,pk,al,ah,amax,c1,c2);
        xnew = x + a*pk;
        x = xnew;
        ah = alpha_bar;
        ii = ii + 1;
    end
    toc
    disp(ii)
end
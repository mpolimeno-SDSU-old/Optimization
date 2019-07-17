function anew = wolfe_strong2(x,pk,al,ah,amax,c1,c2) %strong wolfe conditions
    aii = [al,ah];  
    ii = 1;
    while true
        if (rosenbrock_2Nd(x+aii(2)*pk,0) > rosenbrock_2Nd(x,0) + c1*aii(2)*rosenbrock_2Nd(x,0)) || ((rosenbrock_2Nd(x+aii(2)*pk,0) >= rosenbrock_2Nd(x+aii(1)*pk,0))) && ii > 1
            anew = zoom2(x,pk,aii(1),aii(2),c1,c2);
            break
        end
        if abs(phi_prime(x,pk,aii(2))) <= -c2*phi_prime(x,pk,0)
            anew = aii(2);
            break
        end
        if phi_prime(x,pk,aii(2)) >= 0
            anew = zoom2(x,pk,aii(2),aii(1),c1,c2);
            break
        end
        if (abs(aii(2)-aii(1)) < 1e-14) || (abs(aii(2)) < 1e-14)
            anew = aii(2);
            break
        end
        aii(1) = aii(2);
        aii(2) = (amax+aii(2))/2;
        ii = ii + 1;
    end
end
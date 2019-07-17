function astar = zoom2(x,pk,al,ah,c1,c2)
    go = true;
    while go
        ajj = abs(ah-al)/2;
        if (rosenbrock_2Nd(x+ajj*pk,0) > rosenbrock_2Nd(x,0)+c1*ajj*phi_prime(x,pk,0)) || (rosenbrock_2Nd(x+ajj*pk,0) >= rosenbrock_2Nd(x+al*pk,0))
            ah = ajj;
        else
            if abs(phi_prime(x,pk,ajj)) <= -c2*phi_prime(x,pk,0)
                astar = ajj;
                go = false;
            end
            if (phi_prime(x,pk,ajj)*(ah-al)) >= 0
                ah = al;
            al = ajj;
            end
        end
    end
end

% function akp1 = interp(x,al,ah)  %return result of the interpolation
%     d1 = phi_prime(x,al) + phi_prime(x,ah) - 3*(rosenbrock_2Nd(x+al*dir_newt2(x),0)-rosenbrock_2Nd(x+ah*dir_newt2(x),0))/(al-ah);
%     d2 = sign(ah-al)*(sqrt(d1.^2-phi_prime(x,al).*phi_prime(x,ah)));
%     akp1 = ah - (ah-al)*((phi_prime(x,ah)+d2-d1)./(phi_prime(x,ah)-phi_prime(x,al)+2*d2));
% end

% function phip = phi_prime(x,alpha)
%     phip = dot(rosenbrock_2Nd(x+alpha*dir_newt2(x),1),(dir_newt2(x)));
% end

% function pk = dir_newt2(x)
%     pk = -(rosenbrock_2Nd(x,2)\rosenbrock_2Nd(x,1));
% end
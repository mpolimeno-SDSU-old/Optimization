function akp1 = interp2(x,pk,al,ah)  %return result of the interpolation
    d1 = phi_prime(x,al) + phi_prime(x,ah) - 3*(rosenbrock_2Nd(x+al*pk,0)-rosenbrock_2Nd(x+ah*pk,0))/(al-ah);
    d2 = sign(ah-al)*(sqrt(d1.^2-phi_prime(x,al).*phi_prime(x,ah)));
    akp1 = ah - (ah-al)*((phi_prime(x,ah)+d2-d1)./(phi_prime(x,ah)-phi_prime(x,al)+2*d2));
end
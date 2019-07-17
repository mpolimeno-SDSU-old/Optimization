function phip = phi_prime(x,pk,alpha)
    phip = dot(rosenbrock_2Nd(x+alpha*pk,1),pk);
end
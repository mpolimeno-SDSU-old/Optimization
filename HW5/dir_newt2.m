function pk = dir_newt2(x)
      pk = -(rosenbrock_2Nd(x,2)\rosenbrock_2Nd(x,1));
end
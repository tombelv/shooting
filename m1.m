function m1_sym = m1(in1,sigma)

    m1_sym = f(in1) + sigma*norm(g(in1),1);
end
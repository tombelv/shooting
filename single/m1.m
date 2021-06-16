function m1_sym = m1(x,sigma)

    m1_sym = f(x) + sigma*norm(g(x),1);
    %m1_sym = f(x) + sigma*norm(g(x),1) + sigma*norm(max(h(x), 0),1);
end
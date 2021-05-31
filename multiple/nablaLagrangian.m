function nablaLagrangian_sym = nablaLagrangian(in1,in2,in3)
    [~,nablag_] = g(in1);
    nablaLagrangian_sym = nablaf(in1) + nablag_*in2 + nablah(in1)*in3;
end
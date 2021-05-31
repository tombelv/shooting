function nablaLagrangian_sym = nablaLagrangian(in1,in2,in3)
    nablaLagrangian_sym = nablaf(in1) + nablag(in1)*in2 + nablah(in1)*in3;
end
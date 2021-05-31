function nablag_sym = nablag(in1)
    global state_coeff s_init ns
    [~, nablaR] = single_shooting(in1, s_init);
    nablag_sym = (sqrt(state_coeff) \ nablaR(:,end-ns+1:end).').';
end

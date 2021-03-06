function g_sym = g(in1)
    global state_coeff s_init ns
    [R, ~] = single_shooting(in1, s_init);
    g_sym = sqrt(state_coeff) \ R(end-ns+1:end);
end

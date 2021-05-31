function B_sym = B(in1, ~, ~)
    global s_init
    [~, nablaR] = single_shooting(in1, s_init);
    B_sym = nablaR*nablaR.'; 
end
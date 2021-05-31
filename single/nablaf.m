function nablaf_sym = nablaf(input)
    global s_init
    [R, nablaR] = single_shooting(input, s_init);
    nablaf_sym = nablaR*R;
end
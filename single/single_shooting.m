function [R, nablaR] = single_shooting(u_guess, s0)
    global n_step ns t_init nu n_int int_step state_coeff input_coeff d_step s_fin s_init
    state_coeff_sqrt = sqrt(state_coeff);
    input_coeff_sqrt = sqrt(input_coeff);
    
    R = zeros((n_step+1)*ns + n_step*nu, 1);
    nablaR = zeros(n_step*nu, (n_step+1)*ns + n_step*nu, 1);
    sk = s0;
    R(1:ns) = state_coeff_sqrt * (sk-s_init);
    dxdw = zeros(ns, n_step * nu);
    nablaR(:,1:ns) = dxdw.';
    %increment for linearly growing reference for theta
    delta_theta = (s_fin(2)-s_init(2))/n_step;
    theta_ref = s_init(2);
    
    % single shooting
    for k = 0:(n_step - 1)
        u_init = u_guess(k*nu+1);
        tk = t_init + d_step*k;
        [s,A,B] = expl_rk4(tk,sk,u_init,eye(ns),zeros(ns,nu),int_step,n_int);
        dxdw = A * dxdw + [zeros(ns, k * nu) B zeros(ns, (n_step-k-1) * nu)];
        nablaR((k*nu+1):(k*nu+nu),(k*(ns+nu)+ns+1):(k*(ns+nu)+ns+nu)) = input_coeff_sqrt;
        nablaR(:,((k+1)*(ns+nu)+1):((k+1)*(ns+nu)+ns)) = (state_coeff_sqrt * dxdw).';
        % update s_init and u_init
        sk = s;
        R((k*(ns+nu)+ns+1):(k*(ns+nu)+ns+nu)) = input_coeff_sqrt * u_init;
        theta_ref = theta_ref + delta_theta;
        s_ref = [s_fin(1); theta_ref; s_fin(3:end)];
        R(((k+1)*(ns+nu)+1):((k+1)*(ns+nu)+ns)) = state_coeff_sqrt * (sk-s_ref);
    end
end
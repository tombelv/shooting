commandwindow
clear all

addpath('..')

global n_step d_step t_init s_init s_fin ns nu n_int int_step state_coeff input_coeff
s_init = [0;0;0;0];
s_fin = [0;pi;0;0];
nu = 1;
% problem parameters
M = 1;
m = 1;
l = 1;
gravity = 9.81;
u_max = 20;
state_coeff = [100 0 0 0;
                 0 1 0 0;
               0 0 0.1 0;
                 0 0 0 0.1];
state_coeff = 0.00001*state_coeff;
input_coeff = 0.001;
% optimization horizon
t_init = 0;
t_fin = 2;
n_step = 40;
d_step = (t_fin-t_init)/n_step;

ns = length(s_init);
syms t 
u = sym('u', [nu 1]);
s = sym('s', [ns 1]);

% integration horizon
% each optimization interval do n_int integration steps of duration
% int_step
n_int = 20;
int_step = d_step / n_int;
times = linspace(t_init,t_fin,n_step+1);

iters = 1;

% INPUT
tol = 1e-4;
u_init = 0;
w_init = s_init;
for i=1:n_step/2
    w_init = [w_init; u_init; s_init];
end
for i=(n_step/2+1):n_step
    w_init = [w_init; u_init; s_fin];
end
nw = length(w_init);
sigma_coeff = 2;
sigma_init = 0;
damping_coeff = 0.5;

linesearch = 'MERIT';

w = sym('w', [nw; 1]);
sigma = sym('sigma');

F_sym = [s(3); 
    s(4); 
    (m*l*sin(s(2))*s(4)^2 + m*gravity*cos(s(2))*sin(s(2)) + u) / (M+m-m*(cos(s(2)))^2);
    -((m*l*cos(s(2))*sin(s(2))*s(4)^2 + (M+m)*gravity*sin(s(2)) + u*cos(s(2))) / (l*(M+m-m*(cos(s(2)))^2)))];
dFds_sym = jacobian(F_sym, s);
dFdu_sym = jacobian(F_sym, u);

matlabFunction(F_sym, 'vars', {t,s,u}, 'file', 'dynamics');
matlabFunction(dFds_sym, 'vars', {t,s,u}, 'file', 'dFds');
matlabFunction(dFdu_sym, 'vars', {t,s,u}, 'file', 'dFdu');

% only for expl rk4
k1_sym = dynamics(t,s,u);
k2_sym = dynamics(t + 1/2 * int_step, s + 1/2 * int_step * k1_sym, u);
k3_sym = dynamics(t + 1/2 * int_step, s + 1/2 * int_step * k2_sym, u);
k4_sym = dynamics(t + int_step, s + int_step * k3_sym, u);
s_sym = s + 1/6 * int_step * (k1_sym +2*k2_sym + 2*k3_sym + k4_sym);
dsds_sym = jacobian(s_sym,s);
dsdu_sym = jacobian(s_sym,u);

matlabFunction(s_sym, 'vars', {t,s,u}, 'file', 's');
matlabFunction(dsds_sym, 'vars', {t,s,u}, 'file', 'dsds');
matlabFunction(dsdu_sym, 'vars', {t,s,u}, 'file', 'dsdu');

% set the inequality constraints
w_input = w(ns+1:ns+nu);
for i = 1:n_step-1
    w_input = [w_input; w(i*(ns+nu)+ns+1:i*(ns+nu)+ns+nu)];
end
h_sym = [eye(nu*n_step); -eye(nu*n_step)]*w_input - ones(nu*n_step*2, 1)*u_max ;

% optimization variables and constraints dimensions
ng = ns*(n_step+2);
nh = length(h_sym);
lambda_init = zeros(ng, 1);
mu_init = zeros(nh, 1);
lambda = sym('lambda', [ng; 1]);
mu = sym('mu', [nh; 1]);

nablah_sym = jacobian(h_sym, w).';

% generate the matlab functions
matlabFunction(h_sym, 'vars', {w}, 'file', 'h');
matlabFunction(nablah_sym, 'vars', {w}, 'file', 'nablah');

w_ = w_init;
lambda_ = lambda_init;
mu_ = mu_init;
sigma_ = sigma_init;


B_ = B(w_,lambda_, mu_);
nablaf_ = nablaf(w_);
nablah_ = nablah(w_);
[g_, nablag_] = g(w_);
f_ = f(w_);
h_ = h(w_);
m1_ = m1(w_, sigma_);

nablaLagrangian_ = nablaLagrangian(w_,lambda_, mu_);

kkt_violation = norm([nablaLagrangian_; g_], inf);

w_history = [w_];
kkt_violation_history = [kkt_violation];
alpha_history = [];

while kkt_violation > tol
    
    opts.ConvexCheck = 'off';
    [deltaw_,~,~,~,multipliers_] = quadprog(B_, nablaf_, nablah_.', -h_, nablag_.', -g_, [], [], [], opts);
    lambda_plus = multipliers_.eqlin;
    mu_plus = multipliers_.ineqlin;

    switch linesearch
        case 'MERIT'
            % perform linesearch with merit function
            nablam1_ = nablaf_.' * deltaw_ - sigma_*norm(g_, 1) - sigma_*norm(max(h_, 0), 1);
            alpha = linesearch_merit(w_, sigma_, m1_, nablam1_,deltaw_);
        case 'ARMIJO'
            % perform linesearch with Armijo condition
            alpha = linesearch_armijo(w_, f_, nablaf_,deltaw_);
        otherwise
            % perform linesearch with Armijo condition
            alpha = linesearch_armijo(w_, f_, nablaf_,deltaw_);
    end


    w_ = w_ + alpha*deltaw_;
    lambda_ = (1-alpha)*lambda_ + alpha*lambda_plus;
    mu_ = (1-alpha)*mu_ + alpha*mu_plus;



    
    B_ = B(w_,lambda_, mu_);
    nablaf_ = nablaf(w_);
    nablah_ = nablah(w_);
    [g_, nablag_] = g(w_);
    h_ = h(w_);
    f_ = f(w_);
    nablaLagrangian_ = nablaLagrangian(w_,lambda_, mu_);
%     if (sigma_coeff*norm(lambda_,inf) > sigma_) 
%         sigma_ = sigma_coeff*norm(lambda_,inf);
%     end
    
    if (norm(lambda_,inf) > sigma_) 
        sigma_ = norm(lambda_,inf)+1;
    
    end
    
    m1_ = m1(w_, sigma_);
    kkt_violation = norm([nablaLagrangian_; g_], inf);
    
    
    disp("-------------------------------------------------------------")
    disp("iteration: " + iters)
    disp("KKT violation: " + kkt_violation)
    
    disp("w: ")
    disp(w_)
    disp("lambda: " + lambda_)
    disp("cost: " + f_)
    disp("alpha: " + alpha)
    disp("m1: " + m1_)
    disp("sigma: " + sigma_)
    
    w_history = [w_history, w_];
    kkt_violation_history = [kkt_violation_history, kkt_violation];
    alpha_history = [alpha_history, alpha];
    
    
    iters = iters + 1;
    
end
%%


%extract w_state
state_trajectory = w_(1:ns);
for i = 1:n_step

    state_trajectory = [state_trajectory, w_(i*(ns+nu)+1:i*(ns+nu)+ns)];
end


figure(1)
title("state trajectory")
plot(times, state_trajectory', 'lineWidth', 1.5)
xlabel("Time")
grid on
legend(["$w$", "$\theta$", "$v$", "$\omega$"])

figure(2)
input_trajectory = w_(ns+1:ns+nu);
for i = 1:n_step-1
    input_trajectory = [input_trajectory; w_(i*(ns+nu)+ns+1:(i+1)*(ns+nu))];
end
plot(input_trajectory, 'lineWidth', 1.5, 'Marker', 'x')
xlabel("Iteration")
title("input trajectory")

figure(3)
plot(alpha_history, 'lineWidth', 1.5, 'Marker', 'x')
xlabel("Iteration")
title("alpha history")

figure(4)
semilogy(kkt_violation_history, 'lineWidth', 1.5), grid on
xlabel("Iteration")
title("KKT violation")

commandwindow
clear all

set(groot, 'defaultAxesTickLabelInterpreter','latex'); set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultColorbarTickLabelInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(groot, 'defaultFigureRenderer','painters');

addpath('..')

global n_step d_step t_init s_init s_fin ns nu n_int int_step state_coeff input_coeff
s_init = [0;0;0;0];
s_fin = [2;0;0;0];
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
state_coeff = 0.01*state_coeff;
input_coeff = 0.001;
% optimization horizon
t_init = 0;
t_fin = 2;
n_step = 20;
times = linspace(t_init,t_fin,n_step+1);
d_step = (t_fin-t_init)/n_step;
u_guess = zeros(n_step * nu, 1);
% integration horizon
% each optimization interval do n_int integration steps of duration
% int_step
n_int = 20;
int_step = d_step / n_int;

ns = length(s_init);
syms t 
u = sym('u', [nu 1]);
s = sym('s', [ns 1]);

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

iters = 1;

% INPUT
tol = 1e-4;
w_init = u_guess;
nw = length(w_init);
sigma_coeff = 2;
sigma_init = 0;
damping_coeff = 0.5;

linesearch = 'MERIT';

w = sym('w', [nw; 1]);
sigma = sym('sigma');

% set the equality constraints
h_sym = [eye(nw); -eye(nw)]*w - ones(nw*2, 1)*u_max ;

% optimization variables and constraints dimensions
ng = 4;
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
nablag_ = nablag(w_);
nablah_ = nablah(w_);
g_ = g(w_);
f_ = f(w_);
h_ = h(w_);
m1_ = m1(w_, sigma_);

nablaLagrangian_ = nablaLagrangian(w_,lambda_, mu_);

kkt_violation = norm([nablaLagrangian_; g_; max(zeros(nh, 1), h_)], inf);

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
            %nablam1_ = nablaf_.' * deltaw_ - sigma_*norm(g_, 1);
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
    nablag_ = nablag(w_);
    nablah_ = nablah(w_);
    g_ = g(w_);
    h_ = h(w_);
    f_ = f(w_);
    nablaLagrangian_ = nablaLagrangian(w_,lambda_, mu_);
    
    
%     if (sigma_coeff*lambda_ > sigma_) 
%          sigma_ = sigma_coeff*lambda_;
%     end
    
    sigma_ = sigma_init;
    if (norm(lambda_,inf) > sigma_) 
        sigma_ = norm(lambda_,inf)+0.01;
    
    end
    
    if (alpha<1e-4)
        sigma_ = sigma_init;
    end
    
    m1_ = m1(w_, sigma_);
    kkt_violation = norm([nablaLagrangian_; g_; max(zeros(nh, 1), h_)], inf);
    
    
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

state_trajectory = [s_init];
for k=1:n_step
    tk = t_init + k*d_step;
    [s,~,~] = expl_rk4(tk,state_trajectory(:,k),w_(k),eye(ns),zeros(ns,nu),int_step,n_int);
    state_trajectory = [state_trajectory, s];
end

figure(1)
title("state trajectory")
plot(times, state_trajectory.', 'lineWidth', 1.5)
xlabel("Time")
grid on
legend(["$w$", "$\theta$", "$v$", "$\omega$"])
%saveas(gcf,'1_state','epsc')
    
figure(2)
plot(times(1:end-1), w_history(:,1:end-1), 'lineWidth', 1.5, 'Color',[0, 0.4470, 0.7410]), hold on
plot(times(1:end-1), w_history(:,end), 'lineWidth', 2, 'Marker', 'x', 'Color', [0.8500, 0.3250, 0.0980])
xlabel("Time")
grid on
title("input evolution over the iterations")
%saveas(gcf,'1_input_evolution','epsc')

figure(3)
plot(alpha_history, 'lineWidth', 1.5, 'Marker', 'x')
xlabel("Iteration")
title("$\alpha$ linesearch")
grid on
ylim([0 1])
%saveas(gcf,'1_alpha','epsc')

figure(4)
semilogy(kkt_violation_history, 'lineWidth', 1.5), grid on
xlabel("Iteration")
title("KKT violation")


%%

%{
times_20 = linspace(t_init,t_fin,20+1);
times_40 = linspace(t_init,t_fin,40+1);
times_80 = linspace(t_init,t_fin,80+1);

state_20 = load('N_20', 'state_trajectory').state_trajectory;
state_40 = load('N_40', 'state_trajectory').state_trajectory;
state_80 = load('N_80', 'state_trajectory').state_trajectory;

input_20 = load('N_20', 'w_history').w_history;
input_40 = load('N_40', 'w_history').w_history;
input_80 = load('N_80', 'w_history').w_history;

alpha_20 = load('N_20', 'alpha_history').alpha_history;
alpha_40 = load('N_40', 'alpha_history').alpha_history;
alpha_80 = load('N_80', 'alpha_history').alpha_history;


figure(1)
title("state trajectory")
plot(times_20, state_20.', 'lineWidth', 1.5, 'Color',[0, 0.4470, 0.7410]), hold on
plot(times_40, state_40.', 'lineWidth', 1.5, 'Color',[0.8500, 0.3250, 0.0980])
plot(times_80, state_80.', 'lineWidth', 1.5, 'Color',[0.9290, 0.6940, 0.1250])
xlabel("Time")
legend("20","","","","40","","","","80")
grid on

figure(2)
plot(times_20(1:end-1), input_20(:,end), 'lineWidth', 1.5, 'Color',[0, 0.4470, 0.7410]), hold on
plot(times_40(1:end-1), input_40(:,end), 'lineWidth', 1.5, 'Color',[0.8500, 0.3250, 0.0980])
plot(times_80(1:end-1), input_80(:,end), 'lineWidth', 1.5, 'Color',[0.9290, 0.6940, 0.1250])
xlabel("Time")
grid on
legend("20","40","80")
title("input trajectory")

figure(3)
plot(alpha_20, 'lineWidth', 1.5, 'Marker', 'x', 'Color',[0, 0.4470, 0.7410]), hold on
plot(alpha_40, 'lineWidth', 1.5, 'Marker', 'x', 'Color',[0.8500, 0.3250, 0.0980])
plot(alpha_80, 'lineWidth', 1.5, 'Marker', 'x', 'Color',[0.9290, 0.6940, 0.1250])
xlabel("Iteration")
grid on
legend("20","40","80")
title("$\alpha$ linesearch")


%}




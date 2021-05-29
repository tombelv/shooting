commandwindow

s_init = [0;0;0;0];
s_fin = [2;0;0;0];
nu = 1;
% problem parameters
M = 1;
m = 1;
l = 1;
g = 9.81;
u_max = 20;
% optimization horizon
t_init = 0;
t_fin = 2;
n_step = 20;
d_step = (t_fin-t_init)/n_step;
u_guess = zeros(n_step, nu);
% integration horizon
% each optimization interval do n_int integration steps of duration h
n_int = 10;
h = d_step / n_int;

ns = length(s_init);
syms t 
u = sym('u', [nu 1]);
s = sym('s', [ns 1]);

F_sym = [s(3); 
    s(4); 
    (m*l*sin(s(2))*s(4)^2 + m*g*cos(s(2))*sin(s(2)) + u) / (M+m-m*(cos(s(2)))^2);
    -((m*l*cos(s(2))*sin(s(2))*s(4)^2 + (M+m)*g*sin(s(2)) + u*cos(s(2))) / (l*(M+m-m*(cos(s(2)))^2)))];
dFds_sym = jacobian(F_sym, s);
dFdu_sym = jacobian(F_sym, u);



matlabFunction(F_sym, 'vars', {t,s,u}, 'file', 'dynamics');
matlabFunction(dFds_sym, 'vars', {t,s,u}, 'file', 'dFds');
matlabFunction(dFdu_sym, 'vars', {t,s,u}, 'file', 'dFdu');

% only for expl rk4
k1_sym = dynamics(t,s,u);
k2_sym = dynamics(t + 1/2 * h, s + 1/2 * h * k1_sym, u);
k3_sym = dynamics(t + 1/2 * h, s + 1/2 * h * k2_sym, u);
k4_sym = dynamics(t + h, s + h * k3_sym, u);
s_sym = s + 1/6 * h * (k1_sym +2*k2_sym + 2*k3_sym + k4_sym);
dsds_sym = jacobian(s_sym,s);
dsdu_sym = jacobian(s_sym,u);

matlabFunction(s_sym, 'vars', {t,s,u}, 'file', 's');
matlabFunction(dsds_sym, 'vars', {t,s,u}, 'file', 'dsds');
matlabFunction(dsdu_sym, 'vars', {t,s,u}, 'file', 'dsdu');

%:)
% [s,A,B] = expl_euler(t_init,s_init,u_init,eye(ns),zeros(ns,nu),h,n_int)
% [x,A,B] = expl_euler_numerical_end(t_init,s_init,u_init,eye(ns),zeros(ns,nu),h,n_int)
% 

% write S_bar adding elements
dxdw = zeros(ns, n_step * nu);
T_bar = [];
S_bar = [];
s_hat =[];
% single shooting
for k = 0:(n_step - 1)
    u_init = u_guess(k+1, :);
    tk = t_init + d_step*k;
    [s,A,fu] = expl_rk4(tk,s_init,u_init,eye(ns),zeros(ns,nu),h,n_int);
    T_bar = [T_bar; A];
    dxdw = A * dxdw + [zeros(ns, k * nu) fu zeros(ns, (n_step-k-1) * nu)];
    S_bar = [S_bar; dxdw];
    % update s_init and u_init
    s_init = s;
    s_hat = [s_hat; s];
end
iters = 1;
Q = [100 0 0 0;
    0 1 0 0;
    0 0 0.1 0;
    0 0 0 0.1];
Q = 0.01*Q;
R = 0.001;
Q_bar = [];
for i = 1:n_step
   Q_bar = [Q_bar; zeros(ns, (i-1) * ns) Q zeros(ns, (n_step-i)* ns)];
end
R_bar = [];
for i = 1:n_step
   R_bar = [R_bar; zeros(nu, (i-1) * nu) R zeros(nu, (n_step-i)* nu)];
end
% INPUT
tol = 1e-8;
w_init = u_guess;
nw = length(w_init);
sigma_coeff = 2;
sigma_init = 1;
damping_coeff = 0.5;

hessian_approx = 'EXACT';
linesearch = 'MERIT';

w = sym('w', [nw; 1]);
sigma = sym('sigma');

% set the cost symbolic expression f_sym as a function of w
f_sym = w.'*(S_bar.'*Q_bar*S_bar + R_bar)*w + s_hat.'*(2*Q_bar*S_bar)*w + s_init.'*(2*T_bar.'*Q_bar*S_bar)*w;
% set the equality constraints
g_sym = s_hat(end-ns+1:end,:) + S_bar(end-ns+1:end,:)*w + T_bar(end-ns+1:end,:)*s_init - s_fin;
h_sym = [eye(nw); -eye(nw)]*w - ones(nw*2, 1)*u_max;
% set merit function
m1_sym = f_sym + sigma * norm(g_sym, 1);

% optimization variables and constraints dimensions
ng = length(g_sym);
nh = length(h_sym);
lambda_init = zeros(ng, 1);
mu_init = zeros(nh, 1);
lambda = sym('lambda', [ng; 1]);
mu = sym('mu', [nh; 1]);

% compute Lagrangian and gradients
nablaf_sym = gradient(f_sym, w);
nablag_sym = jacobian(g_sym, w).';
nablah_sym = jacobian(h_sym, w).';

lagrangian_sym = f_sym + lambda.'*g_sym + mu.'*h_sym;
nablaLagrangian_sym = gradient(lagrangian_sym, w);

% compute the hessian B according to the chosen hessian approximation
switch hessian_approx
    case 'EXACT'
        B_sym = jacobian(jacobian(lagrangian_sym,w),w);
    case 'GAUSS_NEWTON'
%         nablaR = jacobian(R_sym,w);
%         B_sym = nablaR*(nablaR.');
    otherwise
%         disp("defaulted to EXACT hessian")
%         B_sym = jacobian(jacobian(lagrangian_sym,x),x);
end


% generate the matlab functions
matlabFunction(f_sym, 'vars', {w}, 'file', 'f');
matlabFunction(g_sym, 'vars', {w}, 'file', 'g');
matlabFunction(h_sym, 'vars', {w}, 'file', 'h');
matlabFunction(m1_sym, 'vars', {w, sigma}, 'file', 'm1');
matlabFunction(nablaf_sym, 'vars', {w}, 'file', 'nablaf');
matlabFunction(nablag_sym, 'vars', {w}, 'file', 'nablag');
matlabFunction(nablah_sym, 'vars', {w}, 'file', 'nablah');
matlabFunction(nablaLagrangian_sym, 'vars', {w, lambda, mu}, 'file', 'nablaLagrangian');
matlabFunction(B_sym, 'vars', {w, lambda, mu}, 'file', 'B');



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

x_history = [w_];
kkt_violation_history = [kkt_violation];
alpha_history = [];

while kkt_violation > tol
    % regularization of the hessian
    B_ = hessian_regularization(nablag_, B_)  
    
    opts.ConvexCheck = 'off';
    [deltaw_,~,~,~,multipliers_] = quadprog(B_, nablaf_, nablah_.', -h_, nablag_.', -g_, [], [], [], opts);
    lambda_plus = multipliers_.eqlin;
    mu_plus = multipliers_.ineqlin;

    switch linesearch
        case 'MERIT'
            % perform linesearch with merit function
            nablam1_ = nablaf_.' * deltaw_ - sigma_*norm(g_, 1);
            alpha = linesearch_merit(x_, sigma_, m1_, nablam1_,deltaw_);
        case 'ARMIJO'
            % perform linesearch with Armijo condition
            alpha = linesearch_armijo(x_, f_, nablaf_,deltaw_);
        otherwise
            % perform linesearch with Armijo condition
            alpha = linesearch_armijo(x_, f_, nablaf_,deltaw_);
    end

    

    w_ = w_ + alpha*deltaw_;
    lambda_ = (1-alpha)*lambda_ + alpha*lambda_plus;
    mu_ = (1-alpha)*mu_ + alpha*mu_plus;


    
    B_ = B(w_,lambda_, mu_);
    nablaf_ = nablaf(w_);
    nablag_ = nablag(w_);
    g_ = g(w_);
    f_ = f(w_);
    nablaLagrangian_ = nablaLagrangian(w_,lambda_, mu_);
    if (sigma_coeff*lambda_ > sigma_) 
        sigma_ = sigma_coeff*lambda_;
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
    
    w_history = [w_history, w_];
    kkt_violation_history = [kkt_violation_history, kkt_violation];
    alpha_history = [alpha_history, alpha];
    
    
    iters = iters + 1;
    
    
end
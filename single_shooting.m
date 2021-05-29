commandwindow

s_init = [2;0;0;0];
nu = 1;
% problem parameters
M = 1;
m = 1;
l = 1;
g = 9.81;
% optimization horizon
t_init = 0;
t_fin = 2;
n_step = 20;
d_step = (t_fin-t_init)/n_step;
u_guess = zeros(nu, n_step);
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
% single shooting
for k = 0:(n_step - 1)
    u_init = u_guess(:, k+1);
    tk = t_init + d_step*k;
    [s,A,B] = expl_rk4(tk,s_init,u_init,eye(ns),zeros(ns,nu),h,n_int);
    T_bar = [T_bar; A];
    dxdw = A * dxdw + [zeros(ns, k * nu) B zeros(ns, (n_step-k-1) * nu)];
    S_bar = [S_bar; dxdw];
    % update s_init and u_init
    s_init = s;
end
commandwindow

s_init = [0;0;0;0];
nu = 1;
% optimization horizon
t_init = 0;
t_fin = 2;
n_step = 20;
d_step = (t_fin-t_init)/n_step;
u_guess = zeros(nu, n_step-1);
% integration horizon
% each optimization interval do n_int integration steps of duration h
n_int = 10;
h = d_step / n_int;

ns = length(s_init);
syms t 
u = sym('u', [nu 1]);
s = sym('s', [ns 1]);

F_sym = [-16*s(1)+12*s(2)+16*cos(t)-13*sin(t)+u; 16*s(1)-9*s(2)-11*cos(t)+9*sin(t)+u];
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
S_bar = zeros(nx, nu);
T_bar = [];

% single shooting
for k = 0:(n_step - 1)
    u_init = u_guess(k, :);
    tk = t_init + d_step*k;
    [~,fx,fu] = F(tk,s_init,u_init);
    [s,A,~] = expl_rk4(tk,s_init,u_init,eye(ns),zeros(ns,nu),h,n_int);
    T_bar = [T_bar; A];
    dxdw = fx * S_bar(end-nx+1:end, :) + [zeros(nx, k * nu); fu; zeros(nx, (n_step-k-1) * nu)];
    S_bar = [S_bar; dxdw];
    % update s_init and u_init
    s_init = s;
end
% remove from S_bar elements relative to x0
S_bar = S_bar(nx + 1:end, :);
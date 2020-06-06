addpath('/Users/adi/Downloads/casadi-osx-matlabR2015a-v3.5.1')

import casadi.*

% Declare model variables

params;

T = 1; % Time horizon
N = 20; % number of control intervals

x1 = SX.sym('x1');
x2 = SX.sym('x2');
x3 = SX.sym('x3');
x4 = SX.sym('x4');
x = [x1; x2; x3; x4];
u = SX.sym('u');

q1 = x(1);
q2 = x(2);
dq1 = x(3);
dq2 = x(4);

D = [ M2*l1^2 + M2*cos(q2)*l1*l2 + (M2*l2^2)/4 + J1 + J2,            (M2*l2^2)/4 + (M2*l1*cos(q2)*l2)/2 + J2;
     (M2*l2^2)/4 + (M2*l1*cos(q2)*l2)/2 + J2,                        (M2*l2^2)/4 + J2];

C = [(M2*g*l2*sin(q1 + q2))/4 + (M1*g*l1*sin(q1))/4 + (M2*g*l1*sin(q1))/2 - (M2*dq2^2*l1*l2*sin(q2))/2 - M2*dq1*dq2*l1*l2*sin(q2);
                                                                              (M2*l2*(2*l1*sin(q2)*dq1^2 + g*sin(q1 + q2)))/4];

d2X = D\(-C+u);

dx = [x(3);x(4);d2X];


% Model equations

xdot = dx;


% Objective term

L = sqrt(((x1-pi)^2)+((x2-pi)^2) + 1e-5);


% Fixed step Runge-Kutta 4 integrator

   M = 4; % RK4 steps per interval
   DT = T/N/M;
   f = Function('f', {x, u}, {xdot, L});
   X0 = MX.sym('X0', 4);
   U = MX.sym('U');
   X = X0;
   Q = 0;
   for j=1:M
       [k1, k1_q] = f(X, U);
       [k2, k2_q] = f(X + DT/2 * k1, U);
       [k3, k3_q] = f(X + DT/2 * k2, U);
       [k4, k4_q] = f(X + DT * k3, U);
       X=X+DT/6*(k1 +2*k2 +2*k3 +k4);
       Q = Q + DT/6*(k1_q + 2*k2_q + 2*k3_q + k4_q);
   end
    F = Function('F', {X0, U}, {X, Q}, {'x0','p'}, {'xf', 'qf'})
    
    
% Start with an empty NLP

w={};
w0 = [];
lbw = [];
ubw = [];
J = 0;
g={};
lbg = [];
ubg = [];


% Formulate the NLP

Xk = [0; 0; 0; 0];
for k=0:N-1
    
    % New NLP variable for the control

    Uk = MX.sym(['U_' num2str(k)]);
    w = {w{:}, Uk};
    lbw = [lbw, -inf];
    ubw = [ubw, inf];
    w0 = [w0, 0];

    % Integrate till the end of the interval
    
    Fk = F('x0', Xk, 'p', Uk);
    Xk = Fk.xf;
    J=J+Fk.qf;

    % Add inequality constraint
    
    g = {g{:}, Xk(1)};
    lbg = [lbg; 0];
    ubg = [ubg; pi];
end

%J = Xk(1);

J = (Xk(1)-pi)^2+(Xk(2)-0)^2;

% Create an NLP solver

prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));
solver = nlpsol('solver', 'ipopt', prob);

% Solve the NLP

sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw, ...
             'lbg', lbg, 'ubg', ubg);
w_opt = full(sol.x);


% Plot the solution

u_opt = w_opt;
x_opt = [0;0;0;0];
for k=0:N-1
    Fk = F('x0', x_opt(:,end), 'p', u_opt(k+1));
    x_opt = [x_opt, full(Fk.xf)];
end
x1_opt = x_opt(1,:);
x2_opt = x_opt(2,:);
tgrid = linspace(0, T, N+1);
% tgrid = tgrid(:,1:end-1);
clf;

% subplot(1,2,1)
% plot(tgrid, x1_opt, '--')
% hold on;
% plot(tgrid, x2_opt, '-')
% hold on;
% stairs(tgrid, [u_opt], '-.')
% xlabel('t')
% legend('x1','x2','u')
% hold off

% [tgrid' x1_opt' x2_opt' u_opt]

% subplot(1,2,2)
animate(tgrid,[x1_opt' x2_opt'])

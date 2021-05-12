%dual implementation of peak estimation of the disturbed flow system

% for debugging purposes only

%compare against meas_class/testing/flow_test_b.m

%% Initiation

%variable declaration
t = sdpvar(1,1); %just to be safe
x = sdpvar(2,1);
w = sdpvar(1,1);
vars = [t; x; w];

%solving parameters
order = 4;
d = 2*order;

%% Build Support Sets

%time t
Tsupp = struct('ineq', t*(1-t), 'eq', 0);
Tmax = 5;

%state x
% Xsupp = struct('ineq', [16-(x(1)+1)^2; 3.5^2 - (x(2)+0.25)^2], 'eq', []);
% Xsupp = struct('ineq', 9-x.^2, 'eq', 0);
box = [-1, 3; -1.5, 2];
[bo,bc,bh]=box_process(2,box);
Xsupp = struct('ineq', bh.^2 - (x-bc).^2, 'eq', []);

% objective = -x(2);
% objective = -sum(x);
objective = -x;

C0 = [1.5; 0];
R0 = 0.4;

X0 = struct('ineq', R0^2 - sum((x-C0).^2), 'eq', []);

%disturbance w
Wsupp = struct('ineq', w*(1-w), 'eq', []);


%dynamics
f0 = Tmax * [x(2); -0.85*x(1)-x(2)+(x(1)^3)/3];
f1 = Tmax * [0; -0.3*x(1)];

%in line with dynamics from most recent presentation

%% Build Polynomials
%auxiliary function and bounds
[v, cv] = polynomial([t;x],d);
gamma = sdpvar(1,1);

cons = [];
coeff_list = [gamma; cv];

%psatz time
v0 = replace(v, t, 0);
[p0, cons0, coeff0] = constraint_psatz(gamma - v0, X0, [x], d);

%lie derivative (ignore switching for now
Xf = struct('ineq', [Tsupp.ineq; Xsupp.ineq; Wsupp.ineq], 'eq', []);
Xall = struct('ineq', [Tsupp.ineq; Xsupp.ineq], 'eq', []);

Lv = jacobian(v, x)*(f0 + w*f1) + jacobian(v, t);
[pf, consf, coefff] = constraint_psatz(-Lv, Xf, [t;x; w], d);

%cost
% [pc, consc, coeffc] = constraint_psatz(v - objective, Xf, [t; x], d);
beta = sdpvar(length(objective), 1);
% [pc, consc, coeffc] = constraint_psatz(v - beta'*objective, Xall, [t; x], d);
[pc, consc, coeffc] = constraint_psatz(v - beta'*objective, Xall, [t; x], d);
consc = [consc; sum(beta)==1; beta >= 0];
coeffc = [coeffc; beta];

    
cons = [cons0; consf; consc];
coeff_list = [coeff_list; coeff0; coefff; coeffc];

opts = sdpsettings('solver', 'mosek');
opts.sos.model = 2;

[sol, monom, Gram, residual] = solvesos(cons, gamma, opts, coeff_list);
peak_val = value(gamma)
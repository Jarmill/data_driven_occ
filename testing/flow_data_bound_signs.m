%dual implementation of peak estimation of the disturbed flow system
%
% use the prior knowledge that w is in [0.85, 1.15] through the Convex
% Duality method. This is a specific instance of data-driven estimation,
% and is distinct from the box method

%the correct signs are used for the dual multipliers

% for debugging purposes only

%% Initiation

%variable declaration
t = sdpvar(1,1); %just to be safe
x = sdpvar(2,1);

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

objective = -x(2);

C0 = [1.5; 0];
R0 = 0.4;

X0 = struct('ineq', R0^2 - sum((x-C0).^2), 'eq', 0);

%dynamics
% f0 = Tmax * [x(2); -0.85*x(1)-x(2)+(x(1)^3)/3];
% f1 = Tmax * [0; -0.3*x(1)];


f0 = Tmax * [x(2); -x(1)-x(2)+(x(1)^3)/3];
f1 = Tmax * [0; x(1)];

w_bound = [-0.15;0.15];
w_lo = w_bound(1);
w_hi = w_bound(2);

%dual constarint
A = [1; -1];
b = [w_lo; w_hi];


%in line with dynamics from most recent presentation

%% Build Polynomials
%auxiliary function and bounds
[v, cv] = polynomial([t;x],d);
gamma = sdpvar(1,1);

%alternative polynomials
d_altern = 2*ceil(d/2 + degree(f1)/2); %figure out degree bounds later
[qlo, cqlo] = polynomial([t; x], d_altern);
[qhi, cqhi] = polynomial([t; x], d_altern);


cons = [];
coeff_list = [gamma; cv; cqlo; cqhi];

%psatz time
v0 = replace(v, t, 0);
[p0, cons0, coeff0] = constraint_psatz(gamma - v0, X0, [x], d);

Xall = struct('ineq', [Tsupp.ineq; Xsupp.ineq], 'eq', []);

%cost
[pc, consc, coeffc] = constraint_psatz(v - objective, Xall, [t; x], d);


%% Dual Constraints



%nonnegativity of alternative functions
[plo, conslo, coefflo] = constraint_psatz(qlo, Xall, [t;x], d_altern);
[phi, conshi, coeffhi] = constraint_psatz(qhi, Xall, [t;x], d_altern);


%lie derivative constraint
Lv0 = jacobian(v, x)*(f0) + jacobian(v, t);
Lv1 = jacobian(v, x)*(f1);


[pf0, consf0, coefff0] = constraint_psatz(Lv0-b'*[qlo;qhi], Xall, [t;x], d);

% consf1 = (coefficients(Lv1, cv) == coefficients(A'*[qlo; qhi], [cqlo;cqhi]));
consf1 = (coefficients(Lv1 + A'*[qlo; qhi], [t;x])==0);

%todo: replace

cons = [cons0; consc; conslo; conshi; consf0; consf1];
coeff_list = [coeff_list; coeffc; coefflo; coeffhi; coeff0; coefff0];

opts = sdpsettings('solver', 'mosek');
opts.sos.model = 2;

[sol, monom, Gram, residual] = solvesos(cons, gamma, opts, coeff_list);
peak_val = value(gamma)


%% Define quantities
t = sdpvar(1, 1);
x = sdpvar(2, 1);

Tmax = 5;

f0 = Tmax * [x(2); -0.85*x(1)-x(2)+(x(1)^3)/3];
f1 = Tmax * [0; -0.3*x(1)];


C0 = [1.5; 0];
R0 = 0.4;


X0 = struct('ineq', R0^2 - sum((x-C0).^2), 'eq', 0);

% objective = -x(2);
objective = x;


%dual constarint Aw <= b
A = [1; -1];
b = [1; 0];

W = struct('A', A, 'b', b);

%% fill in location support
lsupp = loc_sos_options();
lsupp.t = t;
lsupp.x = x;
lsupp = lsupp.set_box([-1, 3; -1.5, 2]);
lsupp.X_init = X0;
lsupp.f0 = f0;
lsupp.fw = f1;
lsupp.W = W;
% lsupp.Tmax = Tmax;
lsupp.Tmax = 1;


lsupp.verbose = 1;


%% start up tester
PM = peak_sos(lsupp, objective);

order = 4;
d = 2*order;

[prog]= PM.make_program(d);
out = PM.solve_program(prog)

% [poly_out, coeff_out] = PM.make_poly(d);
% [coeff_lie, con_lie] = PM.make_cons(d, poly_out)
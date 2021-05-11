
%% Define quantities
t = sdpvar(1, 1);
x = sdpvar(2, 1);

Tmax = 5;

% f0 = [x(2); -x(1)-x(2)+(x(1)^3)/3];
% f0 = [2*x(1); - x(2)];
f0 = [0; -1];

% f1 = {[0; -0.3*x(1)]};


%% Add noise
NOISE = 1;
if NOISE == 2
    %rows: dynamics
    %columns: number of uncertain inputs

    fw = [0 1; x(1) 0];
    W = struct('A', kron(eye(2), [1; -1]), 'b', 0.15*[-1; -1; -10; -10]);
elseif NOISE == 1
    fw = [0; x(1)];
    W = struct('A', [1; -1], 'b', [0.15; -0.15]);
else
    fw = [];
    W = [];
end


C0 = [1.5; 0];
R0 = 0.4;

X0 = struct('ineq', R0^2 - sum((x-C0).^2), 'eq', 0);

objective = -x(2);
% objective = -x;

%% fill in location support
lsupp = loc_sos_options();
lsupp.t = t;
lsupp.x = x;
lsupp = lsupp.set_box([-1, 3; -1.5, 2]);
lsupp.X_init = X0;
lsupp.f0 = f0;
lsupp.fw = fw;
lsupp.W = W;
lsupp.Tmax = Tmax;


lsupp.verbose = 1;


%% start up tester
PM = peak_sos(lsupp, objective);

order = 2;
d = 2*order;

[prog]= PM.make_program(d);
out = PM.solve_program(prog)

% [poly_out, coeff_out] = PM.make_poly(d);
% [coeff_lie, con_lie] = PM.make_cons(d, poly_out)
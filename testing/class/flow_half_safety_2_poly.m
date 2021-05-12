%

clear
%% Define quantities
t = sdpvar(1, 1);
x = sdpvar(2, 1);

Tmax = 5;
dmax = 0.15;
f0 = [x(2); -x(1)-x(2)+(x(1)^3)/3];
% f1 = {[0; -0.3*x(1)]};


%% Add noise
NOISE = 1;
if NOISE == 2
    %rows: dynamics
    %columns: number of uncertain inputs

    fw = [0 1; x(1) 0];
    W = struct('A', kron(eye(2), [1; -1]), 'b', 0.15*[-1; -1; -10; -10]);
elseif NOISE == 1
    fw = [0 0; dmax*x(1) 1];
    W = struct('A', [1; -1], 'b', [1; 1]);
else
    
    fw = [];
    W = [];
end




C0 = [1.5; 0];
R0 = 0.4;
INIT_POINT = 0;
if INIT_POINT
    X0 = C0;
else
    X0 = struct('ineq', R0^2 - sum((x-C0).^2), 'eq', 0);
end


%unsafe set
Cu = [0; -0.5];
%Cu = [2.5; 0];
Ru = 0.5;
c1f = Ru^2 - (x(1) - Cu(1)).^2 - (x(2) - Cu(2)).^2;

theta_c = 5*pi/4;       %p* = -0.1417, beta = [0, 1]
w_c = [cos(theta_c); sin(theta_c)];
c2f = w_c(1)*(x(1) - Cu(1)) + w_c(2) * (x(2) - Cu(2)); 

% objective = -x(2);
objective = [c1f; c2f];
% objective = -x;

%% fill in location support
lsupp = loc_sos_options();
lsupp.t = t;
lsupp.x = x;
lsupp = lsupp.set_box([-1, 3; -1.5, 2]);
% lsupp = lsupp.set_box(3);
lsupp.X_init = X0;
lsupp.f0 = f0;
lsupp.fw = fw;
lsupp.W = W;
lsupp.Tmax = Tmax;


lsupp.verbose = 1;


%% start up tester
PM = peak_sos(lsupp, objective);

order = 4;
d = 2*order;

[prog]= PM.make_program(d);
out = PM.solve_program(prog)

% [poly_out, coeff_out] = PM.make_poly(d);
% [coeff_lie, con_lie] = PM.make_cons(d, poly_out)
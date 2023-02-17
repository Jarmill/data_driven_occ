%dual implementation of peak estimation of the disturbed flow system
%perturbed by semidefinite-constrained noise
% for debugging purposes only

%% Initiation

%variable declaration
t = sdpvar(1,1); %just to be safe
x = sdpvar(2,1);
% w = sdpvar(3,1);
vars = [t; x];
% w = 0.5;

%solving parameters
% order = 1; %1
% order =2; % 1
% order =3; % 0.8952
% order = 4;%0.8477
order = 5;

d = 2*order;

A0 = eye(3);
A = cell(3, 1);
A{1} = [0 1 0; 1 0 0; 0 0 0];
A{2} = [0 0 1; 0 0 0; 1 0 0];
A{3} = [0 0 0; 0 0 1; 0 1 0];

% build the multiplier
Z = zeros(3, 3, 'like',sdpvar);
coeffZ = [];
for i = 1:3
    for j = 1:i
        [pcurr, cpcurr] = polynomial([t; x], d);
        Z(i, j) = pcurr;
        Z(j, i) = pcurr;
        coeffZ = [coeffZ; cpcurr];
    end
end

%% Build Support Sets

%time t
Tsupp = struct('ineq', t*(1-t), 'eq', 0);
Tmax = 5;

%state x
box = [-0.5, 1.75; -1, 0.5];
[bo,bc,bh]=box_process(2,box);
Xsupp = struct('ineq', bh.^2 - (x-bc).^2, 'eq', []);

objective = -x(2);

C0 = [1.25; 0];
X0 = C0;

%disturbance w
% Wsupp = struct('ineq', w.*(1-w), 'eq', 0);
% Wsupp = struct('ineq', 1-w.^2, 'eq', 0);


%dynamics
f0 = Tmax * [x(2); -x(1)-x(2)+(x(1)^3)/3];
% f1 = Tmax * [0; -0.3*x(1)];
% w_scale = 0.25;
% w_scale = 0.02;
% w_scale = 0; %TODO: bugchecker
fw = Tmax*w_scale*[0, 0, 0; x(1), x(1)*x(2), x(2)];

%% Build Polynomials
%auxiliary function and bounds
[v, cv] = polynomial([t;x],d);
gamma = sdpvar(1,1);

Xall = struct('ineq', [Tsupp.ineq; Xsupp.ineq], 'eq', []);


% cons = [];
coeff_list = [gamma; cv];

%psatz time
% v0 = replace(v, t, 0);
% [p0, cons0, coeff0] = constraint_psatz(gamma - v0, X0, [x], d);

%initial condition
coeff0 = [];
v0_eval = replace(v, [t; x], [0; C0]);
cons0 = (gamma >= v0_eval);

%cost
[pc, consc, coeffc] = constraint_psatz(v - objective, Xall, [t; x], d);

%lie derivative
Lv0 = jacobian(v, t) + jacobian(v, x)*(f0);
Lvcon = -Lv0 - sum(A0.*Z, 'all');
[pf, consf, coefff] = constraint_psatz(Lvcon, Xall, [t; x], d);

%equality constraints
con_eq = [];
for i = 1:3
    ccurr = coefficients(sum(A{i}.*Z, 'all') - jacobian(v, x)*fw(:, i), [t; x]);
    con_eq = [con_eq; ccurr ==0];
end
% Xf = struct('ineq', [Tsupp.ineq; Xsupp.ineq; Wsupp.ineq], 'eq', []);

%% generate the psatz for Z
mom_0 = simple_moments(3, order);
v_monom_0 = recovermonoms(mom_0.monom_int, vars);

mom_g = simple_moments(3, order-1);
v_monom_g = recovermonoms(mom_g.monom_int, vars);

Gram = cell(1+length(Xall.ineq), 1);

[Phi, Gram{1}] = make_Gram(1, 3, mom_0, v_monom_0);
Phi_Z = Phi;

for i = 1:length(Xall.ineq)
    [Phicurr, Gram{i+1}] = make_Gram(Xall.ineq(i), 3, mom_g, v_monom_g);
    Phi_Z = Phi_Z + Phicurr;
end

cons_Gram_Z = [];
Gramvar = [];
for i = 1:length(Gram)
    cons_Gram_Z = [cons_Gram_Z; Gram{i}>=0];
    Gramvar = [Gramvar; reshape(Gram{i}, [], 1)];
end

con_eq_Z = (coefficients(Z - Phi_Z, [t; x]) == 0);

   
%% assemble together the constraints
cons = [cons0:'initial'; consf:'Lie Base'; consc:'cost'; con_eq:'Lie Altern'; con_eq_Z:'Eq Z'; cons_Gram_Z:'Gram Z'];
coeff_list = [coeff_list; coeff0; coefff; coeffc; Gramvar; coeffZ];

opts = sdpsettings('solver', 'mosek');
opts.sos.model = 2;

[sol, monom, Gram_solve, residual] = solvesos(cons, gamma, opts, coeff_list);
peak_val = value(gamma);

fprintf('peak bound %0.4f\n', peak_val)

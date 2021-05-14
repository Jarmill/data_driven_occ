%safety margin of a flow system with polytopic uncertainty


SOLVE = 1;
SAMPLE = 0;
PLOT = 0;

if SOLVE
    
%% Define quantities
t = sdpvar(1, 1);
x = sdpvar(2, 1);

Tmax = 5;
dmax = 0.15;
f0 = [x(2); -(1-dmax)*x(1)-x(2)+(x(1)^3)/3];
% f1 = {[0; -0.3*x(1)]};


%% Add noise
NOISE = 1;
if NOISE == 2
    %rows: dynamics
    %columns: number of uncertain inputs

    fw = [0 dmax; -2*dmax*x(1) 0];
    W = struct('A', kron(eye(2), [1; -1]), 'b', [1; 0; 1; 1]);
%     W = struct('A', kron(eye(2), [1; -1]), 'b', 0.15*[-1; -1; -10; -10]);
elseif NOISE == 1
%     fw = [0; x(1)];
%     W = struct('A', [1; -1], 'b', [-0.15; -0.15]);
    fw = [0; -2*dmax*x(1)];
%     W = struct('A', [1; -1], 'b', [1; 1]);
    W = struct('A', [1; -1], 'b', [1; 0]);
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

objective = -x(2);
% objective = [c1f; c2f];
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

% [prog]= PM.make_program(d);
% out = PM.solve_program(prog)
out = PM.run(order);


% [poly_out, coeff_out] = PM.make_poly(d);
% [coeff_lie, con_lie] = PM.make_cons(d, poly_out)
end

if SAMPLE
    % Nsample = 300;
% Nsample = 150;
% Nsample = 100;
Nsample = 2;
% Nsample = 50;
    rng(33, 'twister')
s_opt = sampler_options;

if INIT_POINT
    s_opt.sample.x = @() C0;
else
    s_opt.sample.x = @() sphere_sample(1, 2)'*R0 + C0;
end
if NOISE==1
    s_opt.sample.d = @() rand();
    s_opt.Nd = 1;
end
    
s_opt.Tmax = lsupp.Tmax;
s_opt.parallel = 0;

% s_opt.mu = 0.1;

%% compatibility with old sampler code

% out.dynamics.nonneg = @(t,x,w,d)cell2mat(arrayfun(@(i)out.dynamics.nonneg_val([t(i);x(:,i)]),...
%     (1:length(t)),'UniformOutput',false));
% out.func.vval = @(t,x,w) cell2mat(arrayfun(@(i) out.func.v([t(i), x(i, :)])',(1:length(t)),'UniformOutput',false));
% out.func.beta = out.poly.beta;
% out.func.cost_all = cell(length(objective), 1);
% for i = 1:length(objective)
%     out.func.cost_all{i} = polyval_func(objective(i), x);
% end

tic
out_sim = sampler(out.dynamics, Nsample, s_opt);

out_sim = traj_eval(out, out_sim);

sample_time = toc;




end

if PLOT
    
    PS = peak_sos_plotter(out, out_sim);
    PS.nonneg_zeta();
    PS.nonneg_traj();
    PS.state_plot();
%     nplot = nonneg_plot(out, out_sim);
%     splot = state_plot_2(out, out_sim);
end

% if PLOT_NONNEG
% end
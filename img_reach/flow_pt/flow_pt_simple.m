%data driven reachable set estimation of a flow system

%break up the sections here into functions

PROBLEM = 1;
SOLVE = 1;
SAMPLE = 1;
EVAL = 0;
PLOT = 1;

if PROBLEM
rng(35, 'twister')
%% generate samples
% A_true = [-1 1; -1 -0.3];
f_true = @(t, x) [x(2); -x(1) + (1/3).* x(1).^3 - x(2)];

% 
% Nsample = 100;
% Nsample = 50;
Nsample = 40;
% Nsample = 30;
% Nsample = 20;
% Nsample = 4;
% box_lim = 2;
Tmax = 5;
% epsilon = 2;
% epsilon = [0; 2.5];
% epsilon = [0; 1];
epsilon = [0; 0.5];

box_lim =  [-0.3, 1.75;-1,0.5];
sample = struct('t', Tmax, 'x', @() box_lim*(2*rand(2,1)-1));



%% generate model
t = sdpvar(1, 1);
x = sdpvar(2, 1);

DG = data_generator(sample);

observed = DG.corrupt_observations(Nsample, f_true, epsilon);
% [model, W] = DG.reduced_model(observed, x, 1, 1);
% model = DG.poly_model(vars, 3);
mlist = monolist(x, 3);
model = struct('f0', [x(2);(1/3).* x(1).^3 - x(2)], 'fw', [0; x(1)]);

W = DG.data_cons(model, x, observed);
[model_cheb,W_cheb] = DG.center_cheb(model, W);
W_red = DG.reduce_constraints(W_cheb);


model = model_cheb;
W = W_red;
[w_handle, box]= DG.make_sampler(W_cheb);

end
 
%% Solve SOS program
if SOLVE
    
    %start at a single point
    C0 = [1.5; 0];
%     C0 = [0; 0.5];
%     C0 = [-1; 0];
    R0 = 0.2;
    INIT_POINT = 1;
    if INIT_POINT
        X0 = C0;
    else
        X0 = struct('ineq', R0^2 - sum((x-C0).^2), 'eq', []);
    end
    
    
    lsupp = loc_sos_options();
    lsupp.t = t;
    lsupp.TIME_INDEP = 0;
    lsupp.x = x;
%     lsupp.X = struct('ineq', 2*box_lim^2
    lsupp = lsupp.set_box(box_lim);
%     lsupp.X = struct('ineq', 2*box_lim^2 - sum(x.^2), 'eq', []);
    % lsupp = lsupp.set_box(3);
    lsupp.X_init = X0;
    lsupp.f0 = model.f0;
    lsupp.fw = model.fw;
    lsupp.W = W;
    lsupp.Tmax = Tmax;

    lsupp.verbose = 1;

    box_supp = box_process(2, box_lim);
    mom_handle = @(d) LebesgueBoxMom(d, box_supp', 1);

%     mom_handle = @(d) LebesgueSphereMom(monpowers(2, d), sqrt(2)*box_lim);
    
    
    %% start up tester
    PM = reach_sos(lsupp, mom_handle);

    order = 6;
%     order = 3;
    d = 2*order;

    % [prog]= PM.make_program(d);
    % out = PM.solve_program(prog)
    out = PM.run(order);
    
end

%% Sample trajectories
if SAMPLE
    
    s_opt = sampler_options;
    if INIT_POINT
        s_opt.sample.x = @() X0;
    else
        s_opt.sample.x = @() R0*ball_sample(1,2)'+C0;
    end
    s_opt.sample.d = w_handle;
    s_opt.Nd = size(model.fw, 2);
    
    s_opt.Tmax = lsupp.Tmax;
    s_opt.parallel = 1;
    
%     Nsample_traj = 10;
    Nsample_traj = 100;
    
    tic
    out_sim = sampler(out.dynamics, Nsample_traj, s_opt);

    sample_time = toc;
end


if EVAL || SAMPLE
    out_sim = traj_eval_reach(out, out_sim);
end

%% plot trajectories
if PLOT
    
%     if PLOT
    
    PS = reach_sos_plotter(out, out_sim);
    PS.nonneg_zeta();
%     PS.obj_plot();
    PS.state_plot();
    PS.v_plot();
    PS.nonneg_traj();
    
    PS.state_plot_2(box_lim);
    if INIT_POINT
        scatter(C0(1), C0(2), 200, 'k')
    else
        theta = linspace(0,2*pi, 200);
        plot(R0*cos(theta)+C0(1), R0*sin(theta)+C0(2), 'color', 'k', 'LineWidth', 3);
    end
    
    %observation plot    
    DG.data_plot_2(observed, 0.2);

    
    
end


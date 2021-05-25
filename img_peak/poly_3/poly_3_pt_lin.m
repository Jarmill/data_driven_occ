%data driven peak estimation of a linear system

%break up the sections here into functions

PROBLEM = 1;
SOLVE = 0;
SAMPLE = 0;
EVAL = 0;
PLOT = 0;


%sample data only from initial set
C0 = [-1; 0; 0];
R0 = 0.2;
x0_sample_ball = @() R0*ball_sample(1,3)'+C0;

INIT_SAMPLE_ONLY = 0;

if PROBLEM
rng(33, 'twister')
%% generate samples
% A_true = [-1 4; -1 -0.3];
% A_true = [-1 4 1; -1 -0.3 1; -1 1 -1];
% A_true = [-1 1; -1 -0.3];
% f_true = @(t, x) A_true*x;

A_true = [-1 1 1; -1 0 -1; 0 1 -2];
B_true = [-1 0 -1;
          0 1 1;
          1 1 0]/2;

f_true = @(t,x) A_true*x - B_true*(4*x.^3 - 3*x);



% Nsample = 150;
Nsample = 100;
% Nsample = 50;
% Nsample = 40;
% Nsample = 30;
% Nsample = 20;
% Nsample = 10;
% Nsample = 4;
box_lim = 1;
Tmax = 8;
epsilon = 0.5;
% epsilon = 2;
% epsilon = 2.5;
sample = struct('t', Tmax, 'x', @() box_lim*(2*rand(3,1)-1));

if INIT_SAMPLE_ONLY
    sample.x = x0_sample_ball;
end

% [observed] = corrupt_observations(Nsample,sample, f_true, epsilon);

%% generate model
t = sdpvar(1, 1);
x = sdpvar(3, 1);

DG = data_generator(sample);

observed = DG.corrupt_observations(Nsample, f_true, epsilon);
model.f0 = -B_true*(4*x.^3 - 3*x);
model.fw = kron(eye(3), x');

W_orig = DG.data_cons(model, x, observed);
W = W_orig;

[model_cheb,W_cheb] = DG.center_cheb(model, W);
box_cheb = poly_bounding_box(W_cheb.A, W_cheb.b);
% W_red = DG.reduce_constraints(W_cheb);



W_red = DG.reduce_constraints(W_cheb);
[w_handle, box]= DG.make_sampler(W_red);
% [model, W] = DG.reduced_model(observed, x, 1, 1);
model = model_cheb;
W = W_red;


% [model, W] = DG.reduced_model(observed, x, 1, 1);

% [w_handle, box]= DG.make_sampler(W);

% DG.data_plot_3(observed, 0.2);

end
 
%% Solve SOS program
if SOLVE
    
    %start at a single point

%     C0 = [0; 0.5];

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
    lsupp = lsupp.set_box(box_lim);
    lsupp = lsupp.set_box(1);
%     lsupp.X = struct('ineq', 3*box_lim^2 - sum(x.^2), 'eq', []);
%     lsupp.X = struct('ineq', box_lim^2 - x.^2, 'eq', []);
    lsupp.X = struct('ineq', [1-x(1)^2; 1-x(2)^2; x(3)*(1-x(3))], 'eq', []);
    
    lsupp.X_init = X0;

    lsupp.f0 = model.f0;
    lsupp.fw = model.fw;
    

    lsupp.W = W;
    lsupp.Tmax = Tmax;
    lsupp.verbose = 1;

    objective = x(3);
    
    %% start up tester
    PM = peak_sos(lsupp, objective);

    order = 2;
    d = 2*order;

    % [prog]= PM.make_program(d);
    % out = PM.solve_program(prog)
    out = PM.run(order);
    
end

%% Sample trajectories
if EVAL
    load('poly_3_lin_sim.mat')
elseif SAMPLE
    
    s_opt = sampler_options;
    if INIT_POINT
        s_opt.sample.x = @() X0;
    else
        s_opt.sample.x = x0_sample_ball;
    end
    s_opt.sample.d = w_handle;
    s_opt.Nd = size(model.fw, 2);
    
    s_opt.Tmax = lsupp.Tmax;
    s_opt.parallel = 0;
    
    Nsample_traj = 104;
%     Nsample_traj = 16;
%     Nsample_traj = 2;
    
    tic
    out_sim = sampler(out.dynamics, Nsample_traj, s_opt);

%     out_sim = traj_eval(out, out_sim);

    sample_time = toc;
end

if EVAL || SAMPLE
    out_sim = traj_eval(out, out_sim);
end

%% plot trajectories
if PLOT
    
%     if PLOT
    
    PS = peak_sos_plotter(out, out_sim);
    PS.nonneg_zeta();
    PS.obj_plot();
    PS.state_plot();
    PS.v_plot();
    PS.nonneg_traj();
    
    PS.state_plot_3(box_lim);
    if INIT_POINT
        scatter3(C0(1), C0(2), C0(3), 200, 'ok')
    end
    zlim([0,1])
    view(80, 5)
    
%     viscircles(C0', R0, 'color', 'k', 'LineWidth', 3);
    
%     if ~INIT_POINT
%         theta = linspace(0,2*pi, 200);
%         plot(R0*cos(theta)+C0(1), R0*sin(theta)+C0(2), 'color', 'k', 'LineWidth', 3);
%     end
    
    DG.data_plot_3(observed, 0.2);
%     viscircles(C0', R0, 'color', 'k', 'LineWidth', 3);
    
    %observation plot
    
end


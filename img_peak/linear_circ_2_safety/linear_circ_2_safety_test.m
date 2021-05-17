%data driven peak estimation of a linear system

%break up the sections here into functions

PROBLEM = 0;
SOLVE = 0;
SAMPLE = 1;
PLOT = 1;

if PROBLEM
rng(33, 'twister')
%% generate samples
A_true = [-1 4; -1 -0.3];
% A_true = [-1 1; -1 -0.3];
f_true = @(t, x) A_true*x;

Nsample = 50;
% Nsample = 40;
% Nsample = 30;
% Nsample = 20;
% Nsample = 4;
box_lim = 2;
Tmax = 5;
% epsilon = 2;
epsilon = 2.5;
sample = struct('t', Tmax, 'x', @() box_lim*(2*rand(2,1)-1));

% [observed] = corrupt_observations(Nsample,sample, f_true, epsilon);

%% generate model
t = sdpvar(1, 1);
x = sdpvar(2, 1);

DG = data_generator(sample);

observed = DG.corrupt_observations(Nsample, f_true, epsilon);
[model, W] = DG.reduced_model(observed, x, 1, 1);

[w_handle, box]= DG.make_sampler(W);

end
 
%% Solve SOS program
if SOLVE
    
    %start at a single point

%     C0 = [0; 0.5];
    C0 = [-1; 0];
    R0 = 0.2;
    INIT_POINT = 0;
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
    lsupp.X = struct('ineq', 2*box_lim^2 - sum(x.^2), 'eq', []);
    % lsupp = lsupp.set_box(3);
    lsupp.X_init = X0;
    lsupp.f0 = model.f0;
    lsupp.fw = model.fw;
    lsupp.W = W;
    lsupp.Tmax = Tmax;

    lsupp.verbose = 1;

%     objective = x(1);

    Ru = 0.3;
    Cu = [0; -0.5];
    c1f = Ru^2 - sum((x-Cu).^2);
%     c2f = -diff(x-Cu);

    theta_c = 5*pi/4;
    w_c = [cos(theta_c); sin(theta_c)];
    c2f = w_c(1)*(x(1) - Cu(1)) + w_c(2) * (x(2) - Cu(2)); 

    objective = [c1f; c2f];
    
    %% start up tester
    PM = peak_sos(lsupp, objective);

%     order = 4;
    order = 3;
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

    out_sim = traj_eval(out, out_sim);

    sample_time = toc;
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
    
    PS.state_plot_2(box_lim);
    if INIT_POINT
        scatter(C0(1), C0(2), 200, 'k')
    else
        theta = linspace(0,2*pi, 200);
        plot(R0*cos(theta)+C0(1), R0*sin(theta)+C0(2), 'color', 'k', 'LineWidth', 3);
    end
    
    %plot the unsafe set
    theta_half_range = linspace(theta_c-pi/2, theta_c + pi/2, 200);
    circ_half = [cos(theta_half_range); sin(theta_half_range)];
    Xu = Cu + circ_half* Ru;
    patch(Xu(1, :), Xu(2, :), 'r', 'Linewidth', 3, 'EdgeColor', 'none')
    
    %observation plot    
    DG.data_plot_2(observed);

    
    
end


%plot system trajectories for quadratic rotating system

PROBLEM = 1;
SOLVE = 1;
SAMPLE_TRUE = 0;
SAMPLE = 1;
PLOT = 1;

if PROBLEM
rng(35, 'twister')

f_true = @(t, x) [-0.3*x(1) + 0.8*x(2) + 0.1*x(1)*x(2); -0.75*x(1) - 0.3*x(2)];
beta_x1 = [0; -0.3; 0.8; 0; 0.1; 0];
beta_x2 = [0; -0.75; -0.3; 0; 0; 0];
beta_true = [beta_x1; beta_x2];
C0 = [-1.5; 0];
R0 = 0.4;

Nsample = 40;
box_lim = 2;
Tmax = 10;
epsilon = 0.5;

sample = struct('t', Tmax, 'x', @() box_lim*(2*rand(2,1)-1));


%% generate model
x = sdpvar(2, 1);

DG = data_generator(sample);

observed = DG.corrupt_observations(Nsample, f_true, epsilon);
% [model, W] = DG.reduced_model(observed, x, 1, 1);
% model = DG.poly_model(vars, 3);
mlist = monolist(x,2);
model = struct('f0', [0;0], 'fw', kron(eye(2), mlist'));

% beta_f = model.f0 + model.fw;
W = DG.data_cons(model, x, observed);
[model_cheb,W_cheb] = DG.center_cheb(model, W);
W_red = DG.reduce_constraints(W_cheb);

model = model_cheb;
W = W_red;
[w_handle, box]= DG.make_sampler(W_cheb);
end

INIT_POINT = 1;

if SOLVE
    if INIT_POINT
        X0 = C0;
    else
        X0 = struct('ineq', R0^2 - sum((x-C0).^2), 'eq', []);
    end

    objective = -x(2);

    box_lim = 2;
    lsupp = loc_sos_options();    
    lsupp.x = x;
    lsupp.TIME_INDEP = 1;
    lsupp.DISCRETE_TIME = 1;    
    % lsupp = lsupp.set_box(box_lim*sqrt(2));
    lsupp = lsupp.set_box(box_lim);
    % lsupp.X = struct('ineq', 2*box_lim^2 - sum(x.^2), 'eq', []);
    % lsupp = lsupp.set_box(3);
    lsupp.X_init = X0;
    lsupp.f0 = model.f0;
    lsupp.fw = model.fw;
    lsupp.W = W;
    lsupp.Tmax = Tmax;
    

    lsupp.verbose = 1;

    %% start up tester
    PM = peak_sos(lsupp, objective);


        %point initial set
    % order = 1; %sqrt(8)
    % order = 2; %sqrt(8)
    order = 3; %1.1992365
    % order = 4;  

    d = 2*order;

    % [prog]= PM.make_program(d);
    % out = PM.solve_program(prog)
    out = PM.run(order);
    
else
    out = [];
end

if SAMPLE_TRUE
    s_opt = sampler_options;
    s_opt.Tmax = Tmax;

    if INIT_POINT
        s_opt.sample.x = @() C0;
    else
        s_opt.sample.x = @() R0*ball_sample(1,2)'+C0;
    end

    s_opt.parallel = 0;
    Nsample_traj = 400;

    dynamics_true = struct('f', {@(t, x, w, d, b) f_true(t, x)}, 'discrete', 1, 'cost', @(x) -x(2));

    tic
    out_sim_true = sampler(dynamics_true, Nsample_traj, s_opt);

    sample_time = toc;

end


if SAMPLE
    
    s_opt = sampler_options;
    s_opt.Tmax = Tmax;
    s_opt.Nd = size(model.fw, 2);

    if INIT_POINT
        s_opt.sample.x = @() C0;
    else
        s_opt.sample.x = @() R0*ball_sample(1,2)'+C0;
    end

    s_opt.sample.d = w_handle;
    s_opt.parallel = 0;

    Nsample_traj = 100;

    p0 = polyval_func(model.f0, x); 
    pw = polyval_func(model.fw, x); 

    dynamics = struct('f', {@(t, x, w, d, b) p0(x) + pw(x)*d}, 'discrete', 1, 'cost', @(x) -x(2));

    tic
    out_sim = sampler(dynamics, Nsample_traj, s_opt);

    sample_time = toc;



end

if PLOT
    PS_true = peak_sos_plotter(out, out_sim_true);
    PS = peak_sos_plotter(out, out_sim);
    PS.state_plot_2_discrete(3);

    if INIT_POINT
        scatter(C0(1), C0(2), 200, 'k')
    else
        theta = linspace(0,2*pi, 200);
        plot(R0*cos(theta)+C0(1), R0*sin(theta)+C0(2), 'color', 'k', 'LineWidth', 3);
    end

        DG.data_plot_2_discrete(observed, epsilon);


end
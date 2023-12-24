%data driven peak estimation of a linear system

%break up the sections here into functions

PROBLEM = 1;
SOLVE = 1;
SAMPLE = 0;
PLOT = 0;

if PROBLEM
rng(33, 'twister')
%% generate samples
% A_true = [-1 1; -1 -0.3];
f_true = @(t, x) [x(2); -x(1) + (1/3).* x(1).^3 - x(2)];


% Nsample = 100;
% Nsample = 50;
Nsample = 40;
% Nsample = 30;
% Nsample = 20;
% Nsample = 4;
box_lim = 2;
Tmax = 5;
% epsilon = 2;
% epsilon = [0; 2.5];'
% epsilon = [0; 1];
epsilon = [0; 0.5];
sample = struct('t', Tmax, 'x', @() box_lim*(2*rand(2,1)-1));



%% generate model
t = sdpvar(1, 1);
x = sdpvar(2, 1);
y = sdpvar(2, 1);

DG = data_generator(sample);

observed = DG.corrupt_observations(Nsample, f_true, epsilon);
% [model, W] = DG.reduced_model(observed, x, 1, 1);
% model = DG.poly_model(vars, 3);
mlist = monolist(x, 3);
DATA_DRIVEN = 1;

%data-driven model
if DATA_DRIVEN
    model = struct('f0', [x(2);0], 'fw', [zeros(1, length(mlist)); mlist']);
    W = DG.data_cons(model, x, observed);
    %order1: 1.697859227531732e-05
    %order2: 0.193576151666294
    %order3: 0.200364187374125
    %order4: 0.200965666044131
    %order5: 0.201360421158720 (1629.92 sec)
    
else
% %simple single source of uncertainty
    model = struct('f0', [x(2); -x(1) + (1/3).* x(1).^3 - x(2)], 'fw', [0; -x(1)]);
    W = struct('A', [1; -1], 'b', [1; -1]);
    %order4: 2.2272e-01
end

[model_cheb,W_cheb] = DG.center_cheb(model, W);
W_red = DG.reduce_constraints(W_cheb);


model = model_cheb;
W = W_red;
[w_handle, box]= DG.make_sampler(W);

end
 
%% Solve SOS program
if SOLVE
    
    %start at a single point

%     C0 = [1.5; 0];
    R0 = 0.2;
% C0 = [1.5; 0]; %order  4: 8.6664e-02 
C0 = [1; 0];
% C0 = [1.5; 0.5];
% R0 = 0.4;
    INIT_POINT = 1;
    if INIT_POINT
        X0 = C0;
    else 
        
        X0 = struct('ineq', R0^2 - sum((x-C0).^2), 'eq', []);
    end
    
    %unsafe set
    
%     Ru = 0.3;
%     Cu = [0; -0.5];
% Cu = [0; -0.7];
% Cu = [-0.5; -1];
Cu = [-0.25; -0.7];
Ru = 0.5;
    c1f = Ru^2 - sum((y-Cu).^2);
%     c2f = -diff(x-Cu);

    theta_c = 5*pi/4;
    w_c = [cos(theta_c); sin(theta_c)];
    c2f = w_c(1)*(y(1) - Cu(1)) + w_c(2) * (y(2) - Cu(2)); 

    Xu = struct('ineq', [c1f; c2f], 'eq', []);
    
    
    lsupp = dist_sos_options();
    lsupp.t = t;
    lsupp.TIME_INDEP = 0;
    lsupp.x = x;
    lsupp.y = y;
    lsupp = lsupp.set_box(box_lim);
%     lsupp.X = struct('ineq', 2*box_lim^2 - sum(x.^2), 'eq', []);
    lsupp = lsupp.set_box([-1, 1.25; -1.25, 0.7]); %order 4 7.0404e-04 
    lsupp.X_init = X0;
    lsupp.X_unsafe = Xu;
    lsupp.f0 = model.f0;
    lsupp.fw = model.fw;
    lsupp.Tmax = Tmax;
    lsupp.W = W;

    lsupp.verbose = 1;

%     objective = x(1);

    objective = [c1f; c2f];
    
    %% start up tester
    PM = dist_sos(lsupp);

%     order = 5;
%     order = 4;
order=3; %
%     order = 2; % 
%     order = 1;
    d = 2*order;

    % [prog]= PM.make_program(d);
    % out = PM.solve_program(prog)
    out = PM.run(order);
    dist_close = sqrt(-out.obj);
    
    fprintf('distance cost: %0.4e \n', dist_close);
    
end

%% Sample trajectories
if SAMPLE
    
    s_opt = sampler_options;
    s_opt.mu = 0.2;
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
%     Nsample_traj = 120;
    Nsample_traj = 200;
    
    tic
    out_sim = sampler(out.dynamics, Nsample_traj, s_opt);

%     out_sim = traj_eval(out, out_sim);

    sample_time = toc;
end

%% plot trajectories
if PLOT
    
%     if PLOT
    
    PS = peak_sos_plotter(out, out_sim);
%     PS.nonneg_zeta();
%     PS.obj_plot();
%     PS.state_plot();
%     PS.v_plot();
%     PS.nonneg_traj();
    DG.data_plot_2(observed);
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
    
    Ntheta = 50;
    x_dist = dist_contour(Ntheta, Ru, max(dist_close));
    theta_r = theta_c + pi/2;
    Rot_dist = [cos(theta_r), -sin(theta_r); sin(theta_r), cos(theta_r)];
    
    x_dist = Rot_dist*x_dist + Cu;
    
    plot(x_dist(1, :), x_dist(2, :), 'LineWidth', 3, 'color', 'r');

    
    %observation plot    


    
    
end

function x_dist = dist_contour(Ntheta, R, c)
    %compute a contour at distance c away from the half-circle with N_theta
    %sample points


    theta_q1 = linspace(0, pi/2, Ntheta);
    theta_q2 = linspace(pi/2, pi, Ntheta);
    theta_q34 = linspace(pi, 2*pi, 2*Ntheta);

    %contour level
    

    x_right = [c*cos(theta_q1)+R; c*sin(theta_q1)];
    x_left = [c*cos(theta_q2)-R; c*sin(theta_q2)];
    x_bottom = [(c+R)*cos(theta_q34); (c+R)*sin(theta_q34)];
    x_dist = [x_right, x_left, x_bottom];
end


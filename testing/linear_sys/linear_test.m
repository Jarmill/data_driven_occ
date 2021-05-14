%data driven peak estimation of a linear system

%break up the sections here into functions

PROBLEM = 1;
SOLVE = 0;
SAMPLE = 0;
PLOT = 0;

if PROBLEM
rng(33, 'twister')
%% generate samples
A_true = [-1 3; -1 -0.3];
% A_true = [-1 1; -1 -0.3];
f_true = @(t, x) A_true*x;

% Nsample = 50;
Nsample = 40;
% Nsample = 30;
% Nsample = 20;
% Nsample = 4;
box_lim = 2;
Tmax = 5;
epsilon = 2;
% epsilon = 1;
sample = struct('t', Tmax, 'x', @() box_lim*(2*rand(2,1)-1));

% [observed] = corrupt_observations(Nsample,sample, f_true, epsilon);

%% generate model
t = sdpvar(1, 1);
x = sdpvar(2, 1);

DG = data_generator(sample);

observed = DG.corrupt_observations(Nsample, f_true, epsilon);

model = DG.poly_model(x,1,1);
% model = DG.poly_model(struct('t', t, 'x', x),1,1);

%% create observation constraints Aw <= b
f0_func = polyval_func(f0, [t;x]);
fw_func = polyval_func(fw, [t;x]);

nw = size(fw, 2);   %number of uncertainties
m = nx*Nsample*2;   %number of constraints
A = zeros(m, nw);
b = zeros(m, 1);

counter = 0;
for i = 1:Nsample
    tcurr = observed.t(i);
    xcurr = observed.x(:, i);
    
    xdotcurr = observed.xdot_noise(:, i);
    
    f0_curr = f0_func([tcurr; xcurr]);
    fw_curr = fw_func([tcurr; xcurr]);
    
    %ensure the correct signs over here
    b_pos_curr = epsilon - f0_curr + xdotcurr;
    b_neg_curr = epsilon + f0_curr - xdotcurr;
    
    A(counter+(1:2*nx), :) = [fw_curr; -fw_curr];
    b(counter+(1:2*nx)) = [b_pos_curr; b_neg_curr];
    
    counter = counter + 2*nx;
end

%% Center the polytope

%box
% box = poly_bounding_box(A,b);
% [box_out, box_center, box_half] = box_process(nw, box);
% f0_center = f0 + fw*box_center;
% fw_center = fw*diag(box_half);
% 
% A_scale = A*diag(box_half);
% b_scale = b - A*box_center;

%chebyshev
[c,r] = chebycenter(A,b);
A_scale = A;
b_scale = b - A*c;

f0_center = f0 + fw*c;
fw_center = fw;

%identify redundant constraints
%code from noredund.m by Michael Kleder (2006)
D = A_scale ./ repmat(b_scale,[1 size(A_scale,2)]);
%number of points in convex hull capped at number of constraints
[k, vol] = convhulln(D);
% record which constraints generate points on the convex hull
nr = unique(k(:));
A_scale_orig = A_scale;
b_scale_orig = b_scale;
A_scale=A_scale(nr,:);
b_scale=b_scale(nr);

box_cheb = poly_bounding_box(A_scale,b_scale);
[box_cheb, box_center, box_half] = box_process(nw, box_cheb);

W = struct('A', A_scale, 'b', b_scale);


end
 
%% Solve SOS program
if SOLVE
    
    %start at a single point

    C0 = [0; 0.5];
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
    lsupp = lsupp.set_box(box_lim);
    lsupp.X = struct('ineq', 2*box_lim^2 - sum(x.^2), 'eq', []);
    % lsupp = lsupp.set_box(3);
    lsupp.X_init = X0;
    lsupp.f0 = f0_center;
    lsupp.fw = fw_center;
    lsupp.W = W;
    lsupp.Tmax = Tmax;

    lsupp.verbose = 1;

    objective = x(1);
    
    %% start up tester
    PM = peak_sos(lsupp, objective);

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
        s_opt.sample.x = @() ball_sample(2,1);
    end
    s_opt.sample.d = @() rej_sample_poly(A_scale, b_scale, box_cheb);
    s_opt.Nd = nw;
    
    s_opt.Tmax = lsupp.Tmax;
    s_opt.parallel = 1;
    
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
    PS.cost_plot();
    PS.state_plot();
    PS.nonneg_traj();
    
    PS.state_plot_2();
    
    %observation plot
    figure(2)
    clf
    hold on
    quiver(observed.x(1, :), observed.x(2, :), observed.xdot_true(1, :), observed.xdot_true(2, :))
    quiver(observed.x(1, :), observed.x(2, :), observed.xdot_noise(1, :), observed.xdot_noise(2, :))
    axis square
    xlabel('$x_1$', 'interpreter', 'latex', 'FontSize', PS.FS_axis);
    ylabel('$x_2$', 'interpreter', 'latex', 'FontSize', PS.FS_axis);          
    title(['Noisy Observations with \epsilon=', num2str(epsilon)], 'FontSize', PS.FS_title)

end


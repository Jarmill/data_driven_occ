%data driven peak estimation of a linear system

PROBLEM = 1;
SOLVE = 0;
SAMPLE = 0;
PLOT = 0;

if PROBLEM
rng(33, 'twister')
%% generate samples
f_true = @(t,x) [x(2); (-x(1)-x(2) + (1/3)*x(1)^3)];


% Nsample = 50;
Nsample = 10;
% Nsample = 4;
box_lim = 3;
Tmax = 5;
% epsilon = 2;
epsilon = 0.3;
sample = struct('t', Tmax, 'x', @() box_lim*(2*rand(2,1)-1));

[observed] = corrupt_observations(Nsample,sample, f_true, epsilon);

%% generate model
t = sdpvar(1, 1);
x = sdpvar(2, 1);

dmin = 1;
dmax = 1;
nx = 2;
f0 = [x(2); -x(2) + (1/3)*x(1)^3];
% fw = [0 0; -x(1) 1];
fw = [0; -x(1)];% 
% f0 = zeros(nx, 1);
% fw = kron(eye(nx), monolist(x, dmin, dmax)');
% sdisplay(fw)


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
% D = A_scale ./ repmat(b_scale,[1 size(A_scale,2)]);
% %number of points in convex hull capped at number of constraints
% [k, vol] = convhulln(D);
% % record which constraints generate points on the convex hull
% nr = unique(k(:));
% A_scale_orig = A_scale;
% b_scale_orig = b_scale;
% A_scale=A_scale(nr,:);
% b_scale=b_scale(nr);

box_cheb = poly_bounding_box(A_scale,b_scale);
% [box_cheb, box_center, box_half] = box_process(nw, box_cheb);

W = struct('A', A_scale, 'b', b_scale);


end
 
%% Solve SOS program
if SOLVE
    
    %start at a single point
    X0 = [1.5; 0];
    
    lsupp = loc_sos_options();
    lsupp.t = t;
    lsupp.TIME_INDEP = 0;
    lsupp.x = x;
    lsupp = lsupp.set_box(box_lim);
%     lsupp.X = struct('ineq', 2*box_lim^2 - sum(x.^2), 'eq', []);
%     lsupp = lsupp.set_box(box_lim);;
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
    s_opt.sample.x = @() X0;
    s_opt.sample.d = @() rej_sample_poly(A_scale, b_scale, box_cheb);
    s_opt.Nd = nw;
    
    s_opt.Tmax = lsupp.Tmax;
    s_opt.parallel = 0;
    
    tic
    out_sim = sampler(out.dynamics, Nsample, s_opt);

    out_sim = traj_eval(out, out_sim);

    sample_time = toc;
end

%% plot trajectories
if PLOT
    
%     if PLOT
    
    PS = peak_sos_plotter(out, out_sim);
%     PS.nonneg_zeta();
    PS.state_plot();
    PS.nonneg_traj();
    
end


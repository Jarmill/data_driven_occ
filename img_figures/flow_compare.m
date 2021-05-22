%test out arguments of peak_estimate routine
%clear
mset clear
rng(300, 'twister')

%dynamics and support set
%prajna and rantzer flow
SOLVE = 1;
SAMPLE = 0;
PLOT = 0;

%go back to functions
%and/or figure out how to extract monomials and powers from mpol
mpol('x', 2, 1);
mpol('d', 1, 1);

%dynamics

X = [];


%initial set
%C0 = [1.2; 0];
C0 = [1.5; 0];
R0 = 0.4;

%C0 = [0.8; 0];
%R0 = 0.1;

% X0 = ((x(1)-C0(1))^2 + (x(2)-C0(2))^2 <= R0^2);



rng(50, 'twister')
x0 = C0;
Tmax_sim = 20;

% Nsample = 300;
% Nsample = 150;
% dmax = 0.15;
dmax = 0.5;

%time-dependent
Nsample = 40;

s_opt = sampler_options;
s_opt.sample.x = @() C0;
s_opt.sample.d = @() dmax * (2*rand() -1 );
s_opt.Tmax = Tmax_sim;
s_opt.Nd = 1;
s_opt.parallel = 1;

s_opt.mu = 0.2;

dynamics = struct('Tmax', Tmax_sim, 'discrete', 0);
dynamics.Xval = @(x) all(abs(x) <= 3);
dynamics.event = @(t,x,w) support_event(t, x, dynamics.Xval, ...
                0, Tmax_sim);

dynamics.f = @(t,x,w,d,b)  [x(2); -x(1)*(1+ d*dmax) + (1/3).* x(1).^3 - x(2)];       
            
tic
out_sim_d = sampler(dynamics, Nsample, s_opt);
sample_time = toc;

% time independent
%time-dependent
Nsample = 20;

s_opt = sampler_options;
s_opt.sample.x = @() C0;
s_opt.sample.w = @() dmax * (2*rand() -1 );
s_opt.Tmax = Tmax_sim;
s_opt.Nw = 1;
s_opt.parallel = 1;

s_opt.mu = 0.2;

dynamics = struct('Tmax', Tmax_sim, 'discrete', 0);
dynamics.Xval = @(x) all(abs(x) <= 3);
dynamics.event = @(t,x,w) support_event(t, x, dynamics.Xval, ...
                0, Tmax_sim);

dynamics.f = @(t,x,w,d,b)  [x(2); -x(1)*(1+ w*dmax) + (1/3).* x(1).^3 - x(2)];       
            
tic
out_sim_w = sampler(dynamics, Nsample, s_opt);
sample_time = toc;



%% plot the trajectories
figure(1)
clf
tiledlayout(1,2)
ax1 = nexttile;
    hold on
    for i = 1:Nsample
        if i == 1
            plot(out_sim_w{i}.x(:, 1), out_sim_w{i}.x(:, 2), 'c', 'DisplayName', 'Trajectories');
        else
            plot(out_sim_w{i}.x(:, 1), out_sim_w{i}.x(:, 2), 'c', 'HandleVisibility', 'Off');
        end
    end
    
    scatter(C0(1), C0(2), 200, 'ok')
           
    xlim([-0.3, 1.75])
    ylim([-1,0.5])
    
%     axis square
pbaspect([diff(xlim), diff(ylim), 1])
    xlabel('x_1')
    ylabel('x_2')
%     axis off
title('Time-Independent Uncertainty', 'fontsize', 17)    
hold off
axis off
ax2 = nexttile;
    hold on
    for i = 1:Nsample
        if i == 1
            plot(out_sim_d{i}.x(:, 1), out_sim_d{i}.x(:, 2), 'c', 'DisplayName', 'Trajectories');
        else
            plot(out_sim_d{i}.x(:, 1), out_sim_d{i}.x(:, 2), 'c', 'HandleVisibility', 'Off');
        end
    end
    
    scatter(C0(1), C0(2), 200, 'ok')
           
%     xlim([-0.5, 1.75])
%     ylim([-1, 0.])
    
    axis square
    
    xlabel('x_1')
    ylabel('x_2')
%     axis off
title('Time-Dependent Uncertainty', 'fontsize', 17)
    
hold off
%visualize

% linkaxes([ax1, ax2])
xlim([-0.3, 1.75])
    ylim([-1,0.5])
    
pbaspect([diff(xlim), diff(ylim), 1])
axis off
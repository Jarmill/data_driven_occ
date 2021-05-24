%visualization of corrupted samples
% rng(35, 'twister')
rng(10, 'twister')
%% generate samples
% A_true = [-1 1; -1 -0.3];
f_true = @(t, x) [x(2); -x(1) + (1/3).* x(1).^3 - x(2)];

Nsample = 40;
box_lim = 2;
Tmax = 5;
% epsilon = 2;
% epsilon = [0; 2.5];
% epsilon = [0; 1];
% epsilon = [0; 0.5];
% epsilon = 0.5;
epsilon = 1.5;

% [-0.3, 1.75])
%     ylim([-1,0.5])
sample = struct('t', Tmax, 'x', @() box_lim*(2*rand(2,1)-1));

%% generate model
t = sdpvar(1, 1);
x = sdpvar(2, 1);

DG = data_generator(sample);

observed = DG.corrupt_observations(Nsample, f_true, epsilon);
model = struct('f0', [x(2); (1/3).* x(1).^3 - x(2)], 'fw', [0; -x(1)]);

W = DG.data_cons(model, x, observed);
[model_cheb,W_cheb] = DG.center_cheb(model, W);
W_red = DG.reduce_constraints(W_cheb);

model = model_cheb;
W = W_red;
[w_handle, box]= DG.make_sampler(W_cheb);


rng(10, 'twister')
observed2 = DG.corrupt_observations(Nsample, f_true, [0; epsilon]);
% observed2 = DG.corrupt_observations(Nsample, f_true, epsilon);


%% plot field
scale = 0.25;
figure(2)
tiledlayout(1,2)
ax1 = nexttile;

hold on
quiver(observed.x(1, :), observed.x(2, :), scale*observed.xdot_true(1, :), scale*observed.xdot_true(2, :), 'autoscale','off')
quiver(observed.x(1, :), observed.x(2, :), scale*observed.xdot_noise(1, :), scale*observed.xdot_noise(2, :), 'autoscale','off')
axis square
legend({'Ground Truth', 'Noisy Data'}, 'FontSize', 10, 'location', 'northwest')
xlabel('$x_1$', 'interpreter', 'latex', 'FontSize', 12);
ylabel('$x_2$', 'interpreter', 'latex', 'FontSize', 12);          
title(['40 observations with \epsilon=', num2str(epsilon)], 'FontSize', 16)


ax2 = nexttile;
hold on
quiver(observed2.x(1, :), observed2.x(2, :), scale*observed2.xdot_true(1, :),  scale*observed2.xdot_true(2, :), 'autoscale','off')
quiver(observed2.x(1, :), observed2.x(2, :), scale*observed2.xdot_noise(1, :), scale*observed2.xdot_noise(2, :), 'autoscale','off')
axis square
legend({'Ground Truth', 'Noisy Data'}, 'FontSize', 10, 'location', 'northwest')
xlabel('$x_1$', 'interpreter', 'latex', 'FontSize', 12);
ylabel('$x_2$', 'interpreter', 'latex', 'FontSize', 12);          
title(['40 observations with \epsilon=[0, ', num2str(epsilon), ']'], 'FontSize', 16)

linkaxes([ax1, ax2])
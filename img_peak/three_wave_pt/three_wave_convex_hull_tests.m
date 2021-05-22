%data driven peak estimation of a linear system

%break up the sections here into functions

PROBLEM = 1;
PLOT = 0;


%sample data only from initial set
C0 = [1; 1; 1];
R0 = 0.4;
x0_sample_ball = @() R0*ball_sample(1,3)'+C0;

INIT_SAMPLE_ONLY = 0;

rng(33, 'twister')
%% generate samples
a = 1;
b = 1;
G0 = 2;
f_true = @(t,x) [a*x(1) + b*x(2) + x(3) - 2*x(2)^2;
    a*x(2) - b*x(1) + 2*x(1)*x(2);
    -G0*x(3) - 2*x(1)*x(3)];

Nsample = 150;
% Nsample = 100;
% Nsample = 50;
% Nsample = 40;
% Nsample = 30;
% Nsample = 20;
% Nsample = 10;
% Nsample = 4;
box_lim = [-4, 0.5, 0; 3, 3.6, 4]';
Tmax = 5;
epsilon = 1;
% epsilon = 2;
% epsilon = 2.5;
sample = struct('t', Tmax, 'x', @() box_lim(:, 1) + diff(box_lim, 1, 2).*rand(3, 1));
if INIT_SAMPLE_ONLY
    sample.x = x0_sample_ball;
end

% [observed] = corrupt_observations(Nsample,sample, f_true, epsilon);

%% generate model
t = sdpvar(1, 1);
x = sdpvar(3, 1);

DG = data_generator(sample);


observed = DG.corrupt_observations(Nsample, f_true, epsilon);
% model = DG.poly_model(x, 0,2);
mlist23 = monolist(x, 2, 0);
mlist2 = monolist(x(1:2), 2, 0);
mlist3 = monolist(x([1,3]), 2, 0);
model.f0 = [0;0;0];
model.fw = blkdiag(mlist23, mlist2, mlist3)';
W_orig = DG.data_cons(model, x, observed);
W = W_orig;



[model_cheb,W_cheb] = DG.center_cheb(model, W);
box_cheb = poly_bounding_box(W_cheb.A, W_cheb.b);


if PLOT
    DG.data_plot_3(observed);
end


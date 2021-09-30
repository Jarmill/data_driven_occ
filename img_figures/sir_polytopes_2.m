%data driven reachable set estimation of a linear system

%break up the sections here into functions

PROBLEM = 1;
% SOLVE = 1;
% SAMPLE = 1;
% PLOT = 1;

if PROBLEM
rng(33, 'twister')
%% generate samples
beta_true = 0.4;
gamma_true = 0.1;
f_true = @(t,x) [-beta_true*x(1)*x(2); beta_true*x(1)*x(2)- gamma_true*x(2)];
% A_true = [-1 4; -1 -0.3];
% A_true = [-1 1; -1 -0.3];
% f_true = @(t, x) A_true*x;
Imax0 = 1;
% Imax0 = 0.2;
si_sampler = @() rej_sample_poly([-1 0; 0 -1; 1 1; 0 1], [0;0;1;Imax0], eye(2));


Nsample = 100;
% Nsample = 50;
% Nsample = 40;
% Nsample = 30;
% Nsample = 20;
% Nsample = 10;
% box_lim = 1.25;
box_lim = 1;
Tmax = 40;

% epsilon = 0.05;
epsilon = 1e-1;
% epsilon = 1;
% epsilon = 2.5;
sample = struct('t', Tmax, 'x', si_sampler);

% [observed] = corrupt_observations(Nsample,sample, f_true, epsilon);

%% generate model
t = sdpvar(1, 1);
x = sdpvar(2, 1);

DG = data_generator(sample);

observed = DG.corrupt_observations(Nsample, f_true, epsilon);
% [model, W] = DG.reduced_model(observed, x, 1, 2);

model = struct('f0', [0;0], 'fw', [-x(1)*x(2) 0; x(1)*x(2) -x(2)]);
W_orig = DG.data_cons(model, x, observed);
W = W_orig;
W.A = W.A(:, [2,1]);


[model_cheb,W_cheb] = DG.center_cheb(model, W);
W_red = DG.reduce_constraints(W_cheb);

W_red_2 = DG.reduce_constraints_linprog(W_cheb);


[length(W_red.b), length(W_red_2.b)]

[w_handle, box]= DG.make_sampler(W_red);



% box
% size(W_red.b)
end

%% start the plots
DG.data_plot_2(observed);
legend('location', 'northeast', 'fontsize', 10)
% title('
title([num2str(Nsample), ' SIR Observations'], 'FontSize', 16)
xlabel('$S$', 'interpreter', 'latex', 'FontSize', 16)
ylabel('$I$', 'interpreter', 'latex', 'FontSize', 16)
xlim([-0.05, 1])
ylim([-0.05, 1])
axis square
%% chebyshev center figure
figure(3)
clf
hold on
[c,r] = chebycenter(W.A,W.b);
V = lcon2vert(W.A, W.b);

K = convhull(V);
plot(V(K, 1), V(K, 2), 'k', 'LineWidth', 2)

theta = linspace(0,2*pi, 200);
patch(r*cos(theta)+c(1), r*sin(theta)+c(2), 0.7*[1,1,1]);
scatter(c(1), c(2), 200, '*k')
hold off

ylabel('$\beta$', 'interpreter', 'latex', 'fontsize', 16)
xlabel('$\gamma$', 'interpreter', 'latex', 'fontsize', 16)
% title('Chebyshev Center of SIR', )
ylim([0.37, 0.423])
pbaspect([diff(xlim), diff(ylim), 1])

%% redundant constraints figures
figure(4)
%chebyshev center figure with boundary lines
ax2 = subplot(2,1,1);
% clf
hold on
% xrange = 0.10*[-1,1];
% yrange = 0.05*[-1,1];
xrange = 0.1*[-1,1];
yrange = 0.05*[-1,1];

for i = 1:length(W_cheb.b)    
    p_curr = line_range([W_cheb.A(i, :), -W_cheb.b(i)], xrange, yrange);
    if ~isempty(p_curr)
        if any(i==KD)
            plot([p_curr{1}(1), p_curr{2}(1)], [p_curr{1}(2), p_curr{2}(2)], ':k', 'LineWidth', 3)
        else
            plot([p_curr{1}(1), p_curr{2}(1)], [p_curr{1}(2), p_curr{2}(2)], ':k')
        end
    end
end
plot(V(K, 1)-c(1), V(K, 2)-c(2), 'k', 'LineWidth', 2)
hold off

axis off
title('Primal Polytope', 'Fontsize', 16)


%% compare against the ground set
figure(10)
hold on
plot(V(K, 1)-c(1), V(K, 2)-c(2), 'k', 'LineWidth', 2)
plot([1, 1, -1, -1, 1]*epsilon, [1, -1, -1, 1, 1]*4*epsilon, 'r')
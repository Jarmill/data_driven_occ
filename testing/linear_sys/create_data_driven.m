PLOT = 1;

rng(341, 'twister')
mpol('t', 1,1);
mpol('x', 2, 1);
mpol('w', 4, 1);

%start out with an uncertain linear system

%true dynamics
% A_true = [-0.6 0.5; 0.45 -0.4];
A_true = [-1 3; -1 -0.3];

f_true = @(x) A_true*x;


%generate data
% Nsample = 400;
Nsample = 40;
box_lim = 3;

epsilon = 2;

x_sample = (2*rand(2, Nsample)-1)*box_lim;
xtest = x_sample(:, 1);
noise_sample = (2*rand(2, Nsample)-1)*epsilon;

xdot_true = f_true(x_sample);
xdot_noise = xdot_true + noise_sample;

%data-driven bounds
A_w = kron(x_sample', eye(2));
b_w = reshape(xdot_noise, [], 1);
b_w_true = reshape(xdot_true, [], 1);

w_est_true = A_w \ b_w_true;
w_est = A_w \ b_w;

%chebyshev center
% -epsilon <= (A_w*w - b_w) <= epsilon;

A_w_signed = kron([1; -1], A_w);
b_w_signed = [epsilon + b_w; epsilon - b_w];
%A_w_signed w <= b_w_signed

%c: chebyshev center of polytope
%r: radius of largest hypersphere enclosed by polytope
[c, r] = chebycenter(A_w_signed, b_w_signed);



%% manipulate center

% b_cheb_pos = b_w + epsilon - A_w*c;
% b_cheb_neg = b_w - epsilon - A_w*c;
% b_cheb_stack = [b_cheb_pos b_cheb_neg];
b_w_cheb = [epsilon + b_w - A_w*c; epsilon - b_w + A_w*c];

[V_cheb]=lcon2vert(A_w_signed ,b_w_cheb,[],[]);

box_cheb = poly_bounding_box(A_w_signed ,b_w_cheb);
% b_w_cheb = b_w_signed - kron([1; -1], A_w*c);





sys_true = ss(A_true, [1; 0], [0 1], 0);

X0 = [0; 0.5];
[tp,xp] = ode45(@(t,x)f_true(x), [0, Tmax], X0);


if PLOT
figure(1)
clf
hold on
impulse(sys_true)

figure(2)
clf
hold on
quiver(x_sample(1, :), x_sample(2, :), xdot_true(1, :), xdot_true(2, :))
quiver(x_sample(1, :), x_sample(2, :), xdot_noise(1, :), xdot_noise(2, :))
axis square
end
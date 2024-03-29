%the valid parameter polytope is not centrally symmetric
%this is a central assumption of the measure code which is false
%will need to find a way around this


PLOT = 1;

rng(341, 'twister')
mpol('t', 1,1);
mpol('x', 2, 1);
mpol('w', 2, 1);

%start out with an uncertain linear system

%true dynamics
% A_true = [-0.6 0.5; 0.45 -0.4];

% A_true = [-1 3; -1 -0.3];



f_true = @(x) [x(2); (-x(1)-x(2) + (1/3)*x(1)^3)];

f_model = [0;0]*x(2);
f_model(1) = x(2);
% f_model(2) = -x(1)*(1+w(1))-x(2) + (1/3)*x(1)^3 + w(2);
f_model(2) = -x(1)*w(1)-x(2) + (1/3)*x(1)^3 + w(2);
Nx = length(f_model); %=length(x);
Nw = length(w);

f_box = zeros(length(f_model), Nw+1) * w(1);
f0 = subs(f_model, w, zeros(Nw, 1));                    
f_box(:, 1) = f0;

%each input channel at a time
I = eye(Nw);

for k = 1:Nw
    f_box(:, k+1) = subs(f_model, w, I(:, k)) - f0;                        
end


%generate data
% Nsample = 400;
Nsample = 100;
box_lim = 3;

epsilon = 2;

x_sample = (2*rand(2, Nsample)-1)*box_lim;
xtest = x_sample(:, 1);
noise_sample = (2*rand(2, Nsample)-1)*epsilon;

xdot_true = zeros(Nx, Nsample);
for i = 1:Nsample
    xdot_true(:, i) = f_true(x_sample(:, i));
end
% xdot_true = f_true(x_sample);
xdot_noise = xdot_true + noise_sample;

f_box_eval = zeros(length(f_model)*Nsample, Nw+1);
counter = 0;
for i = 1:Nsample
    for j = 1:(Nw+1)
        f_box_eval(counter + (1:Nx), j)= eval(f_box(:, j), x, x_sample(:, i));
    end
    counter = counter + Nx;
end

%data-driven bounds

A_w = f_box_eval(:, 2:end);
b_w = reshape(xdot_noise, [], 1) - f_box_eval(:, 1);
b_w_true = reshape(xdot_true, [], 1) - f_box_eval(:, 1);


% A_w = kron(x_sample', eye(2));
% b_w = reshape(xdot_noise, [], 1);
% b_w_true = reshape(xdot_true, [], 1);
% 
w_est_true = A_w \ b_w_true;
w_est = A_w \ b_w;

%chebyshev center
% -epsilon <= (A_w*w - b_w) <= epsilon;

A_w_signed = kron([1; -1], A_w);
b_w_signed = [epsilon + b_w; epsilon - b_w];
%A_w_signed w <= b_w_signed

%find vertices
[V,nr,nre]=lcon2vert(A_w_signed ,b_w_signed,[],[]);



%c: chebyshev center of polytope
%r: radius of largest hypersphere enclosed by polytope
[c, r] = chebycenter(A_w_signed, b_w_signed);

%% manipulate center



b_cheb_pos = b_w + epsilon - A_w*c;
b_cheb_neg = epsilon - b_w + A_w*c;
b_cheb_stack = [b_cheb_pos b_cheb_neg];

[V_cheb]=lcon2vert(A_w_signed ,[b_cheb_pos; b_cheb_neg],[],[]);
k_cheb = convhull(V_cheb);
% b_w_cheb = b_w_signed - kron([1; -1], A_w*c);





% sys_true = ss(A_true, [1; 0], [0 1], 0);




if PLOT
    
    FS = 16
figure(2)
clf
hold on
quiver(x_sample(1, :), x_sample(2, :), xdot_true(1, :), xdot_true(2, :))
quiver(x_sample(1, :), x_sample(2, :), xdot_noise(1, :), xdot_noise(2, :))
axis square
xlabel('x_1')
ylabel('y_1')
title('Recorded data from Flow system', 'FontSize', FS)
legend({'True Flow', 'Corrupted Flow'}, 'location', 'northwest')


figure(3)
clf
hold on
% scatter(V(:, 1), V(:, 2))
plot(V_cheb(k_cheb, 1), V_cheb(k_cheb, 2), 'LineWidth', 3)
scatter(V_cheb(:, 1), V_cheb(:, 2), 200, 'ok', 'LineWidth', 2)
title('Chebyshev-Centered Data-Driven Polytope', 'FontSize', FS)
xlabel('w_1')
ylabel('w_2')
grid on
end


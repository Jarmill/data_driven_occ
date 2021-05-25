% 
% a = 1;
% b = 1;
% G0 = 3;

a = 0.5;
b = 0.1;
G0 = 6;
f_true = @(t,x) [a*x(1) + b*x(2) + x(3) - 2*x(2)^2;
    a*x(2) - b*x(1) + 2*x(1)*x(2);
    -G0*x(3) - 2*x(1)*x(3)];

% [tt, xx] 
C0 = [0; 0; 0.3];
Tmax = 100;
[tt, xt] = ode45(f_true, [0, Tmax], C0);
figure(1)
clf
hold on
plot3(xt(:, 1), xt(:, 2), xt(:, 3))
scatter3(xt(1, 1), xt(1, 2), xt(1, 3), 100, 'k')
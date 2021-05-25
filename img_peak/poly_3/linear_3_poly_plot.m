% 
% a = 1;
% b = 1;
% G0 = 3;

a = 0.5;
b = 0.1;
G0 = 6;

% A_true = [-1 4 1; -1 0 1; 0 1 -2];
% % B_true = eye(3);
% B_true = [3 0 -2;
%           3 -3 -1;
%           3 -1 0];
% 
% f_true = @(t,x) A_true*x - B_true*(x.^3);
% 
% B_true = [3 0 -2;
%           3 -3 -1;
%           3 -1 0];
% B_true = eye(3);
A_true = [-1 1 1; -1 0 -1; 0 1 -2];
B_true = [-1 0 -1;
          0 1 1;
          1 1 0]/2;

f_true = @(t,x) A_true*x - B_true*(4*x.^3 - 3*x);




% [tt, xx] 
C0 = [-0.5; 0; 0];
% Tmax = 1000;
Tmax = 5;
[tt, xt] = ode45(f_true, [0, Tmax], C0);
figure(1)
clf
hold on
[min(xt, [], 1);
max(xt, [], 1)]
plot3(xt(:, 1), xt(:, 2), xt(:, 3))
scatter3(xt(1, 1), xt(1, 2), xt(1, 3), 100, 'k')


FS_axis = 12;
xlabel('$x_1$', 'interpreter', 'latex', 'FontSize', FS_axis);
ylabel('$x_2$', 'interpreter', 'latex', 'FontSize', FS_axis);
zlabel('$x_3$', 'interpreter', 'latex', 'FontSize', FS_axis);

view(3)
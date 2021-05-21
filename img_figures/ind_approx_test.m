X = [0, 1];
% 
% XT = [0.1, 0.5;
%       0.8, 0.9;
%       0.95, 1];


XT = [0.1, 0.5;
      0.8, 0.9];

% d = 30;

% dlist = [10, 100];

dlist = [6, 20, 120];


out = cell(length(dlist),1);
for i = 1:length(dlist)
    out{i} = ind_approx(dlist(i), X, XT);
end


%% plot the indicator approximators

colors = linspecer(length(dlist));
xsample = linspace(0,1,200);


figure(5)
clf
hold on
Xx = [0; 0.1;0.1;0.5;0.5;0.8;0.8;0.9;0.9; 1];
Xy = [0;  0;   1;  1;  0;  0;  1;  1;  0; 0]; 
plot(Xx, Xy, 'k', 'LineWidth', 6)
xlabel('$x$', 'interpreter', 'latex', 'fontsize', 14)
ylabel('$\omega(x)$', 'interpreter', 'latex', 'fontsize', 14)
for i = 1:length(dlist)
    wsample = arrayfun(out{i}.w_eval, xsample);
    plot(xsample, wsample, 'color', colors(i, :), 'LineWidth',3)
end
hold off
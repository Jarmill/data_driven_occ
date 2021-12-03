%plot the uncertain flow system under switching uncertainty

load('flow_rand_sample.mat')
Nsample = length(osm);

figure(1)
clf
hold on 
for i = 1:length(osm)
    plot(osm{i}.x(:, 1), osm{i}.x(:, 2), 'c');
%     plot(xlim, -sol.obj_rec*[1, 1], '--r', 'LineWidth', 3)
    plot(xlim, -[0.786230878194936]*[1, 1], 'r', 'LineWidth', 3)
end

    theta = linspace(0, 2*pi, 100);
    circ = [cos(theta); sin(theta)];
 
    X0 = C0 + circ*R0;
    plot(X0(1, :), X0(2, :), 'k', 'LineWidth', 3)

xlim([-0.5, 3])
ylim([-1, 1.5])
pbaspect([diff(xlim), diff(ylim), 1])
axis off
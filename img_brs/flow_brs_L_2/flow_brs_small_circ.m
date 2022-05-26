%data driven region of attraction maximizing control of a flow system

PROBLEM = 1;
SOLVE = 1;
PLOT = 1;

if PROBLEM
f_true = @(t, x) [x(2); -x(1) + (1/3).* x(1).^3 - x(2)];
% f_true = @(t,x) [-x(2); x(1)];
box_lim = [-2,2;-2,2];
% box_lim =  [-0.3, 1.75;-1,0.5];
Tmax = 5;

%% generate model
t = sdpvar(1, 1);
x = sdpvar(2, 1);

model = struct;
model.f0 = f_true(t,x);
model.fw = eye(2);


u_max = 0.1;
u_slant = 0.4;
Asym = [eye(2); -1,1];
bsym = [u_max; u_max; u_slant];

W = struct;
W.A = [Asym; -Asym];
W.b = [bsym; bsym];

end

if SOLVE
    CT = [0.5; 0];
%     CT = [0; 0];
    RT = 0.1;
    XT = struct('ineq', RT^2 - sum((x-CT).^2), 'eq', []);
    
    
    %define options
        
    lsupp = loc_sos_options();
    lsupp.t = t;
    lsupp.TIME_INDEP = 0;
    lsupp.x = x;
    lsupp = lsupp.set_box(box_lim);
    lsupp.X_term = XT;
    lsupp.f0 = model.f0;
    lsupp.fw = model.fw;
    lsupp.W = W;
    lsupp.Tmax = Tmax;
    
    lsupp.verbose = 1;
    
    box_supp = box_process(2, box_lim);
    mom_handle = @(d) LebesgueBoxMom(d, box_supp', 1);

    
    %% start up tester
    PM = brs_sos(lsupp, mom_handle);

%     order = 5;
%     order = 2;
%     order = 3;
    order = 5;
    d = 2*order;

    % [prog]= PM.make_program(d);
    % out = PM.solve_program(prog)
    out = PM.run(order);
    
    fprintf('objective: %0.3f \n', out.obj)
end

if PLOT
    
    figure(11)
    fsurf(@(x1,x2) out.func.w([zeros(size(x1)); x1; x2]), reshape(box_lim', 1, []), 'r',  'LineWidth', 2)
    
    
    figure(10)
    clf
    theta = linspace(0, 2*pi, 201);
    hold on
    circ = [cos(theta); sin(theta)]*RT + CT;

    xlim(box_lim(1, :));
    ylim(box_lim(2, :));

    %plot the roa
    fcontour(@(x1,x2) out.func.w([zeros(size(x1)); x1; x2]), reshape(box_lim', 1, []), ...
        'r', 'LevelList', 1, 'LineWidth', 2, 'Fill','on')
    colormap([1 0 0; 1 1 1]);
    plot(circ(1, :), circ(2, :), 'k', 'LineWidth', 2)
    axis square
    xlabel('$x_1$', 'interpreter', 'latex', 'Fontsize', 14)
    ylabel('$x_2$', 'interpreter', 'latex', 'Fontsize', 14)
    title(sprintf('Order %d Region of Attraction', order), 'fontsize', 16)

%     fcontour(@(x1,x2) out.func.w([zeros(size(x1)); x1; x2])-1, reshape(box_lim', 1, []))
    
end
%set: -umax <= u <= umax,  -umax <= u1 - u2 <= umax





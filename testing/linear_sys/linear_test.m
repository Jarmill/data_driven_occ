%data driven peak estimation of a linear system


rng(33, 'twister')
%% generate samples
A_true = [-1 3; -1 -0.3];
f_true = @(t, x) A_true*x;

Nsample = 50;
% Nsample = 10;
box_lim = 3;

epsilon = 2;
sample = struct('t', 4, 'x', @() 3*(2*rand(2,1)-1));

[observed] = corrupt_observations(Nsample,sample, f_true, epsilon);

%% generate model
t = sdpvar(1, 1);
x = sdpvar(2, 1);

dmin = 1;
dmax = 1;
nx = 2;
f0 = zeros(nx, 1);
fw = kron(eye(nx), monolist(x, dmin, dmax)');
% sdisplay(fw)


%% create observation constraints Aw <= b
f0_func = polyval_func(f0, [t;x]);
fw_func = polyval_func(fw, [t;x]);

nw = size(fw, 2);   %number of uncertainties
m = nx*Nsample*2;   %number of constraints
A = zeros(m, nw);
b = zeros(m, 1);

counter = 0;
for i = 1:Nsample
    tcurr = observed.t(i);
    xcurr = observed.x(:, i);
    
    xdotcurr = observed.xdot_noise(:, i);
    
    f0_curr = f0_func([tcurr; xcurr]);
    fw_curr = fw_func([tcurr; xcurr]);
    
    b_pos_curr = epsilon - f0_curr - xdotcurr;
    b_neg_curr = epsilon + f0_curr + xdotcurr;
    
    A(counter+(1:2*nx), :) = [fw_curr; -fw_curr];
    b(counter+(1:2*nx)) = [b_pos_curr; b_neg_curr];
    
    counter = counter + 2*nx;
end

%% normalize uncertainty to [-1, 1]^nw

box = poly_bounding_box(A,b);
[box_out, box_center, box_half] = box_process(nw, box);

f0_center = f0 + fw*box_center;
fw_center = fw*diag(box_half);

A_scale = A*diag(box_half);
b_scale = b - A*box_center;

 



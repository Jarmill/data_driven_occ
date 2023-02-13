
%% data generation
x = sdpvar(2, 1);
u = sdpvar(1, 1);

f_true = [x(2) - x(1)^3 + x(1)^2; u]; 

%set parameters
Tmax = 5;

R0 = 0.4;           % radius
C0 = [4; 2];        % center

Ru = 0.4;
Cu = [2; 2];

epsilon = 0;

Nsample = 40;

%sampler: [x; u] (because this is controlled)
sample = struct('t', Tmax, ...
    'x', @() [C0 + R0*ball_sample(1,2)']);

DG = data_generator(sample);

observed = DG.corrupt_observations(Nsample, f_true, epsilon);


%% polytope generation
order = 3;
mlist = monolist(x, 3);

model = struct('f0', [0;0], 'fw', diag([1; 1; u])*kron(mlist', eye(3)));
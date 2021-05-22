load('three_wave_data.mat', 'W_cheb');

rng(33, 'twister')
A_scale = W.A;
b_scale = W.b;
d = size(A_scale, 2);

w = sdpvar(d, 1);


objective = w(1);
cons = A_scale*w <= b_scale;

optimize(cons, -objective)

box = poly_bounding_box(A_scale, b_scale);

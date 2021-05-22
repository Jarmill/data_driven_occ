load('three_wave_data.mat', 'W_cheb');

rng(33, 'twister')
A_scale = W_cheb.A;
b_scale = W_cheb.b;
D = A_scale ./ repmat(b_scale,[1 size(A_scale,2)]);

d = size(A_scale, 2);
% max_ind = zeros(d,1);
% min_ind = zeros(d,1);
[Mmax, Imax] = max(D, [], 1);
[Mmin, Imin] = min(D, [], 1);
% for i = 1:d
%     max_ind
% end


[I_oct, ia, ic] = unique([Imax, Imin]);

while length(I_oct) < (d+1)
    w_curr = sphere_sample(1, d);
    [m, I] = max(D*w_curr');
    
    if ~ismember(I, I_oct)
        I_oct = [I_oct, I];
    end
end


%% Attempt to cull points in D
D_oct = D(I_oct, :);

[A,b,Aeq,beq]=vert2lcon(D_oct);


D_slack = b - A*D';
D_slack_cull = all(D_slack > -1e-8, 1);
 
% nnz(D_slack_cull)
ind_cull = logical(ones(length(b_scale), 1));
ind_cull(D_slack_cull) = 0;
ind_cull(I_oct) = 1;
% ind_cull = setdiff(find(D_slack_cull), I_oct);



%this culling did not work



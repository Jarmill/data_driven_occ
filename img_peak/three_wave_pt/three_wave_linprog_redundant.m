load('three_wave_data.mat', 'W_cheb');

rng(33, 'twister')
A_scale = W_cheb.A;
b_scale = W_cheb.b;
tic
[A_red, b_red] = nontrivial_constraints(A_scale, b_scale);
[A_red2, b_red2] = nontrivial_constraints(A_red, b_red);
toc
[A_red3, b_red3] = nontrivial_constraints(A_red2, b_red2);

% [m, d] = size(A_scale);
% 
% w = sdpvar(d, 1);
% cost = sdpvar(d, 1);
% 
% objective = cost'*w;
% cons = A_scale*w <= b_scale;
% 
% Lopt = optimizer(cons, objective, sdpsettings('solver', 'mosek', 'savesolveroutput', 1), cost, w);
% 
% % [w_1, error_opt, status, dual_1] = Lopt([1; zeros(d-1, 1)]);
% 
% % candidate = logical(ones(d, 1));
% nontrivial = logical(zeros(m, 1));
% 
% tic
% tol = 1e-8;
% for i = 1:m
%     %solve linprog in current direction
%     if ~nontrivial(i)
%         [w_curr, error_optA, statusA, dual_curr] = Lopt(-A_scale(i, :)');
%         
%         active_cons = dual_curr{1} > tol;
%         n_active = nnz(active_cons);
%         nontrivial(active_cons) = true;
%     end
% end
% 
% toc
% 
% nnz(nontrivial)
% 
% A_red = A_scale(nontrivial, :);
% b_red = b_scale(nontrivial, :);
% 
% %% postprocess?
% 
% cons_red = A_red*w <= b_red;
% 
% Lopt_red = optimizer(cons_red, objective, sdpsettings('solver', 'mosek', 'savesolveroutput', 1), cost, w);
% 
% % [w_1, error_opt, status, dual_1] = Lopt([1; zeros(d-1, 1)]);
% 
% % candidate = logical(ones(d, 1));
% nontrivial_red = false(length(b_red), 1);
% 
% tic
% tol = 1e-8;
% for i = 1:length(b_red)
%     %solve linprog in current direction
%     if ~nontrivial(i)
%         [w_curr, error_optA, statusA, dual_curr] = Lopt_red(-A_red(i, :)');
%         
%         active_cons = dual_curr{1} > tol;
%         n_active = nnz(active_cons);
%         nontrivial_red(active_cons) = true;
%     end
% end
% 
% toc
% 
% nnz(nontrivial_red)
% 
% A_red2 = A_red(nontrivial_red, :);
% b_red2 = b_red(nontrivial_red, :);
% 
% 
% 

% box = poly_bounding_box(A_scale, b_scale);

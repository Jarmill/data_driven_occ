function [out_sim] = traj_eval(out, out_sim)
%TRAJ_EVAL evaluate functions over all trajectories in out_sim


nonneg_func = @(t,x)cell2mat(arrayfun(@(i)out.func.nonneg([t(i),x(i, :)]),...
    (1:length(t)),'UniformOutput',false));

cost_func = @(t, x) arrayfun(@(i)out.func.cost(x(i,:)),...    
    (1:length(t)));

v_func = @(t,x) arrayfun(@(i)out.func.v([t(i),x(i,:)]),...
    (1:length(t)));


for i = 1:length(out_sim)
    osc = out_sim{i};
    N = length(osc.t);
    
    t_curr = osc.t;
    x_curr = osc.x;
    out_sim{i}.cost = cost_func(t_curr, x_curr);
    out_sim{i}.v = v_func(t_curr, x_curr);
    out_sim{i}.nonneg= nonneg_func(t_curr, x_curr);
%     cost_curr = zeros(1, N);
%     nonneg_curr = zeros(size(out.poly.nonneg(1)), N);
% %     v_curr = zeros(1, N
%     for j = 1:N
%         
%     end
end

end


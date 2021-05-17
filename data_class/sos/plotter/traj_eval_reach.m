function [out_sim] = traj_eval(out, out_sim)
%TRAJ_EVAL evaluate functions over all trajectories in out_sim

Tmax = out.dynamics.Tmax;
% Tmax = 1;

nonneg_func = @(t,x)cell2mat(arrayfun(@(i)out.func.nonneg([t(i)/Tmax,x(i, :)]),...
    (1:length(t)),'UniformOutput',false));

if out.FREE_TERMINAL
w_func = @(t, x) arrayfun(@(i)out.func.w(t(i), x(i,:)),...    
    (1:length(t)));
else
    w_func = @(t, x) arrayfun(@(i)out.func.w(x(i,:)),...    
    (1:length(t)));
end

v_func = @(t,x) arrayfun(@(i)out.func.v([t(i)/Tmax,x(i,:)]),...
    (1:length(t)));




for i = 1:length(out_sim)
    osc = out_sim{i};
    N = length(osc.t);
    
    t_curr = osc.t;
    x_curr = osc.x;
    out_sim{i}.w = w_func(t_curr, x_curr);
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


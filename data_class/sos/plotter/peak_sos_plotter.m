classdef peak_sos_plotter < plotter_sos_interface
    %PEAK_SOS_PLOTTER Summary of this class goes here
    %   Detailed explanation goes here
    
    
    
    methods
        function obj = peak_sos_plotter(out, out_sim)
            %PEAK_SOS_PLOTTER Construct an instance of this class
            %   Detailed explanation goes here
            obj = obj@plotter_sos_interface(out, out_sim);
        end
        
        function F = nonneg_traj(obj)

            F = figure(20);
            
%             N0 = length(obj.manager.opts.get_X_init());
            
            %no other switching present, so there is only one system for
            %the occupation measure
            
            ax_label= {'$\gamma - v(t,x)$', '$v(x) - \beta^T p(x)$', '$-L_{f0} v(t,x) - b^T \zeta(t,x)$'};
            ax_title = {'Bounded v (initial)', 'Dominating cost (peak)', 'Decreasing v (occupation)'};
            
            tiledlayout(3, 1);
            for i = 1:3
                nexttile
                hold on
                
                for j = 1:length(obj.out_sim)
                    osc = obj.out_sim{j};                    
                    if i == 1
                        plot(osc.t, osc.nonneg(1:(end-2), :), 'c');
                    elseif i==2
                        plot(osc.t, osc.nonneg(end-1, :), 'c');
                    else
                        plot(osc.t, osc.nonneg(end, :), 'c');
                    end
                end
                xlabel('time', 'FontSize', obj.FS_axis)
                ylabel(ax_label{i}, 'interpreter', 'latex', 'FontSize', obj.FS_axis);
                title(ax_title{i}, 'FontSize', obj.FS_title);   
                plot(xlim, [0,0], 'k:')
            end
        end

        
        function F = cost_plot(obj)

            F = figure(23);
            hold on
            for j = 1:length(obj.out_sim)
                osc = obj.out_sim{j};                    
                plot(osc.t, osc.cost, 'c');

            end
            plot(xlim, obj.out.obj*[1,1], '--r', 'LineWidth', 3)
            xlabel('time', 'FontSize', obj.FS_axis)
            ylabel('$p(x)$', 'interpreter', 'latex', 'FontSize', obj.FS_axis);
            title('Evaluated Cost', 'FontSize', obj.FS_title);   
        end
        
    end
end


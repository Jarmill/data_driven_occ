classdef peak_sos_plotter < plotter_sos_interface
    %PEAK_SOS_PLOTTER Plotter for peak estimation
    %   Detailed explanation goes here
    
    
    
    methods
        function obj = peak_sos_plotter(out, out_sim)
            %PEAK_SOS_PLOTTER Construct an instance of this class
            %   Detailed explanation goes here
            obj = obj@plotter_sos_interface(out, out_sim);
        end
        
        function F = nonneg_traj(obj)

            F = figure(20);
            clf
%             N0 = length(obj.manager.opts.get_X_init());
            
            %no other switching present, so there is only one system for
            %the occupation measure
            
            ax_label= {'$\gamma - v(t,x)$', '$v(t,x) - \beta^T p(x)$', '$-L_{f0} v(t,x) - b^T \zeta(t,x)$'};
            ax_title = {'Bounded v (initial)', 'Dominating cost (peak)', 'Decreasing v (occupation)'};
            
            tiledlayout(3, 1);
            for i = 1:3
                nexttile
                hold on
                
                Nzeta = length(obj.out.poly.zeta);
                
                for j = 1:length(obj.out_sim)
                    osc = obj.out_sim{j};                    
                    if i == 1
                        plot(osc.t, osc.nonneg(1:(end-2-Nzeta), :), 'c');
                    elseif i==2
                        plot(osc.t, osc.nonneg(end-1-Nzeta, :), 'c');
                    else
                        plot(osc.t, osc.nonneg(end-Nzeta, :), 'c');
                    end
                end
                xlabel('time', 'FontSize', obj.FS_axis)
                ylabel(ax_label{i}, 'interpreter', 'latex', 'FontSize', obj.FS_axis);
                title(ax_title{i}, 'FontSize', obj.FS_title);   
                plot(xlim, [0,0], 'k:')
            end
        end

        function F = state_plot_2(obj, box_lim)
            if nargin == 1
                box_lim = [];
            end
            
            [F, limits] = state_plot_2@plotter_sos_interface(obj, box_lim);
            
            
            %implicit curves
%             syms y [2, 1];
%             cy = obj.out.func.cost(y) - obj.out.obj;
%             fimplicit(cy + 1e-8*sum(y), limits,  '--r', 'HandleVisibility', 'Off', 'LineWidth', 2)
        end

        function F = state_plot_3(obj, box_lim)
            if nargin == 1
                box_lim = [];
            end
            
            [F, limits] = state_plot_3@plotter_sos_interface(obj, box_lim);
            
            
            %implicit curves
            syms y [3, 1];
            cy = obj.out.func.cost(y) - obj.out.obj;
            fimplicit3(cy + 1e-8*sum(y), limits,  'r', 'HandleVisibility', 'Off', 'FaceAlpha', 0.5, 'EdgeColor', 'None')
        end

        
        function F = v_plot(obj)
            %add the gamma upper bound to the peak v plot
            F = v_plot@plotter_sos_interface(obj);
            subplot(2, 1, 1)
            plot(xlim, [1;1]*obj.out.poly.gamma, '--r', 'LineWidth', 3)
        end
        
        
        function F = obj_plot(obj)
            %plot the objective
            F = figure(23);
            clf
            hold on
            for j = 1:length(obj.out_sim)
                osc = obj.out_sim{j};                    
                plot(osc.t, osc.cost, 'c');

            end
            plot(xlim, obj.out.obj*[1,1], '--r', 'LineWidth', 3)
            xlabel('time', 'FontSize', obj.FS_axis)
            ylabel('$p(x)$', 'interpreter', 'latex', 'FontSize', obj.FS_axis);
            title('Evaluated Objective', 'FontSize', obj.FS_title);   
        end
        
    end
end


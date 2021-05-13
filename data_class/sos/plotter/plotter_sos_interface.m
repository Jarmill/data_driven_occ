classdef plotter_sos_interface < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        FS_axis = 14;       %font size of axis
        FS_title = 16;      %font size of title1
        
        manager;
    end
    
    properties(Access=protected)
        out_sim;
    end
    
    methods
        function obj = plotter_sos_interface(manager, out_sim)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
%             obj.Property1 = inputArg1 + inputArg2;
            obj.manager = manager;
            obj.out_sim = out_sim;
        end
        
        function F = state_plot(obj)
            F = figure(19);
            clf
            nx = size(obj.out_sim{1}.x, 2);
            nplt = nx;
            ax = cell(nplt, 1);
            for k = 1:nplt
                ax{k} = nexttile;
                title(['State ', num2str(k), ' vs. Time'], obj.FS_title);   
                ax_loc_curr = ['$x_', num2str(k), '$'];
                xlabel('time', 'FontSize', obj.FS_axis)
                ylabel(ax_loc_curr , 'interpreter', 'latex', 'FontSize', obj.FS_axis);
                hold on
            end
            
            for i = 1:length(obj.out_sim)
                tcurr = obj.out_sim{i}.t;
                xcurr = obj.out_sim{i}.x;
                for k = 1:nplt
                   plot(ax{k}, tcurr, xcurr(:, k), 'c') 
                end
            end
            
        end
        
%         function F = nonneg_traj(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
%             outputArg = obj.Property1 + inputArg;

%             F = figure(20);
%         end
        
        function F = nonneg_zeta(obj)
            %plot the nonnegative slack functions zeta
            F = figure(21);
            
            Nzeta = length(obj.manager.opts.W.b);
            tiledlayout(Nzeta, 1);
            for i = 1:Nzeta
                nexttile
                hold on
                for j = 1:length(obj.out_sim)
                    osc = obj.out_sim{j};                    
                    plot(osc.t, osc.nonneg(end-(Nzeta-i), :), 'c')
                end
                
                ax_loc_curr = ['$\zeta_', num2str(i), '(t,x)$'];
                xlabel('time', 'FontSize', obj.FS_axis)
                ylabel(ax_loc_curr , 'interpreter', 'latex', 'FontSize', obj.FS_axis);
                title(['Constraint ', num2str(i), ' Slack'], 'FontSize', obj.FS_title);   
                plot(xlim, [0,0], 'k:')
            end
            
            
        end
    end
end


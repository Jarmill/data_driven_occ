classdef plotter_sos_interface < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        FS_axis = 14;       %font size of axis
        FS_title = 16;      %font size of title1
        
        
    end
    
    properties(Access=protected)
        out;            %result of optimization
        out_sim;        %trajectories
    end
    
    methods
        function obj = plotter_sos_interface(out, out_sim)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
%             obj.Property1 = inputArg1 + inputArg2;
            obj.out = out;
            obj.out_sim = out_sim;
            
            
        end
        
        function F = v_plot(obj)

            F = figure(18);
            clf
            subplot(2, 1, 1)
            hold on

            for j = 1:length(obj.out_sim)
                osc = obj.out_sim{j};                    
                plot(osc.t, osc.v, 'c');

            end
            xlabel('time', 'FontSize', obj.FS_axis)
            ylabel('$v(t,x)$', 'interpreter', 'latex', 'FontSize', obj.FS_axis);
            title('Auxiliary Function', 'FontSize', obj.FS_title);   
            
            subplot(2,1,2)
            hold on
            for j = 1:length(obj.out_sim)
                osc = obj.out_sim{j};                    
                plot(osc.t(1:end-1), diff(osc.v), 'c');

            end
            xlabel('time', 'FontSize', obj.FS_axis)
            ylabel('$\Delta v(t,x)$', 'interpreter', 'latex', 'FontSize', obj.FS_axis);
            title('Auxiliary Function Decrease', 'FontSize', obj.FS_title);   
            
            
        end
        
        function F = state_plot(obj)
            F = figure(19);
            clf
            nx = size(obj.out_sim{1}.x, 2);
            nplt = nx;
            ax = cell(nplt, 1);
            for k = 1:nplt
                ax{k} = nexttile;
                title(['State ', num2str(k), ' vs. Time'], 'FontSize', obj.FS_title);   
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
        
        function [F, limits] = state_plot_2(obj, box_lim)
            F = figure(30);
            clf
            hold on 
            for i = 1:length(obj.out_sim)
                tcurr = obj.out_sim{i}.t;
                xcurr = obj.out_sim{i}.x;
                plot(xcurr(:, 1), xcurr(:, 2), 'c')
            end
            for i = 1:length(obj.out_sim)
                tcurr = obj.out_sim{i}.t;
                xcurr = obj.out_sim{i}.x;
                scatter(xcurr(1, 1), xcurr(1, 2), 100, 'k')
            end                        
            
            xlabel('$x_1$', 'interpreter', 'latex', 'FontSize', obj.FS_axis);
            ylabel('$x_2$', 'interpreter', 'latex', 'FontSize', obj.FS_axis);
            title('Phase Plane', 'FontSize', obj.FS_title);   
            
            
            if (nargin == 2 ) && ~isempty(box_lim)
                box = box_process(2, box_lim);
                xlim(box(1, :));
                ylim(box(2, :));
                limits = [box(1, :), box(2, :)];
            else
                limits = [xlim, ylim];
            end
            
            pbaspect([diff(xlim), diff(ylim), 1])
            
        end
 
        
        function F = nonneg_zeta(obj)
            %plot the nonnegative slack functions zeta
            F = figure(21);
            
            Nzeta = length(obj.out.poly.zeta);
            hold on
            for j = 1:length(obj.out_sim)
                osc = obj.out_sim{j};                    
                plot(osc.t, osc.nonneg(end-Nzeta:end, :), 'c')
            end

            ax_loc_curr = ['$\zeta(t,x)$'];
            xlabel('time', 'FontSize', obj.FS_axis)
            ylabel(ax_loc_curr , 'interpreter', 'latex', 'FontSize', obj.FS_axis);
            title(['Constraint Slack'], 'FontSize', obj.FS_title);   
            plot(xlim, [0,0], 'k:')
            
            
        end
    end
end


classdef data_generator
    %DATA_GENERATOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        sample;        
        observed;
    end
    
    methods
        function obj = data_generator(sample, vars)
            %DATA_GENERATOR Construct an instance of this class
            %   Detailed explanation goes here
            
            %pass in the sampler in (t, x) to generate datapoints            
            obj.sample = sample;            
        end
        
        %% generate constraints
        function [observed] = corrupt_observations(obj, Np, f, epsilon)
        %Generate Nsample corrupted observations of the system x'(t) = f(t,x)
        %with an L_infinity bounded noise of epsilon
        %Inputs:
        %   Np:         number of sampled points to run
        %   sampler:    struct with fields 't' and 'x' of function_handles
        %   f:          system dynamics
        %   epsilon:    noise bound
        %
        %Output:        struct 'observed'
        %                   t:          time points
        %                   x:          state points
        %                   xdot_true:  noiseless observations
        %                   xdot_noise: noisy observations
        %                   epsilon:    noise bound


        %% initial fill of output


        x0 = pull_x();
        nx = length(x0);
        observed = struct('t', zeros(1, Np), 'x', zeros(nx, Np), ...
            'xdot_true', zeros(nx, Np), 'xdot_noise', zeros(nx, Np), 'epsilon', epsilon);


        noise_sample = (2*rand(nx, Np)-1)*epsilon;


        %% generate sampled points

        if isnumeric(obj.sample.t)
            Tmax = obj.sample.t;
            obj.sample.t = @() rand()*Tmax;
        end


        for i = 1:Np
           tcurr = obj.sample.t();
           xcurr = pull_x();


           xdotcurr = f(tcurr, xcurr);


           observed.t(i) = tcurr;
           observed.x(:, i) = xcurr;
           observed.xdot_true(:, i) = xdotcurr;
        %    observed.xdot_noise(:, i) = xdotcurr + noise_obj.sample(:, i);
        end

        observed.xdot_noise = observed.xdot_true + noise_sample;

            function x0 = pull_x()
                if isa(obj.sample.x, 'function_handle')
                    x0 = obj.sample.x();      
                else
                    %numeric
                    %assume that the dimensions of x and w are compatible if arrays are
                    %passed in.
                    if size(obj.sample.x, 2) == 1
                        %single point
                        x0 = obj.sample.x;
                    else
                        %array, so Ns = size(obj.sampler, 2)
                        x0 = obj.sample.x(:, i);
                    end
                end
            end


        end


        function [model] = poly_model(obj, vars, dmin, dmax)
            %POLY_MODEL generate polynomial model
            %uncertainty based on standard monomials
            %dynamics are x' = f0(t,x) + sum_i w_i fw_i(t,x)
            %assum
            
            if nargin <= 3
                dmax = dmin;
                dmin = 0;
            end
            
            if isstruct(vars)
                nx = length(vars.x);
                vars = [vars.t; vars.x];
            else
                nx = length(vars);
            end
            
%             nx =length(vars);
            f0 = zeros(nx, 1);
            fw = kron(eye(nx), monolist(vars, dmin, dmax)');         
            model = struct('f0', f0, 'fw', fw);
        end
        
        function [W] = data_cons(obj, model, vars, observed)
            %DATA_CONS linear constraint to form agreement with corrupted
            %observed data
            
            %form functions to evaluate dynamics at data 
            if isstruct(vars)
                TIME_INDEP = 0;
                nx = length(vars.x);
                vars = [vars.t; vars.x];
            else
                TIME_INDEP = 1;
                nx = length(vars);
            end
            f0_func = polyval_func(model.f0, vars);
            fw_func = polyval_func(model.fw, vars);
            
            Nsample = length(observed.t);
            nw = size(model.fw, 2);   %number of uncertainties
            m = nx*Nsample*2;   %number of constraints
            A = zeros(m, nw);
            b = zeros(m, 1);

            counter = 0;
            
            %form the linear constraints
            for i = 1:Nsample
                tcurr = observed.t(i);
                xcurr = observed.x(:, i);

                
                                
                xdotcurr = observed.xdot_noise(:, i);

                if TIME_INDEP
                    f0_curr = f0_func([xcurr]);
                    fw_curr = fw_func([xcurr]);
                else
                    f0_curr = f0_func([tcurr; xcurr]);
                    fw_curr = fw_func([tcurr; xcurr]);
                end
                %ensure the correct signs over here
                b_pos_curr = observed.epsilon - f0_curr + xdotcurr;
                b_neg_curr = observed.epsilon + f0_curr - xdotcurr;

                A(counter+(1:2*nx), :) = [fw_curr; -fw_curr];
                b(counter+(1:2*nx)) = [b_pos_curr; b_neg_curr];

                counter = counter + 2*nx;
            end         
            
            W = struct('A', A, 'b', b);
            
        end
        
        %% Process constraints
        function [model_out, W_cheb] = center_cheb(obj, model, W)
            %% Center the polytope
            %box
            % box = poly_bounding_box(A,b);
            % [box_out, box_center, box_half] = box_process(nw, box);
            % f0_center = f0 + fw*box_center;
            % fw_center = fw*diag(box_half);
            % 
            % A_scale = A*diag(box_half);
            % b_scale = b - A*box_center;

            %chebyshev center of polytope
            [c,r] = chebycenter(W.A,W.b);
            A_scale = W.A;
            b_scale = W.b - W.A*c;

            f0_center = model.f0 + model.fw*c;
            fw_center = model.fw;
            
            model_out = struct('f0', f0_center, 'fw', fw_center);
            W_cheb = struct('A', A_scale, 'b', b_scale);

        end
        
        
        function [W_reduced] = reduce_constraints(obj, W_in)
            %identify redundant constraints
            %code from noredund.m by Michael Kleder (2006)
            
            %assume that the constraints have already been
            %chebyshev-centered
            A_scale = W_in.A;
            b_scale = W_in.b;
            D = A_scale ./ repmat(b_scale,[1 size(A_scale,2)]);
            %number of points in dual polytope's convex hull capped at the
            %number of constraints
            [k, vol] = convhulln(D);
            % record which constraints generate points on the convex hull
            nr = unique(k(:));
%             A_scale_orig = A_scale;
%             b_scale_orig = b_scale;
            
            %export reduced constraint set
            W_reduced = struct;
            W_reduced.A=A_scale(nr,:);
            W_reduced.b=b_scale(nr);

        end
        
        function [model_cheb, W_red] = reduced_model(obj, observed, vars, dmin, dmax)
            
            %shortcut many of the functions
            %given polynomial degree bounds and observed data, form reduced
            %constraints
            
            if nargin <= 4
                dmax = dmin;
                dmin = 0;
            end
            
            model = obj.poly_model(vars, dmin, dmax);
            W = obj.data_cons(model, vars, observed);
            [model_cheb,W_cheb] = obj.center_cheb(model, W);
            W_red = obj.reduce_constraints(W_cheb);
        end
        
        function [w_handle, box] = make_sampler(obj, W)
            %use rejection sampling to sample points from a polytope
            box = poly_bounding_box(W.A,W.b);
            w_handle = @() rej_sample_poly(W.A, W.b, box);
        end
        
        %% plot vector fields
        
        function F = data_plot_2(obj, observed)
            F=figure(2);
            clf
            hold on
            quiver(observed.x(1, :), observed.x(2, :), observed.xdot_true(1, :), observed.xdot_true(2, :))
            quiver(observed.x(1, :), observed.x(2, :), observed.xdot_noise(1, :), observed.xdot_noise(2, :))
            axis square
            legend({'Ground Truth', 'Noisy Data'}, 'FontSize', 12, 'location', 'northwest')
            xlabel('$x_1$', 'interpreter', 'latex', 'FontSize', 12);
            ylabel('$x_2$', 'interpreter', 'latex', 'FontSize', 12);          
            title(['Noisy Observations with \epsilon=', num2str(observed.epsilon)], 'FontSize', 16)
        end
        
        
        function F = data_plot_3(obj, observed)
            F=figure(2);
            clf
            hold on
            quiver3(observed.x(1, :), observed.x(2, :),observed.x(3, :), observed.xdot_true(1, :), observed.xdot_true(2, :),observed.xdot_true(3, :))
            quiver3(observed.x(1, :), observed.x(2, :),observed.x(3, :), observed.xdot_noise(1, :), observed.xdot_noise(2, :),observed.xdot_noise(3, :))
            axis square
            legend({'Ground Truth', 'Noisy Data'}, 'FontSize', 12, 'location', 'northwest')
            xlabel('$x_1$', 'interpreter', 'latex', 'FontSize', 12);
            ylabel('$x_2$', 'interpreter', 'latex', 'FontSize', 12);          
            ylabel('$x_3$', 'interpreter', 'latex', 'FontSize', 12);          
            view(3)
            title(['Noisy Observations with \epsilon=', num2str(observed.epsilon)], 'FontSize', 16)
        end
    end
end


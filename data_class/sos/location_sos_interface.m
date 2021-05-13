classdef location_sos_interface < handle
    %LOCATION_SOS_INTERFACE Abstract wrapper for location methods for the
    %sum-of-squares optimization problems (in yalmip)    
    
    properties
        %% properties of run
        
%         supp; %support set

        opts = []; %includes support set
        
        vars = struct('t', [], 'x', []);
        
    end
    
    methods
        function obj = location_sos_interface(opts)
            %LOCATION_SOS_INTERFACE Construct an instance of this class
            %   Detailed explanation goes here
            obj.opts = opts;
            obj.vars = struct('t', opts.t, 'x', opts.x);
            
            if isempty(obj.opts.t)
                obj.opts.t = sdpvar(1, 1);
            end
            
            if ~isempty(obj.opts.W)
                obj.opts.W.b = reshape(obj.opts.W.b, [], 1);
            end
            
            
            
        end
        
        %% Set up the program
        function [prog]= make_program(obj, d)
            %Form the SOS program associated with this problem
            [poly_var, coeff_var] = obj.make_poly(d);
            
%             nonneg = obj.form_nonneg(poly_var);                       
            [coeff, con, nonneg] = obj.make_cons(d, poly_var);    
            
            poly_var.nonneg = nonneg;
            
            [objective] = obj.get_objective(poly_var);            
            
            %add coefficients of variable polynomials terms to coefficients
            coeff = [coeff_var; coeff];
            
            prog = struct('poly', poly_var, ...
                'objective', objective, 'coeff', coeff, 'con', con,... 
                'order', d/2);
        end
        
        function [poly_out, coeff_out] = make_poly(obj,d)
            %MAKE_POLY Create polynomials v and zeta for SOS programs
            %subclasses will define further polynomials
            
            t = obj.opts.t;
            x = obj.opts.x;
            [v, cv] = polynomial([t; x], d);

            %every application includes an initial set valuation of v
            v0 = replace(v, t, 0);
            
            %The higher programs may create other variables
            %   gamma:          peak estimation
            %   w(x):           reachable set
            %   gamma and w(x): distance estimation
            %             gamma = sdpvar(1,1);

            %the number of half-space constraints in the uncertainty
            
            
            zeta = [];
            coeff_zeta = [];
            
            if ~isempty(obj.opts.W)

                m = length(obj.opts.W.b);
                d_altern = d;
                for i = 1:size(obj.opts.fw, 1)
                    d_altern = max((2*ceil(d/2 + degree(obj.opts.fw(i, :))/2-1)), d_altern); %figure out degree bounds later
                end

                for i = 1:m
                    [pzeta, czeta] = polynomial([t; x], d_altern);
                    zeta = [zeta; pzeta];
                    coeff_zeta = [coeff_zeta; czeta];
                end
            end
            
            poly_out=struct('v', v, 'zeta', zeta, 't', t, 'x', x, 'v0', v0);
            coeff_out = [cv; coeff_zeta];
        end
        
        function [coeff_lie, cons_lie, nonneg_lie] = make_lie_con(obj, d, poly)        
            %make constraints for the SOS program
            %in this interface, only perform the Lie derivative <= 0
            %decomposition                        
            
            %polynomials
            v = poly.v;
            zeta = poly.zeta;
            
            %variables
            t = obj.vars.t;
            x = obj.vars.x;
            
            if obj.opts.scale
                scale_weight = obj.opts.Tmax;
                obj.opts.f0 = replace(obj.opts.f0, t, scale_weight*t);
                obj.opts.fw = replace(obj.opts.fw, t, scale_weight*t);
            else
                scale_weight = 1;
            end
            

            %region of validity [0, T] times X for dynamics
            Xall = obj.opts.get_all_supp();
            
            
            %the original constraint is that Lie v <= 0
            %This implies that the set Lie v > 0 is empty for (t,x,w)
            %The resultant theorem of alternatives yield a condition with
            %Lie v - (terms) >= 0
            %be careful of the signs

            %lie derivative with no uncertainty
            dvdx = jacobian(v, x);                        
            
            if isempty(obj.opts.fw)
                %no uncertainty at all in dynamics
                Lv0 = dvdx*(scale_weight*obj.opts.f0) + jacobian(v, t);
                [cons_lie, coeff_lie] =  obj.make_psatz(d, Xall, -Lv0, [t;x]);
%                 constraint_psatz(-Lv0, Xall, [t;x], d);
                cons_lie = cons_lie:'Lie standard (no input)';
                
                nonneg_lie = -Lv0;
            else
                %there is uncertainty
               
                
                %uncertainty constraint Aw <= b

                A = obj.opts.W.A;
                b = obj.opts.W.b;
                [m, num_input] = size(A); %number of constraints
            
                
                if isempty(obj.opts.f0)
                    coeff_lie = [];
                    cons_lie = [];
                    nonneg_lie = [];
                else
                    Lv0 = dvdx*(scale_weight*obj.opts.f0) + jacobian(v, t);
%                       %Aw >= b
%                     [consf0, coefff0] =  obj.make_psatz(d, Xall, Lv0-b'*zeta, [t;x]);

                    %Aw <= b
%                     [consf0, coefff0] =  obj.make_psatz(d, Xall, Lv0+b'*zeta, [t;x]);
                    [consf0, coefff0] =  obj.make_psatz(d, Xall, -Lv0-b'*zeta, [t;x]);


                    %output
                    coeff_lie = coefff0;
                    cons_lie = consf0:'Lie base duality';
                    nonneg_lie = -Lv0 - b'*zeta;
                end
                %terms with uncertainty
                %these are equality constraints in polynomial coefficients
                %TODO: implement sparsity
                for j = 1:num_input
                    %current dynamics
                    Lvj = dvdx*(scale_weight*obj.opts.fw(:, j));

                    %current linear constraint
                    Aj = A(:, j);

                    %equality constraint
%                     %Aw >= b
%                     equal_j = Lvj + Aj'*zeta;                
                    
                    %Aw <= b
                    equal_j = Lvj - Aj'*zeta;      
                    cons_equal = (coefficients(equal_j, [t; x]) == 0);
                    cons_lie = [cons_lie; cons_equal:['Lie input ', num2str(j), ' duality']];
                end



                %sos constraints on zeta
                %zeta are >=0 over the dynamics support region Xall 
                for j = 1:m
                    zeta_curr = zeta(j);
                    d_zeta = degree(zeta_curr);
                    
                    [conszeta, coeffzeta] =  obj.make_psatz(d_zeta, Xall, zeta_curr, [t;x]);
%                     [pzeta, conszeta, coeffzeta] = constraint_psatz(zeta_curr, Xall, [t;x], d_zeta);

                    coeff_lie = [coeff_lie; coeffzeta];
                    cons_lie = [cons_lie; conszeta:['zeta ', num2str(j), ' sos']];
                    nonneg_lie = [nonneg_lie; zeta_curr];
                end
            end
            
            %output [coeff_lie, cons_lie]
        end
        
        
        
        
        
        %% Solve the program        
        function [out] = solve_program(obj, prog)
            %solve SOS program in YALMIP, return solution    
            
            sdp_opts = sdpsettings('solver', obj.opts.solver, 'verbose', obj.opts.verbose);
            sdp_opts.sos.model = 2;
            
            [sol, monom, Gram, residual] = solvesos(prog.con, prog.objective, sdp_opts, prog.coeff);
            
            out = struct('poly', [], 'problem', sol.problem, 'sol', [], 'block', [], 'func', []);
            if sol.problem == 0
                %the sets X0 and X1 are disconnected in time range [0, T]
                [out.poly, out.func] = obj.recover_poly(prog.poly);
                out.sol = sol;   
                out.block = struct;
                out.block.monom = monom;
                out.block.Gram = Gram;
                out.block.residual = residual;     
                out.obj = value(prog.objective);
                out.order = prog.order;                
                
            end
            
            out.dynamics = obj.package_dynamics(out.func);
        end
        
        
        function out = run(obj, order)
            %the main routine, run the SOS program at the target order
            d = 2*order;
            prog = obj.make_program(d);
            out= solve_program(obj, prog);
        end
        
        %% Recover from the solution
        function [poly_eval, func_eval] = recover_poly(obj, poly_var, nonneg)
        %recover polynomials from SOS certificate
        
            t = poly_var.t;
            x = poly_var.x;
            n = length(poly_var.x);
            
            %solved coefficients of v and zeta
            [cv,mv] = coefficients(poly_var.v,[poly_var.t; poly_var.x]);
            v_eval = value(cv)'*mv;                                   
            v_eval = replace(v_eval, poly_var.t, poly_var.t/obj.opts.Tmax);
            %remember to scale by time 
            
            [cz, mz] = coefficients(poly_var.zeta,[poly_var.t; poly_var.x]);
            if n == 1
                zeta_eval = value(cz)'*mz;
            else
                zeta_eval = value(cz)*mz;
            end
            zeta_eval = replace(zeta_eval, poly_var.t, poly_var.t/obj.opts.Tmax);
            
            %evaluations of v at initial and terminal times
            v0 = replace(v_eval, t, 0);            
            
            %terminal set
%             if obj.opts.scale
%                 v1 = replace(v_eval, t, 1);
%             else
%                 v1 = replace(v_eval, t, obj.opts.Tmax);
%             end

            %nonnegative evaluation
            [cnn,mnn] = coefficients(poly_var.nonneg,[poly_var.t; poly_var.x]);
            nn_eval = value(cnn)*mnn;
            nn_eval = replace(nn_eval, poly_var.t, poly_var.t/obj.opts.Tmax);

                        
            poly_eval = struct('v', v_eval, 'zeta', zeta_eval, 'v0', v0, 'nonneg', nn_eval);
            
            % form functions using helper function 'polyval_func'
            func_eval = struct;
            func_eval.v = polyval_func(v_eval, [t; x]);
            func_eval.zeta = polyval_func(zeta_eval, [t; x]);
            
%             func_eval.v0 = polyval_func(v0, [x]);            
            func_eval.nonneg = polyval_func(nn_eval, [t;x]);            
        end
    
        function dynamics = package_dynamics(obj, func_in)
            %package up dynamics for use in the (old) sampler
            %peak/sampler
            
            
            dynamics = struct('Tmax', obj.opts.Tmax, 'discrete', 0);
            
            %TODO: replace with time independence
            dynamics.time_indep = 0;
            %iterate through points array and evaluate nonnegative
            %functions
            dynamics.nonneg_val = func_in.nonneg;
            %truly terrible code
            %evaluate nonnegative functions at each point [t, x]
            %functional design pattern
            dynamics.nonneg = @(t,x,w,d) cell2mat(arrayfun(@(i) dynamics.nonneg_val([t(i); x(:, i)]),(1:length(t)), 'UniformOutput', false))';
            
            
            
            %This is old sampler code as a kludge. 
            
            %event handle
            dynamics.Xval = constraint_func(obj.opts.X, obj.vars.x);
            dynamics.event = {@(t,x,w) support_event(t, x, dynamics.Xval, ...
                0, obj.opts.Tmax)};
            
            
            %dynamics handle (time-varying polytopic uncertainty called 'd'
            %here)
            dynamics.f0= polyval_func(obj.opts.f0, [obj.vars.t; obj.vars.x]);
            dynamics.fw = polyval_func(obj.opts.fw, [obj.vars.t; obj.vars.x]);
            dynamics.f = {@(t,x,w,d,b) dynamics.f0([t; x]) + dynamics.fw([t; x])*d};
            
        end
        
        %% Helper functions
        function X_cell = prep_space_cell(obj, X)
            %PREP_SPACE_CELL: Split up X into cells
            if iscell(X)
                %X is already a cell, nothing to be done
                X_cell = cell(length(X), 1);
                for i = 1:length(X)
                    X_cell{i} = fill_constraint(X{i});
                end
            else
                if isstruct(X)
                    %wrap up the structure in a cell
                    X_cell = {X};
                else
                    %X is supported at discrete points
                    %Each column of X is a possible origin datapoint
                    if size(X, 2) == 1
                        X_cell = {X};
                    else
                        Xt = X';
                        Xt_cell = mat2cell(Xt, ones(size(X, 2), 1));
                        X_cell = cellfun(@(x) x', Xt_cell, 'UniformOutput', false);
                    end
                end
            end
        end
        
        
        function [cons, coeff] = make_psatz(obj, d, X, f, vars)
            %MAKE_PSATZ a positivestellensatz expression to impose that
            %f(vars) >= 0 on the region X
            if isstruct(X)
                %X1 is a set, nonnegativity of region
%                 set1 = true;
                X = fill_constraint(X);

                [p, cons, coeff] = constraint_psatz(f, X, vars, d);        
            else  
                %X is a point, nonnegativity of evaluation
                f_pt = replace(f, vars, X);
                cons = f_pt>= 0; 
                coeff = [];
            end
            
        end                                
        
    end
    
    
    methods(Abstract)
        
        get_objective(obj, poly)
        %fetch the SOS objective to minimize
        
        make_cons(obj, poly)
        %make constraints for the SOS program
        
        make_init_con(obj, poly)
        %Constraint on initial measure
        
        make_term_con(obj, poly)
        %constraint on terminal measure
    end
end


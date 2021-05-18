classdef reach_sos < location_sos_interface
    %REACH_SOS location handle for reachable set estimation
    %includes fixed-time and free-time reachable set
    
    
    properties
        mom_handle;     %a function mom_handle(d) computing moments of the 
                        %lebesgue distribution on X (or [0,T]times X for 
                        %free terminal time) from order 0...d
        FREE_TIME;      %free terminal time (1) or fixed terminal time Tmax (0)
        
    end
    
    methods
        function obj = reach_sos(opts, mom_handle, FREE_TIME)
            %REACH_SOS Construct an instance of this class
            %opts:      loc_sos_options data structure
            %mom_handle:function handle for Lebesgue moments
            %FREE_TIME: fixed or free terminal time
            obj@location_sos_interface(opts);
            obj.mom_handle = mom_handle;
            if nargin == 2
                obj.FREE_TIME = 0;
            end
        end
        
        function [coeff, con, nonneg] = make_cons(obj, d, poly)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            
            %lie derivative constraint
            [coeff_lie, con_lie, nonneg_lie] = obj.make_lie_con(d, poly);
            
            
            %objective constraint (terminal)
            [coeff_term, con_term, nonneg_term] = obj.make_term_con(d, poly);
            
            %initial constraint in v0
            [coeff_init, con_init, nonneg_init] = obj.make_init_con(d, poly);
            
            
            %package up the constraints
            coeff = [coeff_init; coeff_lie; coeff_term];
            con = [con_init; con_lie; con_term];
            nonneg = [nonneg_init; nonneg_term; nonneg_lie];
        end
        
        function [poly_out, coeff_out] = make_poly(obj,d)
            %MAKE_POLY Create polynomials v and zeta for SOS programs
            %includes 'w' for reachable set estimation
            %w: outer approximation to indicator function on reach. set
            %   sublevel set {x | w(x) >= 1}            
            [poly_out, coeff_out] = make_poly@location_sos_interface(obj,d);
            
            
            %now define new terms
            %polynomial w
            t = obj.opts.t;
            x = obj.opts.x;
            
            if obj.FREE_TIME
                [w, cw] = polynomial([t;x],d);
            else
                [w, cw] = polynomial(x, d);
            end                                    
            
            %add terms to polynomial storage structure
            poly_out.w = w;            
            poly_out.cw = cw;
            coeff_out = [coeff_out; cw];
            
        end
        
        function [coeff, con, nonneg] = make_term_con(obj, d, poly)
            %Constraint on terminal measure
            %w(x) is an indicator function
            %w(x) + v(T,x) >= 1 on X
            
            
            if obj.FREE_TIME
                Xall = obj.opts.get_all_supp();
                vT = replace(poly.v, obj.vars.t, 1);
                term_pos = vT + poly.w - 1;
                [con_term, coeff_term] = obj.make_psatz(d, Xall, term_pos, [poly.t; poly.x]);
                [con_w, coeff_w] = obj.make_psatz(d, Xall, poly.w, [poly.t; poly.x]);
            else
                Xall = obj.opts.get_X();
                term_pos = poly.v + poly.w - 1;
                [con_term, coeff_term] = obj.make_psatz(d, Xall, term_pos, [poly.t; poly.x]);
                [con_w, coeff_w] = obj.make_psatz(d, Xall, poly.w, [poly.x]);
            end                        
            
            
            
            coeff = [coeff_term; coeff_w];
            con = [con_term:'terminal w', con_w:'nonneg w'];            
            
            nonneg = [term_pos; poly.w];
        end
    
        
        function [coeff, con, nonneg] = make_init_con(obj, d, poly)
            %Constraint on initial measure
            %auxiliary function v(0,x) is nonpositive on X0
            
            X0 = prep_space_cell(obj.opts.X_init);
            
            coeff = [];
            con = [];
            nonneg = [];
            init_pos = -poly.v0;
            for i = 1:length(X0)
                [con_curr, coeff_curr] = obj.make_psatz(d, X0{i}, init_pos, poly.x);
                
                if length(X0) == 1
                    init_tag = 'initial';
                else
                    init_tag = ['initial ', num2str(i)];
                end
                
                coeff = [coeff; coeff_curr];
                con   = [con; con_curr:init_tag];
%                 nonneg = [nonneg; init_pos];

                %the nonnegative function should hold for all time
                nonneg = [nonneg; -poly.v];
            end
            
        end
        
        function objective = get_objective(obj, poly_var, d)
            %fetch the SOS objective to minimize
            %volume of reachable set
  
        %integral of w(x) w.r.t. lebesgue measure
            mom_leb = obj.mom_handle(d);
            
%             [cw, mw] = coefficients(poly_var.w, poly_var.x);
            
            
            objective = poly_var.cw'*mom_leb;                        
        end    
        
        function [poly_val, func_eval] = recover_poly(obj, poly_var)
            %recover polynomials from SOS certificate            
            [poly_val, func_eval] = recover_poly@location_sos_interface(obj, poly_var);
               
            func_eval.FREE_TIME = obj.FREE_TIME;
            if obj.FREE_TIME
                [cw,mw] = coefficients(poly_var.w,[poly_var.t; poly_var.x]);                                
            else
                [cw,mw] = coefficients(poly_var.w, poly_var.x);
            end
            w_eval = value(cw)'*mw;    
            
            poly_val.w = w_eval;
            func_eval.w = polyval_func(w_eval, [poly_var.t; poly_var.x]);               
        end
        
    end
    
end


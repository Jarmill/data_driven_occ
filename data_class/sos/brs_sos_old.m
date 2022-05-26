classdef brs_sos_old < location_sos_interface
    %BRS_SOS location handle for backwards reachable (BRS) set estimation 
    %given an external control process. Equivalent to control design for 
    %region of attraction maximization (ROA)
    
    %free terminal time is NOT included here.
    
    %this probably could have been designed through inheritence with 
    %brs_sos < reach_sos. oops.
    
    
    properties
        mom_handle;     %a function mom_handle(d) computing moments of the 
                        %lebesgue distribution on X from order 0...d        
    end
    
    methods
        function obj = brs_sos_old(opts, mom_handle)
            %BRS_SOS Construct an instance of this class
            %opts:      loc_sos_options data structure
            %mom_handle:function handle for Lebesgue moments            
            obj@location_sos_interface(opts);
            obj.mom_handle = mom_handle;
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
            
            [w, cw] = polynomial(x, d);
                                 
            
            %add terms to polynomial storage structure
            poly_out.w = w;            
            poly_out.cw = cw;
            coeff_out = [coeff_out; cw];
            
        end
        
        function [coeff, con, nonneg] = make_init_con(obj, d, poly)
            %Constraint on initial measure
            %w(x) is an indicator function
            %w(x) >= 1 + v(0,x)on X
            
            

            Xall = obj.opts.get_X();
            v0 = replace(poly.v, obj.vars.t, 0);
            init_pos = -v0 + poly.w - 1;
            [con_init, coeff_init] = obj.make_psatz(d, Xall, init_pos, [poly.x]);
            [con_w, coeff_w] = obj.make_psatz(d, Xall, poly.w, [poly.x]);                   


            
            coeff = [coeff_init; coeff_w];
            con = [con_init:'initial w', con_w:'nonneg w'];            
            
            nonneg = [init_pos; poly.w];
        end
    
        
        function [coeff, con, nonneg] = make_term_con(obj, d, poly)
            %Constraint on terminal measure
            %auxiliary function v(T,x) is nonpositive on X0
            
            XT = prep_space_cell(obj.opts.X_term);
            
            coeff = [];
            con = [];
            nonneg = [];
%             init_pos = -poly.vT;
            vT = replace(poly.v, obj.vars.t, 1);
            for i = 1:length(XT)
                [con_curr, coeff_curr] = obj.make_psatz(d, XT{i}, vT, poly.x);
                
                if length(XT) == 1
                    term_tag = 'terminal';
                else
                    term_tag = ['terminal ', num2str(i)];
                end
                
                coeff = [coeff; coeff_curr];
                con   = [con; con_curr:term_tag];
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
            
            
            objective = poly_var.cw'*mom_leb;                        
        end    
        
        function [poly_val, func_eval] = recover_poly(obj, poly_var)
            %recover polynomials from SOS certificate            
            [poly_val, func_eval] = recover_poly@location_sos_interface(obj, poly_var);
               
           [cw,mw] = coefficients(poly_var.w, poly_var.x);
            w_eval = value(cw)'*mw;    
            
            poly_val.w = w_eval;
            func_eval.w = polyval_func(w_eval, [poly_var.t; poly_var.x]);               
        end
        
    end
    
end


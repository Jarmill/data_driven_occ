classdef peak_sos < location_sos_interface
    %LOCATION_PEAK_SOS location handle for peak estimation
    % includes maximin objectives for safety margins
    
    properties
        p; %objective for peak estimation (possibly maximin if p is an array)
    end
    
    methods
        function obj = peak_sos(opts, p)
            %LOCATION_PEAK_SOS Construct an instance of this class
            %opts:  loc_sos_options data structure
            %p:     an array of sdpvars defining the objective to maximize
            obj@location_sos_interface(opts);
            obj.p = p;
        end
        
        function [poly_out, coeff_out] = make_poly(obj,d)
            %MAKE_POLY Create polynomials v and zeta for SOS programs
            %includes 'gamma' and 'beta' for peak estimation
            %gamma: objective bound
            %beta:  maximin peak constraint terms
            
            [poly_out, coeff_out] = make_poly@location_sos_interface(obj,d);
            
            
            %now define new terms
            gamma = sdpvar(1,1);
            
            if length(obj.p)>1
                beta = sdpvar(length(obj.p), 1);
            else
                beta = [];
            end
            
            %add terms to polynomial storage structure
            poly_out.gamma = gamma;
            poly_out.beta = beta;
            
            coeff_out = [coeff_out; gamma; beta];
            
        end
        
        function [coeff, con] = make_cons(obj, d, poly)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            
            %lie derivative constraint
            [coeff_lie, con_lie] = obj.make_lie_con(d, poly);
            
            
            %objective constraint (terminal)
            [coeff_term, con_term] = obj.make_term_con(d, poly);
            
            %initial constraint in v0
            [coeff_init, con_init] = obj.make_init_con(d, poly);
            
            
            %package up the constraints
            coeff = [coeff_lie; coeff_term; coeff_init];
            con = [con_lie; con_term; con_init];
        end
        
        function [coeff, con] = make_init_con(obj, d, poly)
            %Constraint on initial measure
            %auxiliary function v(t,x) upper bounded by gamma
            
            X0 = prep_space_cell(obj.opts.X_init);
            
            coeff = [];
            con = [];
            
            init_pos = poly.gamma - poly.v0;
            for i = 1:length(X0)
                [con_curr, coeff_curr] = obj.make_psatz(d, X0{i}, init_pos, poly.x);
                
                if length(X0) == 1
                    init_tag = 'initial';
                else
                    init_tag = ['initial ', num2str(i)];
                end
                
                coeff = [coeff; coeff_curr];
                con   = [con; con_curr:init_tag];
            end
        end
        
        function [coeff, con] = make_term_con(obj, d, poly)
            %Constraint on terminal measure
            %auxiliary function v(t, x) upper bounds objective beta'p(x)
            
%             XT = prep_space_cell(obj.opts.X_term);
            
            %peak estimation is over all time
            Xall = obj.opts.get_all_supp();
            
            coeff = [];
            con = [];
            
            
            if isempty(poly.beta)
                %only one objective
                term_pos = poly.v - obj.p;
                con_beta = [];
            else
                term_pos = poly.v - poly.beta'*obj.p;
                con_beta = [sum(poly.beta)==1; poly.beta >= 0]:'maximin';
            end
            
            [con_term, coeff_term] = obj.make_psatz(d, Xall, term_pos, [poly.t; poly.x]);
            
            coeff = [coeff_term];
            con = [con_term:'peak'; con_beta];
            
        end
        
        
        function objective = get_objective(obj, poly_var)
            %fetch the SOS objective to minimize
            objective = poly_var.gamma;
        end
        
        function [poly_val, func_eval] = recover_poly(obj, poly_var)
            %recover polynomials from SOS certificate
            
            [poly_val, func_eval] = recover_poly@location_sos_interface(obj, poly_var);
            
            poly_val.gamma = value(poly_var.gamma);
            poly_val.beta = value(poly_var.beta);
        end
        
    end
end

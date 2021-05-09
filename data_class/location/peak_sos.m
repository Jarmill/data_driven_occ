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
        
        function [coeff, con] = obj.make_cons(d, poly)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            
            %lie derivative constraint
            [coeff_lie, con_lie] = obj.make_lie_con(d, poly);
            
            
            %objective constraint (terminal)
            
            
            %initial constraint in v0
            
            
            
            %package up the constraints
            coeff = [coeff_lie];
            con = [con_lie];
        end
        
        function objective = get_objective(obj, poly_var)
            %fetch the SOS objective to minimize
            objective = poly_var.gamma;
        end
        
        function [poly_val, func_eval] = recover_poly(obj, poly_var)
            %recover polynomials from SOS certificate
            
            [poly_val, func_eval] = recover_poly@location_sos_interface(obj, poly_var);
            
            poly_val.gamma = value(gamma);
            poly_val.beta = value(beta);
        end
        
    end
end


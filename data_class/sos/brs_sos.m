classdef brs_sos < reach_sos
    %BRS_SOS location handle for backwards reachable (BRS) set estimation 
    %given an external control process. Equivalent to control design for 
    %region of attraction maximization (ROA)
    
    %free terminal time is NOT included here.
    
    %this probably could have been designed through inheritence with 
    %brs_sos < reach_sos. oops.
    
    
    methods
        function obj = brs_sos(opts, mom_handle)
            %BRS_SOS Construct an instance of this class
            %opts:      loc_sos_options data structure
            %mom_handle:function handle for Lebesgue moments            
            obj@reach_sos(opts, mom_handle, 0);
%             obj.mom_handle = mom_handle;
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
                nonneg = [nonneg; poly.v];
            end
            
        end       
        
    end
    
end


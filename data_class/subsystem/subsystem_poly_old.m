classdef subsystem_poly_old < subsystem_interface
    %SUBSYSTEM_POLY A subsystem x'=f(t, x, w) of a dynamical system where
    %the uncertainty w is constrained to a polytope. This class is meant
    %may be applied to data-driven measure analysis
    
    properties
        
        %additional measures for polytopic uncertainty
        meas_pos = {};
        meas_neg = {};
        meas_slack = {};
        
        varnames = {'t', 'x'};
        
        f_box = [];     %polytopic decomposition of dynamics 
                        %[no input, input 1, input 2, ...]
                        
        A = [];
        b = [];
    end
    
    methods
    function obj = subsystem_poly_old(loc_supp, f, sys_id, loc_id)
            %SUBSYSTEM_POLY Construct a continuous system with polytopic
            %uncertainty. fill in information
            
            %process input
            if nargin < 3
                sys_id = 1;
            end
            
            if nargin < 4
                loc_id = [];
            end
            
            %superclass constructor
            obj@subsystem_interface(loc_supp, f, sys_id, loc_id, @meas_base);
%             obj.meas_type = @meas_uncertain;
            
            obj.dual = struct('v', 0, 'Lv', 0, 'Lv_box', 0, 'zeta', 0, 'nn', 0);
            
            Nconstraints = length(loc_supp.W);
            Nw = length(obj.vars.w);
            
                        
            obj.meas_pos = cell(Nw, 1);
            obj.meas_neg = cell(Nw, 1);
            obj.meas_comp = cell(Nconstraints, 1);
            for i = 1:Nw
                %box measure
                
                if Nw == 1
                    suffix_add = [];
                else
                    suffix_add = ['_', num2str(i)];
                end
                
                %positive input
                obj.meas_pos{i}  = obj.meas_def({'t', 'x'}, ['_pos', suffix_add], obj.supp);


                %negative input
                obj.meas_neg{i} = obj.meas_def({'t', 'x'}, ['_neg', suffix_add], obj.supp);
            end
            
            for i = 1:Nconstraints
                %negative input
                if Nw == 1
                    suffix_add = [];
                else
                    suffix_add = ['_', num2str(i)];
                end
                
                obj.meas_comp{i} = obj.meas_def({'t', 'x'}, ['_comp', suffix_add], obj.supp);
            end
            
            %process the dynamics f in terms of box dynamics
            %f_poly: [no input, input 1, input 2, ...]
            obj.f_poly = zeros(length(obj.f), Nw+1) * obj.vars.w(1);
            f0 = subs(obj.f, obj.vars.w, zeros(Nw, 1));                    
            obj.f_poly(:, 1) = f0;

            %each input channel at a time
            I = eye(Nw);

            for k = 1:Nw
                obj.f_poly(:, k+1) = subs(obj.f, obj.vars.w, I(:, k)) - f0;                        
            end                                
        end      
        
        function vars_out = get_vars(obj)
            %GET_VARS_BOX include box variables b
            vars_out = [obj.vars.t; obj.vars.x];
        end
        
        function vars_out = get_vars_box(obj)
            %GET_VARS_BOX include box variables b
            vars_out = [obj.vars.t; obj.vars.x; obj.vars.w];
        end
        
        
        function Ay = cons_liou(obj, d)
            function Ay = cons_liou(obj, d)
            %CONS_LIOU Liouville Equation includes an affine combination of
            %Lie derivatives (continuous systems only)
            
            if isempty(obj.vars.w)
                %no box inputs, simple to perform
                 Ay = obj.meas_occ.mom_lie(d, obj.get_vars, obj.f);
            else
                %non-trivial box inputs, more involved processing
                
                Nw = length(obj.vars.w);
                %base occupation measure (with no box disturbance)
%                 vars_
                Ay = obj.meas_occ.mom_lie(d, obj.get_vars, obj.f_poly(:, 1));
                
                %each input channel at a time
%                 I = eye(Nb);
                
                for k = 1:Nb
%                     fk = subs(obj.f, obj.vars.b, I(:, k)) - f0;
%                     Ay_curr = obj.meas_box{k}.mom_lie(d, obj.get_vars, obj.f_poly(:, k+1), 0);
                    Ay_pos = obj.meas_pos{k}.mom_lie(d, obj.get_vars, obj.f_poly(:, k+1), 0);
                    Ay_neg= obj.meas_neg{k}.mom_lie(d, obj.get_vars, obj.f_poly(:, k+1), 0);


                    %add contribution to lie derivative
                    Ay = Ay + Ay_pos - Ay_neg;
                end
                
            end
            
        function Ay = abscont_box(obj, d)
            %ABSCONT_BOX absolute continuity constraints of each box+complement with 
            %respect to the occupation measure
            Ay = [];
            
            %moments of each measure
            mom_occ = obj.meas_occ.mom_monom(d);
            for i = 1:length(obj.vars.b)
                mom_box  = obj.meas_box{i}.mom_monom(d);
                mom_comp = obj.meas_comp{i}.mom_monom(d);
                
                %absolute continuity constraint
                Ay_curr = -mom_occ + mom_box + mom_comp;
                Ay = [Ay; Ay_curr]; 
            end
        end  
            
        
        end
    end
    end

end


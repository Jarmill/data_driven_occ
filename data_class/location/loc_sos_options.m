classdef loc_sos_options
    %LOC_SOS_OPTIONS options data structure of pure SOS implementation of
    %system analysis. Meant only for polytopic input sets W
    %   this is a very bare-bones implementation
    %   only includes polytopic time-varying uncertainty in W: Aw >= b
    %   no switching, no box (box is a special type of polytopic)
    
    properties
        
        %% properties of run
        %terminal time   
        Tmax(1,1) double{mustBePositive}  = 5;           
        
        %% Variables and descriptors
        %variables defining sets (array of symbolic variables)
        
        t = []; %time
        x = []; %state
        
        verbose = 0; %solver output: https://yalmip.github.io/faq/runinsilent/
        
        %% support sets
        X = [];         %valid set
        X_init = [];    %initial set
        X_term = [];    %terminal set
        
        solver = 'mosek';
        
        %% Dynamics        
        f0 = []; %dynamics without uncertainty
        fw = {}; %elements are dynamics with each uncertainty
        
        %uncertainty set
        W = struct('A', [], 'b', []);
        
        %scale time from [0, Tmax] to [0, 1]?
        scale = 1;
        
    end
    
    methods
        
        %% support set getters
        function Tsupp = get_t_supp(obj)
            %get time support
            Tsupp = struct('ineq', [], 'eq', []);
            t = obj.t;
            if obj.scale
                Tsupp.ineq = t*(1-t);
            else
                Tsupp.ineq = t*(obj.Tmax - t);
            end
        end
        
        function Xsupp = get_X(obj)
            Xsupp = obj.X;
        end
        
        function allsupp = get_all_supp(obj)
            Tsupp = obj.get_t_supp();
            Xsupp = obj.get_X();
            
            allsupp = struct('ineq', [Tsupp.ineq; Xsupp.ineq], 'eq', Xsupp.eq);
        end
        
        function obj = set_box(obj, bounding_box)
            %set bounding box for X
            nx = length(obj.vars.x);
            [box, box_center, box_half] = box_process(nx, bounding_box);
            
            X_box = struct('ineq', box_half.^2 - (obj.vars.x - box_center).^2, 'eq', []);
            obj.X = X_box;                        
        end
        
%         function W = bounded_noise(t, x, epsilon)
%             %produce polytopic constraints describing system with bounded
%             %noise epsilon in the H infinity sense
%             
%         end
    end
    
end


classdef dist_sos_options < loc_sos_options
    %DIST_SOS_OPTIONS options data structure of pure SOS implementation of
    %system analysis. Meant only for polytopic input sets W
    %   this is a very bare-bones implementation
    %   only includes polytopic time-varying uncertainty in W: Aw >= b
    %   no switching, no box (box is a special type of polytopic)
    
    properties
        y = [];    %coordinates on unsafe set 
        X_unsafe = [];   %unsafe set
        c = [];    %distance function
    end
    
    methods
        function Xsupp = get_Xu(obj)
            Xsupp = obj.X_unsafe;
        end
        
        function Xsupp = get_X_Xu(obj)
            Xsupp = struct('ineq', [obj.X.ineq; obj.X_unsafe.ineq], ...
                'eq', [obj.X.eq; obj.X_unsafe.eq]);                       
        end
        
        function c_out = dist(obj)
            if isempty(obj.c)
                c_out = sum((obj.x-obj.y).^2);
            else
                c_out = obj.c;
            end
        end
    end
end
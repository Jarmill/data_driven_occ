function pfunc = polyval_func(poly, vars)
%return a function that evaluates the polynomial poly
%string-hacking to do this. For YALMIP SDPVAR only
%
%Input:
%   poly:   A polynomial in terms of sdpvar
%   vars:   The sdpvars of interest (important for namespace)
%
%Output:
%   pfunc:  A function that takes a point pt -> poly(pt) by evaluation

    
    if isempty(poly)
        pfunc = @(vars) [];
    else
        
        %find names of variables to replace
        s = sdisplay(poly);
        varnames_old = sdisplay(vars);
        varnames_new = cell(length(vars), 1);
        for i = 1:length(vars)
            varnames_new{i} = ['varnew(', num2str(i), ')'];
        end

        sjoin = ['[', sprintf('%s; ',s{1:end-1}),s{end}, ']'];
        var_pre = '@(varnew) ';
        sjoin_old = sjoin;
        
        for i = 1:length(vars)
            sjoin = strrep(sjoin, varnames_old{i}, varnames_new{i});
        end
        
        pfunc = eval([var_pre, sjoin]);
    end
    %
end
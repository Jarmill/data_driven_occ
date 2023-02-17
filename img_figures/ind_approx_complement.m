function [out] = ind_approx_complement(d, box, XT)
%IND_APPROX_COMPLEMENT approximate the indicator function on the sets XT
%from below, yielding an inner approximation
%univariate
x = sdpvar(1, 1);

%form sets
if length(box)==1
    box = box*[-1,1];
end

X = struct('ineq', (x-box(1))*(box(2)-x), 'eq', []);

Nset = size(XT, 1);
XT_cell = cell(Nset, 1);
for i = 1:Nset
    XT_cell{i} = struct('ineq', (x-XT(i, 1))*(XT(i,2)-x), 'eq', []);
end

% XTC = [[box(1), XT]; [XT, box(2)]];
XTC = [[box(1); XT(:, 2)], [XT(:, 1); box(2)]];
XTC_cell = cell(Nset+1, 1);
% XT_cell = cell(Nset, 1);
for i = 1:Nset+1
    XTC_cell{i} = struct('ineq', (x-XTC(i, 1))*(XTC(i,2)-x), 'eq', []);
end

leb = LebesgueBoxMom( d, box', 1);

%form polynomials
[w, cw] = polynomial(x, d);

objective = -leb'*cw;

cons = [];
coeff = [cw];

for i = 1:Nset+1
    [p0, cons0, coeff0] = constraint_psatz(-w, XTC_cell{i}, x, d);     
    cons = [cons; cons0:['Region Comp.', num2str(i)]];
    coeff = [coeff; coeff0];
end
for i = 1:Nset
    [pcurr, conscurr, coeffcurr] = constraint_psatz(1-w, XT_cell{i}, x, d);
    cons = [cons; conscurr:['Region ', num2str(i)]];
    coeff = [coeff; coeffcurr];
end

%sos program
sdp_opts = sdpsettings('solver', 'mosek', 'verbose', 1);
sdp_opts.sos.model = 2;
            
[sol, monom, Gram, residual] = solvesos(cons, objective, sdp_opts, coeff);
out = struct;
out.sol = sol;
out.monom = monom;
out.Gram = Gram;
out.residual = residual;


%extract polynomial w
% cw_eval = value(cw);
[cwe,mwe] = coefficients(w,x);
out.w = value(cwe)'*mwe;                 
out.w_eval = polyval_func(out.w, x);           
out.volume = -value(objective);



end


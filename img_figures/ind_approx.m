function [out] = ind_approx(d, box, XT)
%IND_APPROX approximate the indicator function on the sets XT
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

leb = LebesgueBoxMom( d, box', 1);

%form polynomials
[w, cw] = polynomial(x, d);

objective = leb'*cw;

[p0, cons0, coeff0] = constraint_psatz(w, X, x, d);     
cons = [cons0:'Nonneg'];
coeff = [cw; coeff0];

for i = 1:Nset
    [pcurr, conscurr, coeffcurr] = constraint_psatz(w-1, XT_cell{i}, x, d);
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
out.volume = value(objective);



end


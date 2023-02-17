function [putinar_curr, Gram] = make_Gram(g, q, M, monom)
    %create a Gram matrix for constraint g(vars)>=0
    %where q is the block size of the PSD constraint
    %mom: moment matrix indexing

    N = M.N;    
    Gram = sdpvar(q*N, q*N, 'symmetric');
    
    gmonom = g*monom;
    putinar_curr = zeros(q, q, 'like', sdpvar);
    for qi = 1:q
        ind_i = N *(qi-1) + (1:N);
        for qj = 1:q
            ind_j = N *(qj-1) + (1:N);
            G_accum = accumarray(M.M, reshape(Gram(ind_i, ind_j), [], 1));
            putinar_curr(qi, qj) = (G_accum'*gmonom);
       end
    end
end
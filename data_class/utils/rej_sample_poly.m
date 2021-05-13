function w = rej_sample_poly(A, b, box)
%REJ_SAMPLE_POLY Rejection sampling to generate uniformly-distributed 
%data inside a convex polytope
%
%go to hit and run later

[m, n] = size(A);

done = 0;
while ~done
    w = rand(n, 1)*(box(:, 2)-box(:, 1)) + box(:, 1);
    
    if all(A*w <= b)
        done = 1;
    end
end


end


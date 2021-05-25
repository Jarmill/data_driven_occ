function xdot = poly_3(t, x)
%POLY_3 Summary of this function goes here
%   Detailed explanation goes here
A_true = [-1 1 1; -1 0 -1; 0 1 -2];
B_true = [-1 0 -1;
          0 1 1;
          1 1 0]/2;

xdot =  A_true*x - B_true*(4*x.^3 - 3*x);

end


function R_f = Rfine(x,h_f)
%RCOARSE Summary of this function goes here
%   fine cheese-model response
if nargin < 2
    h_f = 1;
end
l_f = x(1);
f1 = x(2);
f2 = x(3);
w_f = x(4);
w_f1 = x(5);
w_f2 = x(6);
R_f(1) = (l_f*w_f*h_f);
R_f(2) = (l_f-f1)*w_f1*h_f;
R_f(3) = (l_f-f2)*w_f2*h_f;
end


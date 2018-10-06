function Rc = Rcoarse(x,h_c)
%RCOARSE Summary of this function goes here
%   coarse cheese model response
if nargin < 2
    h_c = 1;
end
l_c = x(1);
c = x(2);
w_c = x(3);
Rc(1) = l_c*w_c*h_c;
Rc(2) = (l_c-c)*w_c*h_c;
Rc(3) = (l_c-c)*w_c*h_c;
end


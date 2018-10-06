function [lc, opt_lc] = optim_sug(Raim,xin,dR)
%OPTIM_SUG Summary of this function goes here
%   optimise surrogate model w.r.t the length parameter
if nargin < 3
    dR = 0;
end

if nargin == 0
    % test-values
    xin = [11 2 0.8728];
    Raim = [10 10 10];
    dR = [1.399,-0.855,-0.855];
end

l_c = xin(1); c = xin(2); w_c = xin(3); 
% fun = @(l)abs((Raim - (Rcoarse([l, c, w_c])+ dR) ));
fun = @(l)abs((Raim -  Rsurrogate([l,c,w_c],dR) ));
options = optimoptions('fminimax','AbsoluteMaxObjectiveCount',3);
opt_lc = fminimax(fun,0,[],[],[],[],[],[],[],options);
lc = opt_lc;

%% quasi-newton method
% opt_lc = fminunc(fun,xin(1));
%% simulatedannealing
% opt_lc = simulannealbnd(fun,xin(1),[],[]);
%% pso
% opt_lc = particleswarm(fun,1,[],[]);

end



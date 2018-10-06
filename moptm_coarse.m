function [opt_lc, xe] = moptm_coarse(aim,xin,dr)
%OPTIMCOARSE Summary of this function goes here
%   optimize coarse model w.r.t the length parameter
if nargin < 3
    dr = 0;
end

if nargin == 0
    % test-values
    xin = [1 2 1];
    aim = [10 10 10];
end

fun = @(l)abs((aim - Rcoarse([l,xin(2),xin(3)])+dr));
% numel fun
% Minimize absolute values
%% minimax method
options = optimoptions('fminimax','AbsoluteMaxObjectiveCount',3);
opt_lc = fminimax(fun,xin(1),[],[],[],[],[],[],[],options);
xe = opt_lc - aim(1);

%% quasi-newton method
% opt_lc = fminunc(fun,xin(1));
%% simulatedannealing
% opt_lc = simulannealbnd(fun,xin(1),[],[]);
%% pso
% opt_lc = particleswarm(fun,1,[],[]);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              f_kstest
% author - Brad Rauscher (created 2024)
% 
% Calculates significance between each cell in 'X' using the two-sample,
% two-sided Kolmogorov-Smirnov test with significance value 'alpha'.
% 
% INPUTS: f_kstest(X, alpha)
%   X: cell array of distributions.
%   alpha: significance value
% 
% OUTPUTS:
%   h: hypothesis test values (1 if null hypothesis is rejected)
%   p: p-values
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [h, p] = f_kstest(X, alpha)

N = numel(X);

h = zeros(N);
p = zeros(N);

for hI = 1 : N
    for wI = 1 : N
        [h(hI, wI), p(hI, wI)] = kstest2(X{hI}, X{wI}, 'Alpha', alpha);
    end
end

end
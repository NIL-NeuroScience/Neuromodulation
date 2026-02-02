%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               f_lme
% author - Brad Rauscher (created 2024)
% 
% Calculates significance between each row of 'g1' and 'g2' after removing
% differences between mice using a linear mixed-effects model.
% 
% INPUTS: f_lme(MouseID, g1, g2, alpha)
%   MouseID: ID for each row of 'g1' and 'g2'
%   g1: group 1
%   g2: group 2
%   alpha: significance value
% 
% OUTPUTS:
%   h: hypothesis test values (1 if null hypothesis is rejected)
%   p: p-values
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [h, p] = f_lme(MouseID, g1, g2, alpha)

N = numel(MouseID);

col_MouseID = categorical([MouseID, MouseID])';
col_Group = categorical([repmat({'g1'}, 1, N), repmat({'g2'}, 1, N)])';
col_Data = cat(1, g1, g2);

G = size(col_Data, 2);

p = zeros(G, 1);

for i = 1 : G
    Data = col_Data(:, i);
    T = table(col_MouseID, col_Group, Data);
    
    lme = fitlme(T, 'Data ~ col_Group + (1|col_MouseID)');
    coeff = lme.Coefficients;
    p(i) = coeff.pValue(strcmp(coeff.Name, 'col_Group_g2'));
end

h = p < alpha;

end

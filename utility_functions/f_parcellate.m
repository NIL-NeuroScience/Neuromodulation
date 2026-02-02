%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                            f_parcellate
% author - Brad Rauscher (created 2024)
% 
% Calculates the spatial average of 'sig' within each mask in 'masks'.
% 
% INPUTS: f_parcellate(sig, masks)
%   sig: signal to parcellate (H x W x T)
%   masks: masks to parcellate sig (H x W x N)
% 
% OUTPUTS:
%   P: parcellated signal (T x N)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function P = f_parcellate(sig, masks)

% process 'masks'
masks = reshape(masks, [], size(masks, 3));
masks(isnan(masks)) = 0;

dim = size(sig);

sig = reshape(sig, [], dim(3))';
sig(isnan(sig)) = 0;

% add up 'sig' within each mask in 'masks'
P = sig * masks;

% normalize for mask size
P = P ./ sum(masks, 1);

end
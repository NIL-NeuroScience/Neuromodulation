%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           f_hemCoherence
% author - Brad Rauscher (created 2025)
% 
% Calculates coherence between each pixel in 'sig1' and 'sig2' using the
% coherencyc function in the Chronux toolbox. Averages the calculated
% coherence values within each 'mask'. 
% 
% INPUTS: f_hemCoherence(sig1, sig2, fs, masks, tapers)
%   sig1: first signal (H x W x T)
%   sig2: second signal  (H x W x T)
%   fs: sampling rate of signals (Hz)
%   mask: mask to average within (H x W)
%   tapers: multitaper parameters [time-bandwidth product, N tapers]
%       (default = [5, 9])
% 
% OUTPUTS:
%   C: coherence value(s)
%   phi: phase value(s) for C
%   f: frequency value(s) for C
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [C, phi, f] = f_hemCoherence(sig1, sig2, fs, mask, tapers)

% process inputs
dim = size(sig1);
nanIdx = isnan(mask);
sig1 = reshape(sig1, dim(1) * dim(2), dim(3));
sig2 = reshape(sig2, dim(1) * dim(2), dim(3));
sig1 = sig1(~nanIdx, :)';
sig2 = sig2(~nanIdx, :)';

% calculate coherence of each column in 'sig1' and 'sig2'
params = struct;
params.Fs = fs;
params.tapers = tapers;
params.trialave = 0;

[C, phi, ~, ~, ~, f] = coherencyc(sig1, sig2, params);

% calculate mean C and phi
C = mean(C, 2);
phi = mean(phi, 2);
f = f';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             f_hemSpectra
% author - Brad Rauscher (created 2025)
% 
% Calculates power spectral density for each pixel in 'sig' using the
% mtspectrumc function in the Chronux toolbox. Averages the calculated
% power spectral desnity values within 'mask'.
% 
% INPUTS: f_maskSpectra(sig, fs, masks, _)
%   sig: signal (H x W x T)
%   fs: sampling rate of signals (Hz)
%   mask: mask to average within (H x W)
%   tapers: multitaper parameters [time-bandwidth product, N tapers]
%       (default = [5, 9])
% 
% OUTPUTS:
%   S: spectral power density value(s)
%   f: frequency value(s) for S
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [S, f] = f_hemSpectra(sig, fs, fpass, mask, tapers)

% process inputs
dim = size(sig);
sig = reshape(sig, dim(1) * dim(2), dim(3));
nanIdx = isnan(mask);
sig = sig(~nanIdx, :);
sig = sig';

% calculate spectra of each column in 'sig'
params.Fs = fs;
params.fpass = fpass;
params.trialave = 0;
params.tapers = tapers;

[S, f] = mtspectrumc(sig, params);

% calculate mean S
S = mean(S, 2);
f = f';




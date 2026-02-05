%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             MAIN_plotFigures
% author - Brad Rauscher (created 2024)
% 
% Plots figures for data generated from main.m
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% setup paths
parentDir = f_path();
addPaths(parentDir);

dataPath = fullfile(parentDir, 'analysis');
savePath = fullfile(parentDir, 'results');

% load and organize data
[log, order, settings, Fig1, Fig2, Fig3, Behavior, Hb_model, ...
    unfiltered, shuffled, NE_reg, spectra, FC_fr, GRAB_FC] = ...
    f_organizeData(dataPath);

% load reference allen atlas
load(fullfile(parentDir, 'Figures/plot_types/refAllen.mat'));

% register brain masks to allen atlas
BM = f_regImages(settings.brain_mask, refParcellation, ...
    settings.allen_masks, 1);

% sort NE and ACh trials
M = numel(order);

NE_Idx = strcmp({order.GRAB}, 'GRAB_NE');
ACh_Idx = strcmp({order.GRAB}, 'GRAB_ACh');

M_NE = sum(NE_Idx);
M_ACh = sum(ACh_Idx);

% initialize variables
subAvg = struct;

plotBM = refBM;
plotBM(:, 1 : 300) = NaN;

%% plot Main figures
plotFig1
plotFig2
plotFig3

%% plot Extended Data figures
plotFigE1
plotFigE2
plotFigE3
plotFigE4
plotFigE5
plotFigE6
plotFigE7
plotFigE8
plotFigE9
plotFigE10

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EXTRA FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% adds paths of all directories under a given directory
function addPaths(mainDir)
    addDirs = genpath(mainDir);
    addDirs = strsplit(addDirs, pathsep);
    addDirs = addDirs(~cellfun(@(x) any(startsWith(strsplit(x, ...
        filesep), '.')), addDirs));
    addDirs = strjoin(addDirs, pathsep);
    addpath(addDirs);
end
%% organize data

parentDir = f_path();
addPaths(parentDir);
dataPath = [parentDir, '/analysis_final'];

[log, order, settings, Fig1, Fig2, Fig3, Behavior, Hb_model, ...
    unfiltered, shuffled, NE_reg, spectra, FC_fr, GRAB_FC] = ...
    f_organizeData(dataPath);

load(fullfile(parentDir, 'Figures/plot_types/refAllen.mat'));

M = numel(order);

BM = f_regImages(settings.brain_mask, refParcellation, ...
    settings.allen_masks, 1);
NE_Idx = strcmp({order.GRAB}, 'GRAB_NE');
ACh_Idx = strcmp({order.GRAB}, 'GRAB_ACh');

set(0, 'defaultfigurecolor', [1, 1, 1]);

M_NE = sum(NE_Idx);
M_ACh = sum(ACh_Idx);

savePath = fullfile(parentDir, 'results');

%% plot Main figures
plotFig1
plotFig2

%% plot Supplementary figures


%% EXTRA FUNCTIONS

function addPaths(mainDir)
    addDirs = genpath(mainDir);
    addDirs = strsplit(addDirs, pathsep);
    addDirs = addDirs(~cellfun(@(x) any(startsWith(strsplit(x, ...
        filesep), '.')), addDirs));
    addDirs = strjoin(addDirs, pathsep);
    addpath(addDirs);
end
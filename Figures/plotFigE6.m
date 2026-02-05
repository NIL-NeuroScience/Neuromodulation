%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              plotFigE6
% author - Brad Rauscher (created 2024)
% 
% Plots figure panels for Extended Data Figure 6. Must run 
% MAIN_plotFigures.m first.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

% load data
parentDir = f_path();
filePath = fullfile(parentDir, ...
    'Figures/results/Blocker/blocker_metadata.mat');
load(filePath);

load(fullfile(parentDir, 'Figures/plot_types/refAllen.mat'));
savePath = fullfile(parentDir, 'results');

fig_savePath = fullfile(savePath, 'ExtDataFig6');
[~, ~, ~] = mkdir(fig_savePath);

% order metadata
[E6_sessions, E6_ImagingOrder, E6_log] = f_orderMeta(metadata);

% create placeholder matrices

n_fr = numel(metadata.Thy1_302.d24_11_18.Run01.fr);
IRF_width = 101;

comb = struct;

comb.BM = cell(numel(E6_sessions), 1);
comb.allen = cell(numel(E6_sessions), 1);

comb.IRFx1.perf = cell(numel(E6_sessions), 1);
comb.IRFx1_var.perf = cell(numel(E6_sessions), 1);

comb.IRFx1.IRF = NaN(IRF_width, numel(E6_sessions));
comb.IRFx1_var.IRF = NaN(IRF_width, 12, numel(E6_sessions));

comb.SPC.HbT = NaN(n_fr, numel(E6_sessions));
comb.COH.Ca_HbT = NaN(n_fr, numel(E6_sessions));
comb.XC.Ca_HbT = NaN(101, numel(E6_sessions));

comb.COH.Ca_HbT_allen = NaN(n_fr, 12, numel(E6_sessions));

mice = fieldnames(metadata);

M_E6 = numel(mice);
runningIdx = 1;

for mIdx = 1 : M_E6
    date = fieldnames(metadata.(mice{mIdx}));
    d = numel(date);

    for dIdx = 1 : d
        run = fieldnames(metadata.(mice{mIdx}).(date{dIdx}));
        
        tmpDate = date{dIdx}(2 : end);
        tmpDate = strrep(tmpDate, '_', '-');

        r = numel(run);

        for rIdx = 1 : r
            interRun = ...
                numel(metadata.(mice{mIdx}).(date{dIdx}).(run{rIdx}));

            for irIdx = 1 : interRun

                runData = ...
                    metadata.(mice{mIdx}).(date{dIdx}).(run{rIdx})(irIdx);
                fr = runData.fr;

                comb.BM{runningIdx} = runData.brain_mask;

                comb.allen{runningIdx} = runData.parcellation;

                comb.IRFx1.perf{runningIdx} = runData.IRFx1.inv.perf;
                comb.IRFx1_var.perf{runningIdx} = runData.IRFx1.var.perf;

                comb.IRFx1.IRF(:, runningIdx) = runData.IRFx1.inv.IRF;

                comb.SPC.HbT(:, runningIdx) = runData.SPC.overall.HbT / ...
                    sum(runData.SPC.overall.HbT);
                comb.COH.Ca_HbT(:, runningIdx) = ...
                    runData.COH.overall.C.Ca_HbT;
                comb.XC.Ca_HbT(:, runningIdx) = runData.XC.overall.Ca_HbT;
                comb.COH.Ca_HbT_allen(:, :, runningIdx) = ...
                    runData.COH.allen.C.Ca_HbT';

                for i = 1 : 12
                    comb.IRFx1_var.IRF(:, i, runningIdx) = squeeze( ...
                        mean(runData.IRFx1.var.IRF .* ...
                        runData.allen(:, :, i), [1, 2], 'omitnan'));
                end

                runningIdx = runningIdx + 1;
            end
        end
    end
end

% register to allen atlas
comb.BM = f_regImages(comb.BM, refParcellation, comb.allen, 1);
E6_BM = comb.BM;

comb.IRFx1.perf = f_regImages(comb.IRFx1.perf, refParcellation, ...
    comb.allen, 0) .* E6_BM;
comb.IRFx1_var.perf = f_regImages(comb.IRFx1_var.perf, refParcellation, ...
    comb.allen, 0) .* E6_BM;

% combine mice
subAvg.FigE6 = struct;
subAvg.FigE6.IRFx1.perf = NaN([size(refBM), M_E6]);
subAvg.FigE6.IRFx1_var.perf = NaN([size(refBM), M_E6]);

subAvg.FigE6.IRFx1.IRF = NaN(IRF_width, M_E6);
subAvg.FigE6.IRFx1_var.IRF = NaN(IRF_width, 12, M_E6);

subAvg.FigE6.SPC.HbT = NaN(n_fr, M_E6);
subAvg.FigE6.COH.Ca_HbT = NaN(n_fr, M_E6);
subAvg.FigE6.COH.Ca_HbT_allen = NaN(n_fr, 12, M_E6);
subAvg.FigE6.XC.Ca_HbT = NaN(101, M_E6);

remove = 1 : 2;
for idx = 1 : M_E6
    tIdx = zeros(numel(E6_log), 1);
    tIdx(E6_ImagingOrder(idx).runs) = 1;
    tIdx = logical(tIdx .* (~ismember([E6_log.Run], remove))');
    subAvg.FigE6.IRFx1.perf(:, :, idx) = mean( ...
        comb.IRFx1.perf(:, :, tIdx) .* comb.BM(:, :, tIdx), 3, 'omitnan');
    subAvg.FigE6.IRFx1_var.perf(:, :, idx) = mean( ...
        comb.IRFx1_var.perf(:, :, tIdx) .* comb.BM(:, :, tIdx), 3, ...
        'omitnan');
    
    subAvg.FigE6.IRFx1.IRF(:, idx) = mean(comb.IRFx1.IRF(:, tIdx), 2);
    subAvg.FigE6.IRFx1_var.IRF(:, :, idx) = mean( ...
        comb.IRFx1_var.IRF(:, :, tIdx), 3);

    subAvg.FigE6.SPC.HbT(:, idx) = mean(comb.SPC.HbT(:, tIdx), 2);
    subAvg.FigE6.COH.Ca_HbT(:, idx) = mean(comb.COH.Ca_HbT(:, tIdx), 2);
    subAvg.FigE6.COH.Ca_HbT_allen(:, :, idx) = mean( ...
        comb.COH.Ca_HbT_allen(:, :, tIdx), 3);
    subAvg.FigE6.XC.Ca_HbT(:, idx) = mean(comb.XC.Ca_HbT(:, tIdx), 2);
end

subAvg.post = subAvg.FigE6;

remove = 3 : 8;
for idx = 1 : M_E6
    tIdx = zeros(numel(E6_log), 1);
    tIdx(E6_ImagingOrder(idx).runs) = 1;
    tIdx = logical(tIdx .* (~ismember([E6_log.Run], remove))');
    subAvg.FigE6.IRFx1.perf(:, :, idx) = mean( ...
        comb.IRFx1.perf(:, :, tIdx) .* comb.BM(:, :, tIdx), 3, 'omitnan');
    subAvg.FigE6.IRFx1_var.perf(:, :, idx) = mean( ...
        comb.IRFx1_var.perf(:, :, tIdx) .* comb.BM(:, :, tIdx), 3, ...
        'omitnan');
    
    subAvg.FigE6.IRFx1.IRF(:, idx) = mean(comb.IRFx1.IRF(:, tIdx), 2);
    subAvg.FigE6.IRFx1_var.IRF(:, :, idx) = mean( ...
        comb.IRFx1_var.IRF(:, :, tIdx), 3);

    subAvg.FigE6.SPC.HbT(:, idx) = mean(comb.SPC.HbT(:, tIdx), 2);
    subAvg.FigE6.COH.Ca_HbT(:, idx) = mean(comb.COH.Ca_HbT(:, tIdx), 2);
    subAvg.FigE6.COH.Ca_HbT_allen(:, :, idx) = mean( ...
        comb.COH.Ca_HbT_allen(:, :, tIdx), 3);
    subAvg.FigE6.XC.Ca_HbT(:, idx) = mean(comb.XC.Ca_HbT(:, tIdx), 2);
end

subAvg.baseline = subAvg.FigE6;

%% Extended Data Fig 6B

signalPath = fullfile(parentDir, 'Figures/results/Blocker/rawData.mat');
signals = load(signalPath);

Ca = signals.data1.rfp_HD(:, 5);
HbT = signals.data1.HbT(:, 5);
Pupil = signals.data1.Pupil;
Whisking = signals.data1.Whisking;
Movement = signals.data1.Accelerometer;

t = 0.1 : 0.1 : 600;

f = figure;
tiledlayout(4, 1);

nexttile;
hold on;
plot(t, Ca, color = c_Ca);
plot(t, HbT - 10, color = c_HbT);
plot([60, 60], [0, 10], '-k', lineWidth = 2);
xlim([50, 250]);
box off;

nexttile;
hold on;
plot(t, Pupil, color = c_pupil);
plot([60, 60], [0, 0.4], '-k', lineWidth = 2);
xlim([50, 250]);
box off;

nexttile;
plot(t, Whisking, color = [0, 0.7, 0.7]);
xlim([50, 250]);
box off;

nexttile;
plot(t, Movement, color = [0, 0, 0]);
xlim([50, 250]);
box off;

idx = 500 : 2500;

saveas(f, fullfile(fig_savePath, 'ExtDataFig6_B.svg'));

T = table(t(idx)', Ca(idx), HbT(idx), Pupil(idx), Whisking(idx), ...
    Movement(idx), VariableNames = {'Time', 'Ca', 'HbT', 'Pupil', ...
    'Whisking', 'Accelerometer'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig6_B.csv'));

%% Extended Data Fig 6C

signalPath = fullfile(parentDir, 'Figures/results/Blocker/rawData.mat');
signals = load(signalPath);

Ca = signals.data2.rfp_HD(:, 5);
HbT = signals.data2.HbT(:, 5);
Pupil = signals.data2.Pupil;
Whisking = signals.data2.Whisking;
Movement = signals.data2.Accelerometer;

t = 0.1 : 0.1 : 600;

f = figure;
tiledlayout(4, 1);

nexttile;
hold on;
plot(t, Ca, color = c_Ca);
plot(t, HbT - 10, color = c_HbT);
plot([260, 260], [0, 10], '-k', lineWidth = 2);
xlim([250, 450]);
box off;

nexttile;
hold on;
plot(t, Pupil, color = c_pupil);
plot([260, 260], [0, 0.4], '-k', lineWidth = 2);
xlim([250, 450]);
box off;

nexttile;
plot(t, Whisking, color = [0, 0.7, 0.7]);
xlim([250, 450]);
box off;

nexttile;
plot(t, Movement, color = [0, 0, 0]);
xlim([250, 450]);
box off;

idx = 2500 : 4500;

saveas(f, fullfile(fig_savePath, 'ExtDataFig6_C.svg'));

T = table(t(idx)', Ca(idx), HbT(idx), Pupil(idx), Whisking(idx), ...
    Movement(idx), VariableNames = {'Time', 'Ca', 'HbT', 'Pupil', ...
    'Whisking', 'Accelerometer'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig6_C.csv'));

%% Extended Data Fig 6D

sig = subAvg.baseline.COH.Ca_HbT;
meanSig1 = mean(sig(2 : end, :), 2);
error1 = std(sig(2 : end, :), 0, 2) / sqrt(M_E6);

sig = subAvg.post.COH.Ca_HbT;
meanSig2 = mean(sig(2 : end, :), 2);
error2 = std(sig(2 : end, :), 0, 2) / sqrt(M_E6);

tmpfr = fr(2 : end);

f = figure;
f_plotLineError(tmpfr, meanSig1, error1, ...
    color = c_Orange, ...
    log = 1);
f_plotLineError(tmpfr, meanSig2, error2, ...
    color = c_darkCyan, ...
    log = 1, ...
    xlim = [0.05, 5], ...
    xlabel = 'Frequency (Hz)', ...
    ylabel = 'C', ...
    ylim = [0, 0.8]);
set(gca, YScale = 'linear');

saveas(f, fullfile(fig_savePath, 'ExtDataFig6_D.svg'));
T = table(tmpfr', meanSig1, meanSig2, error1, error2, ...
    VariableNames = {'F (Hz)', 'mean_baseline', 'mean_PPA', ...
    'SEM_baseline', 'SEM_PPA'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig6_D.csv'));

%% Extended Data Fig 6E

[~, fr_Range] = min(abs(fr' - [0, 0.1]));

sig = squeeze(mean( ...
    subAvg.baseline.COH.Ca_HbT_allen(fr_Range(1) : fr_Range(2), :, :), 1));

f = figure;
f_plotAllenMap(mean(sig, 2), ...
    cmp = cmpinf, ...
    mask = plotBM, ...
    clim = [0, 1]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'ExtDataFig6_E_baseLow.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);


sig = squeeze(mean( ...
    subAvg.post.COH.Ca_HbT_allen(fr_Range(1) : fr_Range(2), :, :), 1));

f = figure;
f_plotAllenMap(mean(sig, 2), ...
    cmp = cmpinf, ...
    mask = plotBM, ...
    clim = [0, 1]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'ExtDataFig6_E_postLow.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);


sig = squeeze(mean( ...
    subAvg.post.COH.Ca_HbT_allen(fr_Range(1) : fr_Range(2), :, :), ...
    1)) - squeeze(mean( ...
    subAvg.baseline.COH.Ca_HbT_allen(fr_Range(1) : fr_Range(2), :, :), 1));

f = figure;
f_plotAllenMap(mean(sig, 2), ...
    cmp = cmpbbr, ...
    mask = plotBM, ...
    clim = 0.3 * [-1, 1]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'ExtDataFig6_E_diffLow.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);


[~, fr_Range] = min(abs(fr' - [0.1, 0.5]));

sig = squeeze(mean( ...
    subAvg.baseline.COH.Ca_HbT_allen(fr_Range(1) : fr_Range(2), :, :), 1));

f = figure;
f_plotAllenMap(mean(sig, 2), ...
    cmp = cmpinf, ...
    mask = plotBM, ...
    clim = [0, 1]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'ExtDataFig6_E_baseHigh.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);


sig = squeeze(mean( ...
    subAvg.post.COH.Ca_HbT_allen(fr_Range(1) : fr_Range(2), :, :), 1));

f = figure;
f_plotAllenMap(mean(sig, 2), ...
    cmp = cmpinf, ...
    mask = plotBM, ...
    clim = [0, 1]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'ExtDataFig6_E_postHigh.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);


sig = squeeze(mean( ...
    subAvg.post.COH.Ca_HbT_allen(fr_Range(1) : fr_Range(2), :, :), ...
    1)) - squeeze(mean( ...
    subAvg.baseline.COH.Ca_HbT_allen(fr_Range(1) : fr_Range(2), :, :), 1));

f = figure;
f_plotAllenMap(mean(sig, 2), ...
    cmp = cmpbbr, ...
    mask = plotBM, ...
    clim = 0.3 * [-1, 1]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'ExtDataFig6_E_diffHigh.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);

%% Extended Data Fig 6F

sig = subAvg.baseline.XC.Ca_HbT;
meanSig1 = mean(sig, 2);
error1 = std(sig, 0, 2) / sqrt(M_E6);

sig = subAvg.post.XC.Ca_HbT;
meanSig2 = mean(sig, 2);
error2 = std(sig, 0, 2) / sqrt(M_E6);

t = 5:-0.1:-5;

f = figure;
f_plotLineError(t, meanSig1, error1, ...
    color = c_Orange);
f_plotLineError(t, meanSig2, error2, ...
    color = c_darkCyan, ...
    xlim = [-5, 5], ...
    ylim = [-0.2, 0.4], ...
    xlabel = 'Time (s)', ...
    ylabel = 'r');

saveas(f, fullfile(fig_savePath, 'ExtDataFig6_F.svg'));
T = table(t', meanSig1, meanSig2, error1, error2, ...
    VariableNames = {'Time', 'mean_baseline', 'mean_PPA', 'SEM_baseline', ...
    'SEM_PPA'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig6_F.csv'));

%% Extended Data Fig 6G

sig = squeeze(mean(subAvg.baseline.IRFx1_var.IRF(:, [4, 5], :), 2));
meanSig1 = mean(sig, 2);
error1 = std(sig, 0, 2) / sqrt(M_E6);

t = 0 : 0.1 : 10;

f = figure;
f_plotLineError(t, meanSig1, error1, ...
    color = c_Orange);

sig = squeeze(mean(subAvg.post.IRFx1_var.IRF(:, [4, 5], :), 2));
meanSig2 = mean(sig, 2);
error2 = std(sig, 0, 2) / sqrt(M_E6);

f_plotLineError(t, meanSig2, error2, ...
    color = c_darkCyan, ...
    xlim = [0, 7], ...
    xlabel = 'Time (s)', ...
    ylabel = 'a.u.', ...
    ylim = [-0.04, 0.1]);

plot([0, 7], [0, 0], '-k', lineWidth = 2);

saveas(f, fullfile(fig_savePath, 'ExtDataFig6_G1.svg'));
T = table(t', meanSig1, meanSig2, error1, error2, ...
    VariableNames = {'Time', 'mean_baseline', 'mean_PPA', ...
    'SEM_baseline', 'SEM_PPA'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig6_G1.csv'));


sig = squeeze(mean(subAvg.baseline.IRFx1_var.IRF(:, 2, :), 2));
meanSig1 = mean(sig, 2);
error1 = std(sig, 0, 2) / sqrt(M_E6);

t = 0 : 0.1 : 10;

f = figure;
f_plotLineError(t, meanSig1, error1, ...
    color = c_Orange);

sig = squeeze(mean(subAvg.post.IRFx1_var.IRF(:, 2, :), 2));
meanSig2 = mean(sig, 2);
error2 = std(sig, 0, 2) / sqrt(M_E6);

f_plotLineError(t, meanSig2, error2, ...
    color = c_darkCyan, ...
    xlim = [0, 7], ...
    xlabel = 'Time (s)', ...
    ylabel = 'a.u.', ...
    ylim = [-0.04, 0.1]);

plot([0,7],[0,0],'-k',lineWidth=2);

saveas(f, fullfile(fig_savePath, 'ExtDataFig6_G2.svg'));
T = table(t', meanSig1, meanSig2, error1, error2, ...
    VariableNames = {'Time', 'mean_baseline', 'mean_PPA', ...
    'SEM_baseline', 'SEM_PPA'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig6_G2.csv'));

%% Extended Data Fig 6H

f = figure;
f_plotMap(mean(subAvg.baseline.IRFx1_var.perf, 3, 'omitnan') .* plotBM, ...
    cmp = cmpvir, ...
    clim = [0, 1])
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'ExtDataFig6_H_base.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);

f = figure;
f_plotMap(mean(subAvg.post.IRFx1_var.perf, 3, 'omitnan') .* plotBM, ...
    cmp = cmpvir, ...
    clim = [0, 1])
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'ExtDataFig6_H_post.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);

f = figure;
f_plotMap(mean(subAvg.post.IRFx1_var.perf, 3, 'omitnan') - ...
    mean(subAvg.baseline.IRFx1_var.perf, 3, 'omitnan') .* plotBM, ...
    cmp = cmpbbr, ...
    clim = 0.5 * [-1, 1])
colorbar off;

exportgraphics(f, fullfile(fig_savePath,'ExtDataFig6_H_diff.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);

%% Extended Data Fig 6I

tmpMask = refParcellation.Masks(:, :, 3, 2) .* refBM;
tmpMask(tmpMask == 0) = NaN;

post = squeeze(mean( ...
    subAvg.post.IRFx1_var.perf .* tmpMask, [1, 2], 'omitnan'));
baseline = squeeze(mean( ...
    subAvg.baseline.IRFx1_var.perf .* tmpMask, [1, 2], 'omitnan'));

barData = {};
barData{1} = baseline;
barData{2} = post;

f = figure;
[dataMean, dataSEM] = f_plotBar(barData, ...
    colors = [c_Orange; c_darkCyan], ...
    ylabel = 'r', ...
    ylim = [0, 1]);

[h, p] = f_kstest(barData, 0.05);

saveas(f, fullfile(fig_savePath, 'ExtDataFig6_I.svg'));
T = table(mice, baseline, post, ...
    VariableNames = {'Mouse', 'baseline', 'post'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig6_I.csv'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EXTRA FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% order blocker metadata
function [Order,sorted,log] = f_orderMeta(metadata)
    mice = fieldnames(metadata);
    Order = struct;
    log = struct;
    m = numel(mice);
    for mIdx = 1:m
        date = fieldnames(metadata.(mice{mIdx}));
        d = numel(date);
    
        for dIdx = 1:d
            run = fieldnames(metadata.(mice{mIdx}).(date{dIdx}));
            dateSingle = date{dIdx};
            dateSingle = strrep(dateSingle(2:end),'_','-');

            r = numel(run);
    
            for rIdx = 1:r
                runNum = metadata.(mice{mIdx}).(date{dIdx}).(run{rIdx}).Run;
                interRuns = numel(metadata.(mice{mIdx}).(date{dIdx}).(run{rIdx}));
                GRAB = metadata.(mice{mIdx}).(date{dIdx}).(run{rIdx})(1).GRAB;
                for interIdx = 1:interRuns
                    Order(end+1).mouse = mice{mIdx};
                    log(end+1).mouse = mice{mIdx};
                    log(end).date = dateSingle;
                    log(end).Run = runNum;
                    log(end).interRun = interIdx;
                    log(end).GRAB = GRAB;
                end
            end
        end
    end
    Order(1) = [];
    log(1) = [];
    sorted = struct;
    for mIdx = 1:m
        sorted(mIdx).mouse = mice{mIdx};
        sorted(mIdx).runs = find(strcmp({Order.mouse},mice{mIdx}));
    end
end
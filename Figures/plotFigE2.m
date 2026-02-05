%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              plotFigE2
% author - Brad Rauscher (created 2024)
% 
% Plots figure panels for Extended Data Figure 2. Must run 
% MAIN_plotFigures.m first.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

% calculate subject averages
NE_order = order(NE_Idx);

tmp_data.Beh.R_rfp_HD_gfp_HD = ...
    f_regImages(Behavior.R.rfp_HD_low_gfp_HD_low, refParcellation, ...
    settings.allen_masks, 0) .* BM;
tmp_data.Beh.NE_IRF_perf = ...
    f_regImages(Behavior.NE_IRF_perf, refParcellation, ...
    settings.allen_masks, 0) .* BM;
tmp_data.Beh.SPG_rfp_HD = f_adjust_SPG(Behavior.SPG.rfp_HD, 4097, 1);
tmp_data.Beh.SPG_gfp_HD = f_adjust_SPG(Behavior.SPG.gfp_HD, 4097, 1);
tmp_data.Beh.SPG_HbT = f_adjust_SPG(Behavior.SPG.HbT, 4097, 1);
tmp_data.Beh.COH_rfp_HD_gfp_HD = ...
    f_adjust_SPG(Behavior.COH.rfp_HD_gfp_HD, 4097, 0);
tmp_data.Beh.COH_rfp_HD_HbT = ...
    f_adjust_SPG(Behavior.COH.rfp_HD_HbT, 4097, 0);
tmp_data.Beh.COH_gfp_HD_HbT = ...
    f_adjust_SPG(Behavior.COH.gfp_HD_HbT, 4097, 0);
tmp_data.Beh.PHI_rfp_HD_gfp_HD = ...
    f_adjust_SPG(Behavior.PHI.rfp_HD_gfp_HD, 4097, 0);
tmp_data.Beh.PHI_rfp_HD_HbT = ...
    f_adjust_SPG(Behavior.PHI.rfp_HD_HbT, 4097, 0);
tmp_data.Beh.PHI_gfp_HD_HbT = ...
    f_adjust_SPG(Behavior.PHI.gfp_HD_HbT, 4097, 0);

subAvg.FigE2.XC_gfp_HD_HbT = NaN(201, M);
subAvg.FigE2.XC_rfp_HD_HbT = NaN(201, M);
subAvg.FigE2.XC_rfp_HD_gfp_HD = NaN(201, M);
subAvg.FigE2.R_rfp_HD_gfp_HD = NaN(500, 600, M);
subAvg.FigE2.R_behavior = NaN(6, 6, M);
subAvg.FigE2.SPG_rfp_HD = NaN(4097, M);
subAvg.FigE2.SPG_gfp_HD = NaN(4097, M);
subAvg.FigE2.SPG_HbT = NaN(4097, M);
subAvg.FigE2.COH_rfp_HD_gfp_HD = NaN(4097, M);
subAvg.FigE2.COH_rfp_HD_HbT = NaN(4097, M);
subAvg.FigE2.COH_gfp_HD_HbT = NaN(4097, M);
subAvg.FigE2.PHI_rfp_HD_gfp_HD = NaN(4097, M);
subAvg.FigE2.PHI_rfp_HD_HbT = NaN(4097, M);
subAvg.FigE2.PHI_gfp_HD_HbT = NaN(4097, M);
subAvg.FigE2.NE_IRF_perf = NaN(500,600, M);
subAvg.FigE2.NE_IRF = NaN(151, M);
subAvg.FigE2.GRAB_conn = NaN(12, 12, M);

for i = 1:M
    subAvg.FigE2.XC_gfp_HD_HbT(:, i) = ...
        mean(cat(2, Behavior.XC.gfp_HD_HbT{order(i).Runs}), 2);
    subAvg.FigE2.XC_rfp_HD_HbT(:, i) = ...
        mean(cat(2, Behavior.XC.rfp_HD_HbT{order(i).Runs}), 2);
    subAvg.FigE2.XC_rfp_HD_gfp_HD(:, i) = ...
        mean(cat(2, Behavior.XC.rfp_HD_gfp_HD{order(i).Runs}), 2);
    subAvg.FigE2.R_rfp_HD_gfp_HD(:, :, i) = mean( ...
        tmp_data.Beh.R_rfp_HD_gfp_HD(:, :, order(i).Runs), 3, 'omitnan');
    subAvg.FigE2.R_behavior(:, :, i) = mean( ...
        cat(3, Behavior.R.signals{order(i).Runs}), 3,'omitnan');
    subAvg.FigE2.SPG_rfp_HD(:, i) = mean( ...
        cat(2, tmp_data.Beh.SPG_rfp_HD{order(i).Runs}), 2);
    subAvg.FigE2.SPG_gfp_HD(:, i) = mean( ...
        cat(2, tmp_data.Beh.SPG_gfp_HD{order(i).Runs}), 2);
    subAvg.FigE2.SPG_HbT(:, i) = mean( ...
        cat(2, tmp_data.Beh.SPG_HbT{order(i).Runs}), 2);
    subAvg.FigE2.COH_rfp_HD_gfp_HD(:, i) = mean( ...
        cat(2, tmp_data.Beh.COH_rfp_HD_gfp_HD{order(i).Runs}), 2);
    subAvg.FigE2.COH_rfp_HD_HbT(:, i) = mean( ...
        cat(2, tmp_data.Beh.COH_rfp_HD_HbT{order(i).Runs}), 2);
    subAvg.FigE2.COH_gfp_HD_HbT(:, i) = mean( ...
        cat(2, tmp_data.Beh.COH_gfp_HD_HbT{order(i).Runs}), 2);
    subAvg.FigE2.PHI_rfp_HD_gfp_HD(:, i) = mean( ...
        cat(2, tmp_data.Beh.PHI_rfp_HD_gfp_HD{order(i).Runs}), 2);
    subAvg.FigE2.PHI_rfp_HD_HbT(:, i) = mean( ...
        cat(2, tmp_data.Beh.PHI_rfp_HD_HbT{order(i).Runs}), 2);
    subAvg.FigE2.PHI_gfp_HD_HbT(:, i) = mean( ...
        cat(2, tmp_data.Beh.PHI_gfp_HD_HbT{order(i).Runs}), 2);
    subAvg.FigE2.NE_IRF_perf(:, :, i) = mean( ...
        tmp_data.Beh.NE_IRF_perf(:, :, order(i).Runs), 3, 'omitnan');
    subAvg.FigE2.GRAB_conn(:, :, i) = mean( ...
        cat(3, GRAB_FC.FC_detrend{order(i).Runs}), 3);
    try 
        subAvg.FigE2.NE_IRF(:, i) = mean( ...
        cat(2, Behavior.NE_IRF_IRF{order(i).Runs}), 2);
    end
end

fr = Behavior.SPG.fr;

fig_savePath = fullfile(savePath, 'ExtDataFig2');
[~, ~, ~] = mkdir(fig_savePath);

%% Extended Data Fig 2A

f = figure;
meanSig1 = mean(subAvg.FigE2.SPG_rfp_HD, 2);
error1 = std(subAvg.FigE2.SPG_rfp_HD, 0, 2) / sqrt(M);
f_plotLineError(fr, meanSig1, error1, ...
    color = c_Ca, ...
    log = 1);

meanSig2 = mean(subAvg.FigE2.SPG_gfp_HD(:, ACh_Idx), 2);
error2 = std(subAvg.FigE2.SPG_gfp_HD(:, ACh_Idx), 0, 2) / sqrt(M_ACh);
f_plotLineError(fr, meanSig2, error2, ...
    color = c_Orange, ...
    log = 1);

meanSig3 = mean(subAvg.FigE2.SPG_gfp_HD(:, NE_Idx), 2);
error3 = std(subAvg.FigE2.SPG_gfp_HD(:, NE_Idx), 0, 2) / sqrt(M_NE);
f_plotLineError(fr, meanSig3, error3, ...
    color = c_GRAB, ...
    log = 1);

meanSig4 = mean(subAvg.FigE2.SPG_HbT, 2);
error4 = std(subAvg.FigE2.SPG_HbT, 0, 2) / sqrt(M);
f_plotLineError(fr, meanSig4, error4, ...
    color = c_HbT, ...
    log = 1);

xlim([0.05, 5]);
xlabel('F (Hz)');
ylabel('Normalized Power');

set(gcf, 'Renderer','painters');
saveas(f, fullfile(fig_savePath, 'ExtDataFig2_A.svg'));
T = table(fr, meanSig1, meanSig2, meanSig3, meanSig4, error1, error2, ...
    error3, error4, VariableNames = {'F (Hz)', 'mean_Ca', 'mean_ACh', ...
    'mean_NE', 'mean_HbT', 'SEM_Ca', 'SEM_ACh', 'SEM_NE', 'SEM_HbT'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig2_A.csv'));

%% Extended Data Fig 2B

t = 10 : -0.1 : -10;

f = figure;
meanSig1 = mean(subAvg.FigE2.XC_rfp_HD_gfp_HD(:, NE_Idx), 2);
error1 = std(subAvg.FigE2.XC_rfp_HD_gfp_HD(:, NE_Idx), 0, 2) / sqrt(M_NE);
f_plotLineError(t, meanSig1, error1, ...
    color = c_GRAB);

meanSig2 = mean(subAvg.FigE2.XC_rfp_HD_gfp_HD(:, ACh_Idx), 2);
error2 = std(subAvg.FigE2.XC_rfp_HD_gfp_HD(:, ACh_Idx), 0, 2) / ...
    sqrt(M_ACh);
f_plotLineError(t, meanSig2, error2, ...
    color = c_Orange, ...
    xlim = [-5, 5], ...
    xlabel = 'Time (s)', ...
    ylabel = 'r', ...
    title = 'x vs. Ca^2^+', ...
    legend = {'NE', 'ACh'});

saveas(f, fullfile(fig_savePath, 'ExtDataFig2_B.svg'));
T = table(t', meanSig1, meanSig2, error1, error2, ...
    VariableNames = {'Lag', 'mean_NE', 'mean_ACh', 'SEM_NE', 'SEM_ACh'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig2_B.csv'));

%% Extended Data Fig 2C

R_Beh = mean(subAvg.FigE2.R_behavior, 3);
gfp_R = mean(subAvg.FigE2.R_behavior(:, :, ACh_Idx), 3);
R_Beh(2, :) = gfp_R(2, :);
R_Beh(:, 2) = gfp_R(2, :);
gfp_R = mean(subAvg.FigE2.R_behavior(:, :, NE_Idx), 3);
R_Beh(4 : 7, :) = R_Beh(3 : 6, :);
R_Beh(:, 4 : 7) = R_Beh(:, 3 : 6);
R_Beh(3, [1, 3, 4, 5, 6, 7]) = gfp_R(2, :);
R_Beh([1, 3, 4, 5, 6, 7], 3) = gfp_R(2, :);
R_Beh(2, 3) = 0;
R_Beh(3, 2) = 0;

barData = {};
barData{1} = squeeze(subAvg.FigE2.R_behavior(2, 4, ACh_Idx));
barData{2} = squeeze(subAvg.FigE2.R_behavior(2, 5, ACh_Idx));
barData{3} = squeeze(subAvg.FigE2.R_behavior(2, 6, ACh_Idx));
barData{4} = squeeze(subAvg.FigE2.R_behavior(2, 4, NE_Idx));
barData{5} = squeeze(subAvg.FigE2.R_behavior(2, 5, NE_Idx));
barData{6} = squeeze(subAvg.FigE2.R_behavior(2, 6, NE_Idx));
barData{7} = squeeze(subAvg.FigE2.R_behavior(1, 4, :));
barData{8} = squeeze(subAvg.FigE2.R_behavior(1, 5, :));
barData{9} = squeeze(subAvg.FigE2.R_behavior(1, 6, :));
barData{10} = squeeze(subAvg.FigE2.R_behavior(3, 4, :));
barData{11} = squeeze(subAvg.FigE2.R_behavior(3, 5, :));
barData{12} = squeeze(subAvg.FigE2.R_behavior(3, 6, :));

f = figure;
[dataMean, dataSEM] = f_plotBar(barData, ...
    colors = repmat([c_pupil; 0, 0.7, 0.7; 0, 0, 0], 4, 1), ...
    legend = {'Pupil diameter', 'Whisking', 'Movement'}, ...
    ylabel = 'r', ...
    title = 'Model Performance Comparison', ...
    ylim = [0, 0.8]);

saveas(f, fullfile(fig_savePath, 'ExtDataFig2_C1.svg'));

labels = {'ACh_pupil', 'ACh_whisking', 'ACh_movement', ...
    'NE_pupil', 'NE_whisking', 'NE_movement', ...
    'Ca_pupil', 'Ca_whisking', 'Ca_movement', ...
    'HbT_pupil', 'HbT_whisking', 'HbT_movement'};
T = table({order.Mouse}', VariableNames = {'Mouse'});

for i = 1 : 3
    tmp_data = NaN(numel({order.Mouse}), 1);
    tmp_data(ACh_Idx) = barData{i};
    T.(labels{i}) = tmp_data;
end
for i = 4 : 6
    tmp_data = NaN(numel({order.Mouse}), 1);
    tmp_data(NE_Idx) = barData{i};
    T.(labels{i}) = tmp_data;
end
for i = 7 : 12
    T.(labels{i}) = barData{i};
end

writetable(T, fullfile(fig_savePath, 'ExtDataFig2_C1.csv'));

f = figure;
f_plotMap(R_Beh, ...
    cmp = cmpbbr, ...
    clim = [-1, 1]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'ExtDataFig2_C2.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);
writetable(table(R_Beh), fullfile(fig_savePath, 'ExtDataFig2_C2.csv'));

%% Extended Data Fig 2D

f = figure;
f_plotMap(mean(subAvg.FigE2.R_rfp_HD_gfp_HD(:, :, ACh_Idx), 3, ...
    'omitnan') .* plotBM, ...
    cmp = cmpvir, ...
    clim = [0, 1]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'ExtDataFig2_D1.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);

barData = {};
barData{1} = squeeze(mean(subAvg.FigE2.R_rfp_HD_gfp_HD(:, :, ACh_Idx) ...
    .* plotBM, [1, 2], 'omitnan'));

f = figure;
[dataMean, dataSEM] = f_plotBar(barData, ...
    colors = c_darkCyan, ...
    ylabel = 'r', ...
    ylim = [0, 1]);

saveas(f, fullfile(fig_savePath, 'ExtDataFig2_D2.svg'));
T = table({order(ACh_Idx).Mouse}', barData{1}, ...
    VariableNames = {'Mouse', 'r'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig2_D2.csv'));

%% Extended Data Fig 2E

f = figure;
meanSig1 = mean(subAvg.FigE2.COH_rfp_HD_gfp_HD(: ,ACh_Idx), 2);
error1 = std(subAvg.FigE2.COH_rfp_HD_gfp_HD(:, ACh_Idx), 0, 2) / ...
    sqrt(M_ACh);
f_plotLineError(fr, meanSig1, error1, ...
    color = c_Orange, ...
    log = 1, ...
    ylim = [0, 1], ...
    ylabel = 'Coherence');
set(gca, YScale = 'linear');

yyaxis right;
meanSig2 = mean(subAvg.FigE2.PHI_rfp_HD_gfp_HD(:, ACh_Idx), 2);
error2 = std(subAvg.FigE2.PHI_rfp_HD_gfp_HD(:, ACh_Idx), 0, 2) / ...
    sqrt(M_ACh);
f_plotLineError(fr, meanSig2, error2, ...
    color=[0,0,0], ...
    log=1, ...
    ylim = pi * [-1, 1], ...
    xlim = [0.05, 5], ...
    xlabel = 'F (Hz)', ...
    ylabel = 'Phi (rad)');
set(gca, YScale = 'linear');

saveas(f, fullfile(fig_savePath, 'ExtDataFig2_E.svg'));
T = table(fr, meanSig1, meanSig2, error1, error2, ...
    VariableNames = {'F (Hz)', 'mean_coherence', 'mean_phase', ...
    'SEM_coherence', 'SEM_phase'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig2_E.csv'));

%% Extended Data Fig 2F

run = 67;

Ca = Fig1.Ca_allen{run};
ACh = GRAB_FC.GRAB{run};

t = 0.1 : 0.1 : 600;

f = figure;
hold on;
plot(t, 2 * Ca(:, 5), color = c_Ca);
plot(t, ACh(:, 5) - 10, color = c_Orange);
plot(t, 2 * Ca(:, 3) - 20, color = c_Ca);
plot(t, ACh(:, 3) - 30, color = c_Orange);
plot([110, 110], [0, 10], '-k');
xlim([100, 220]);

idx = 1000 : 2200;

saveas(f, fullfile(fig_savePath, 'ExtDataFig2_F.svg'));
T = table(t(idx)', Ca(idx, 5), ACh(idx, 5), Ca(idx, 3), ACh(idx, 3), ...
    VariableNames = {'Time', 'Ca_SSpll', 'ACh_SSpll', 'Ca_SSpbfd', ...
    'ACh_SSpbfd'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig2_F.csv'));

%% Extended Data Fig 2G

run = 67;

Ca = Fig1.Ca_allen{run};
ACh = GRAB_FC.GRAB{run};

lm = fitlm(ACh(:, 5), Ca(:, 5));
lm = table2array(lm.Coefficients);
f_corr(ACh(:, 5), Ca(:, 5), 1)

f = figure;
hold on;
scatter(ACh(:, 5), Ca(:, 5), 75, 'filled', ...
    MarkerFaceAlpha = 0.1, ...
    MarkerFaceColor = c_Orange);
plot([-30, 30], [-30, 30] * lm(2, 1) + lm(1, 1), '-k');
ylim([-10, 15]);
xlim([-30, 30]);

saveas(f, fullfile(fig_savePath, 'ExtDataFig2_G1.svg'));
T = table(Ca(idx, 5), ACh(idx, 5), ...
    VariableNames = {'Ca_SSpll', 'ACh_SSpll'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig2_G1.csv'));


run = 129;
runData = load(...
    'sub-Thy1-332_ses-25-02-04_run-01_irun-01_behavior+ophys.mat');

Ca = runData.Fig1.Ca_allen;
NE = GRAB_FC.GRAB{run};

lm = fitlm(NE(:, 5), Ca(:, 5));
lm = table2array(lm.Coefficients);
f_corr(NE(:, 5), Ca(:, 5), 1)

f = figure;
hold on;
scatter(NE(:, 5), Ca(:, 5), 75, 'filled', ...
    MarkerFaceAlpha = 0.1, ...
    MarkerFaceColor = c_GRAB);
plot([-6, 8], [-6, 8] * lm(2, 1) + lm(1, 1), '-k');
ylim([-10, 15]);
xlim([-6, 8]);

saveas(f, fullfile(fig_savePath, 'ExtDataFig2_G2.svg'));
T = table(Ca(idx, 5), NE(idx, 5), ...
    VariableNames = {'Ca_SSpll', 'NE_SSpll'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig2_G2.csv'));

%% Extended Data Fig 2H

f = figure;
f_plotMap(mean(subAvg.FigE2.R_rfp_HD_gfp_HD(:, :, NE_Idx), 3, ...
    'omitnan') .* plotBM, ...
    cmp = cmpvir, ...
    clim = [0, 1]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'ExtDataFig2_H1.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);

barData = {};
barData{1} = squeeze(mean(subAvg.FigE2.R_rfp_HD_gfp_HD(:, :, NE_Idx) .* ...
    plotBM, [1, 2], 'omitnan'));

f = figure;
[dataMean, dataSEM] = f_plotBar(barData, ...
    colors = c_darkCyan, ...
    ylabel = 'r', ...
    ylim = [0, 1]);

saveas(f, fullfile(fig_savePath, 'ExtDataFig2_H2.svg'));
T = table({order(NE_Idx).Mouse}', barData{1}, ...
    VariableNames = {'Mouse', 'r'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig2_H2.csv'));

%% Extended Data Fig 2I

f = figure;
meanSig1 = mean(subAvg.FigE2.COH_rfp_HD_gfp_HD(:, NE_Idx), 2);
error1 = std(subAvg.FigE2.COH_rfp_HD_gfp_HD(:, NE_Idx), 0, 2) / sqrt(M_NE);
f_plotLineError(fr, meanSig1, error1, ...
    color = c_GRAB, ...
    log = 1, ...
    ylim = [0, 1], ...
    ylabel = 'Coherence');
set(gca, YScale = 'linear');

yyaxis right;
meanSig2 = mean(subAvg.FigE2.PHI_rfp_HD_gfp_HD(:, NE_Idx), 2);
error2 = std(subAvg.FigE2.PHI_rfp_HD_gfp_HD(:, NE_Idx), 0, 2) / sqrt(M_NE);
f_plotLineError(fr, meanSig2, error2, ...
    color = [0, 0, 0], ...
    log = 1, ...
    ylim = pi * [-1, 1], ...
    xlim = [0.05, 5], ...
    xlabel = 'F (Hz)', ...
    ylabel = 'Phi (rad)');
set(gca, YScale = 'linear');

saveas(f, fullfile(fig_savePath, 'ExtDataFig2_I.svg'));
T = table(fr, meanSig1, meanSig2, error1, error2, ...
    VariableNames = {'F (Hz)', 'mean_coherence', 'mean_phase', ...
    'SEM_coherence', 'SEM_phase'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig2_I.csv'));

%% Extended Data Fig 2J

meanSig = mean(subAvg.FigE2.NE_IRF(:, NE_Idx), 2);
error = std(subAvg.FigE2.NE_IRF(:, NE_Idx), 0, 2) / sqrt(M_NE);

t = -5 : 0.1 : 10;

f = figure;
f_plotLineError(t, meanSig, error, ...
    color = c_darkCyan, ...
    xlim = [-3, 7], ...
    xlabel = 'Time (s)', ...
    ylabel = 'a.u.', ...
    title = 'NE IRF');

saveas(f, fullfile(fig_savePath, 'ExtDataFig2_J1.svg'));
T = table(t', meanSig, error, ...
    VariableNames = {'Time', 'mean', 'SEM'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig2_J1.csv'));

f = figure;
f_plotMap(mean(subAvg.FigE2.NE_IRF_perf(:, :, NE_Idx), 3, 'omitnan') ...
    .* plotBM, ...
    cmp = cmpvir, ...
    clim = [0, 1]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'ExtDataFig2_J2.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);

barData = {};
barData{1} = squeeze(mean(subAvg.FigE2.NE_IRF_perf(:, :, NE_Idx) .* ...
    plotBM, [1, 2], 'omitnan'));

f = figure;
[dataMean, dataSEM] = f_plotBar(barData, ...
    colors = c_darkCyan, ...
    ylabel = 'r', ...
    ylim = [0, 1]);

saveas(f, fullfile(fig_savePath, 'ExtDataFig2_J3.svg'));
T = table({order(NE_Idx).Mouse}', barData{1}, ...
    VariableNames = {'Mouse', 'r'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig2_J3.csv'));

%% Extended Data Fig 2K

f = figure;
f_plotFC(mean(subAvg.FigE2.GRAB_conn(:, :, NE_Idx), 3), 1, ...
    cmp = cmpvir, ...
    clim = [0, 1]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'ExtDataFig2_K.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);
writetable(table(mean(subAvg.FigE2.GRAB_conn(:, :, NE_Idx), 3)), ...
    fullfile(fig_savePath, 'ExtDataFig2_K.csv'));

%% Extended Data Fig 2L

f = figure;
meanSig1 = mean(subAvg.FigE2.COH_rfp_HD_HbT, 2);
error1 = std(subAvg.FigE2.COH_rfp_HD_HbT, 0, 2) / sqrt(M);
f_plotLineError(fr, meanSig1, error1, ...
    color = c_Ca, ...
    log = 1, ...
    ylim = [0, 1], ...
    ylabel = 'Coherence');
set(gca, YScale = 'linear');

yyaxis right;
meanSig2 = mean(subAvg.FigE2.PHI_rfp_HD_HbT, 2);
error2 = std(subAvg.FigE2.PHI_rfp_HD_HbT, 0, 2) / sqrt(M);
f_plotLineError(fr, meanSig2, error2, ...
    color = [0, 0, 0], ...
    log = 1, ...
    ylim = pi * [-1, 1], ...
    xlim = [0.05, 5], ...
    xlabel = 'F (Hz)', ...
    ylabel = 'Phi (rad)');
set(gca, YScale = 'linear');

saveas(f, fullfile(fig_savePath, 'ExtDataFig2_L1.svg'));
T = table(fr, meanSig1, meanSig2, error1, error2, ...
    VariableNames = {'F (Hz)', 'mean_coherence', 'mean_phase', ...
    'SEM_coherence', 'SEM_phase'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig2_L1.csv'));


f = figure;
meanSig1 = mean(subAvg.FigE2.COH_gfp_HD_HbT(:, NE_Idx), 2);
error1 = std(subAvg.FigE2.COH_gfp_HD_HbT(:, NE_Idx), 0, 2) / sqrt(M_NE);
f_plotLineError(fr, meanSig1, error1, ...
    color = c_GRAB, ...
    log = 1, ...
    ylim = [0, 1], ...
    ylabel = 'Coherence');
set(gca, YScale = 'linear');

yyaxis right;
meanSig2 = mean(subAvg.FigE2.PHI_gfp_HD_HbT(:, NE_Idx), 2);
error2 = std(subAvg.FigE2.PHI_gfp_HD_HbT(:, NE_Idx), 0, 2) / sqrt(M_NE);
f_plotLineError(fr, meanSig2, error2, ...
    color = [0, 0, 0], ...
    log = 1, ...
    ylim = pi * [-1, 1], ...
    xlim = [0.05, 5], ...
    xlabel = 'F (Hz)', ...
    ylabel = 'Phi (rad)');
set(gca, YScale = 'linear');

saveas(f, fullfile(fig_savePath, 'ExtDataFig2_L2.svg'));
T = table(fr, meanSig1, meanSig2, error1, error2, ...
    VariableNames = {'F (Hz)', 'mean_coherence', 'mean_phase', ...
    'SEM_coherence', 'SEM_phase'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig2_L2.csv'));


f = figure;
meanSig1 = mean(subAvg.FigE2.COH_gfp_HD_HbT(:, ACh_Idx), 2);
error1 = std(subAvg.FigE2.COH_gfp_HD_HbT(:, ACh_Idx), 0, 2) / sqrt(M_ACh);
f_plotLineError(fr, meanSig1, error1, ...
    color = c_Orange, ...
    log = 1, ...
    ylim = [0, 1], ...
    ylabel = 'Coherence');
set(gca, YScale = 'linear');

yyaxis right;
meanSig2 = mean(subAvg.FigE2.PHI_gfp_HD_HbT(:, ACh_Idx), 2);
error2 = std(subAvg.FigE2.PHI_gfp_HD_HbT(:, ACh_Idx), 0, 2) / sqrt(M_ACh);
f_plotLineError(fr, meanSig2, error2, ...
    color = [0, 0, 0], ...
    log = 1, ...
    ylim = pi * [-1, 1], ...
    xlim = [0.05, 5], ...
    xlabel = 'F (Hz)', ...
    ylabel = 'Phi (rad)');
set(gca, YScale = 'linear');

saveas(f, fullfile(fig_savePath, 'ExtDataFig2_L3.svg'));
T = table(fr, meanSig1, meanSig2, error1, error2, ...
    VariableNames = {'F (Hz)', 'mean_coherence', 'mean_phase', ...
    'SEM_coherence', 'SEM_phase'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig2_L3.csv'));

%% Extended Data Fig 2M

t = 10 : -0.1 : -10;

f = figure;
meanSig1 = mean(subAvg.FigE2.XC_gfp_HD_HbT(:, NE_Idx), 2);
error1 = std(subAvg.FigE2.XC_gfp_HD_HbT(:, NE_Idx), 0, 2) / sqrt(M_NE);
f_plotLineError(t, meanSig1, error1, ...
    color = c_GRAB);

meanSig2 = mean(subAvg.FigE2.XC_gfp_HD_HbT(:, ACh_Idx), 2);
error2 = std(subAvg.FigE2.XC_gfp_HD_HbT(:, ACh_Idx), 0, 2) / sqrt(M_ACh);
f_plotLineError(t, meanSig2, error2, ...
    color = c_Orange);

meanSig3 = mean(subAvg.FigE2.XC_rfp_HD_HbT, 2);
error3 = std(subAvg.FigE2.XC_rfp_HD_HbT, 0, 2) / sqrt(M);
f_plotLineError(t, meanSig3, error3, ...
    color = c_Ca, ...
    xlim = [-5, 5], ...
    xlabel = 'Time (s)', ...
    ylabel = 'r', ...
    title = 'x vs. HbT', ...
    legend = {'NE', 'ACh', 'Ca^2^+'});

saveas(f, fullfile(fig_savePath, 'ExtDataFig2_M.svg'));
T = table(t', meanSig1, meanSig2, meanSig3, error1, error2, error3, ...
    VariableNames = {'Lag', 'mean_NE', 'mean_ACh', 'mean_Ca', ...
    'SEM_NE', 'SEM_ACh', 'SEM_Ca'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig2_M.csv'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EXTRA FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% adjust SPG length (if applicable
function spg = f_adjust_SPG(data, length, norm)
    N = numel(data);
    spg = data;
    for i = 1 : N
        if numel(spg{i}) ~= length
            spg{i} = movmean(spg{i}, 2);
            spg{i} = spg{i}(1 : 2 : end);
        end
        if norm
            spg{i} = spg{i} / mean(spg{i}) / 5;
        end
    end
end
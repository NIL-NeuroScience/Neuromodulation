%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               plotFig2
% author - Brad Rauscher (created 2024)
% 
% Plots figure panels for Figure 2. Must run MAIN_plotFigures.m first.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

% calculate subject averages
NE_order = order(NE_Idx);

tmp_data = struct;
tmp_data.s_LR_perf = f_regImages(shuffled.LR_perf, refParcellation, ...
    settings.allen_masks, 0) .* BM;
tmp_data.s_IRFx2_perf = f_regImages(shuffled.IRFx2_perf, ...
    refParcellation, settings.allen_masks, 0) .* BM;
tmp_data.g_LR_perf = f_regImages(Fig2.LR_perf, refParcellation, ...
    settings.allen_masks, 0) .* BM;
tmp_data.g_IRFx2_perf = f_regImages(Fig2.IRFx2_perf, refParcellation, ...
    settings.allen_masks, 0) .* BM;
tmp_data.LR_A = f_regImages(Fig2.LR_A, refParcellation, ...
    settings.allen_masks, 0) .* BM;
tmp_data.LR_B = f_regImages(Fig2.LR_B, refParcellation, ...
    settings.allen_masks, 0) .* BM;
tmp_data.IRFx2_A = f_regImages(Fig2.IRFx2_A, refParcellation, ...
    settings.allen_masks, 0) .* BM;
tmp_data.IRFx2_B = f_regImages(Fig2.IRFx2_B, refParcellation, ...
    settings.allen_masks, 0) .* BM;

subAvg.Fig2.s_LR_perf = NaN(500, 600, M_NE);
subAvg.Fig2.s_IRFx2_perf = NaN(500, 600, M_NE);
subAvg.Fig2.g_LR_perf = NaN(500, 600, M_NE);
subAvg.Fig2.g_IRFx2_perf = NaN(500, 600, M_NE);
subAvg.Fig2.LR_A = NaN(500, 600, M_NE);
subAvg.Fig2.LR_B = NaN(500, 600, M_NE);
subAvg.Fig2.IRFx2_A = NaN(500, 600, M_NE);
subAvg.Fig2.IRFx2_B = NaN(500, 600, M_NE);
subAvg.Fig2.tA = NaN(M_NE, 1);
subAvg.Fig2.tB = NaN(M_NE, 1);
subAvg.Fig2.IRFx2_IRF = NaN(151, 2, M_NE);

for i = 1 : M_NE
    subAvg.Fig2.s_LR_perf(:, :, i) = ...
        mean(tmp_data.s_LR_perf(:, :, NE_order(i).Runs), 3, 'omitnan');
    subAvg.Fig2.s_IRFx2_perf(:, :, i) = ...
        mean(tmp_data.s_IRFx2_perf(:, :, NE_order(i).Runs), 3, 'omitnan');
    subAvg.Fig2.g_LR_perf(:, :, i) = ...
        mean(tmp_data.g_LR_perf(:, :, NE_order(i).Runs), 3, 'omitnan');
    subAvg.Fig2.g_IRFx2_perf(:, :, i) = ...
        mean(tmp_data.g_IRFx2_perf(:, :, NE_order(i).Runs), 3, 'omitnan');
    subAvg.Fig2.LR_A(:, :, i) = ...
        mean(tmp_data.LR_A(:, :, NE_order(i).Runs), 3, 'omitnan');
    subAvg.Fig2.LR_B(:, :, i) = ...
        mean(tmp_data.LR_B(:, :, NE_order(i).Runs), 3, 'omitnan');
    subAvg.Fig2.IRFx2_A(:, :, i) = ...
        mean(tmp_data.IRFx2_A(:, :, NE_order(i).Runs), 3, 'omitnan');
    subAvg.Fig2.IRFx2_B(:, :, i) = ...
        mean(tmp_data.IRFx2_B(:, :, NE_order(i).Runs), 3, 'omitnan');
    subAvg.Fig2.LR_tA(i) = mean([Fig2.LR_tA{NE_order(i).Runs}]);
    subAvg.Fig2.LR_tB(i) = mean([Fig2.LR_tB{NE_order(i).Runs}]);
    subAvg.Fig2.IRFx2_IRF(:, :, i) = ...
        mean(cat(3, Fig2.IRFx2_IRF{NE_order(i).Runs}), 3);
end

fig_savePath = fullfile(savePath, 'Figure2');
[~, ~, ~] = mkdir(fig_savePath);

%% Fig 2A

run = 80;
runData = ...
    load('sub-Thy1-301_ses-24-08-01_run-03_irun-01_behavior+ophys.mat');

Ca = spectra.Ca{run};
HbT = f_bpf(spectra.HbT{run}, [0, 0.5], 10);
NE = GRAB_FC.GRAB_global{run};
Ca_LR = runData.Fig2.LR_Ca;
NE_LR = runData.Fig2.LR_NE;
Ca_IRFx2 = runData.Fig2.IRFx2_Ca;
NE_IRFx2 = runData.Fig2.IRFx2_NE;

Ca = Ca ./ std(Ca, 0);
HbT = HbT ./ std(HbT, 0);
NE = NE ./ std(NE, 0);
t = 0.1 : 0.1 : 600;

f = figure;
tiledlayout(3, 1);

nexttile();
hold on;
plot(t, Ca(:, 12), color=c_Ca);
plot(t, NE - 6, color=c_GRAB);
plot(t, HbT(:, 12) - 12, color=c_HbT);
plot([5, 5], [-3, 3], '-k', LineWidth = 2);
plot([5, 5], [-9, -3], '-k', LineWidth = 2);
plot([5, 5], [-15, -9], '-k', LineWidth = 2);
xlim([0, 300]);
ylim([-20, 5]);
axis off;

nexttile();
hold on;
plot(t, Ca_LR(:, 12), color = c_Ca);
plot(t, NE_LR(:, 12) - 3, color = c_GRAB);
plot(t, HbT(:, 12) - 8, color = c_HbT);
plot(t, Ca_LR(:, 12) + NE_LR(:, 12) - 8, color = [0, 0.7, 0.7]);
xlim([0, 300]);
ylim([-20, 5]);
axis off;

f_corr(HbT(:, 12), Ca_LR(:, 12) + NE_LR(:, 12), 1);

nexttile();
hold on;
plot(t, Ca_IRFx2(:, 12), color = c_Ca);
plot(t, NE_IRFx2(:, 12) - 3, color = c_GRAB);
plot(t, HbT(:, 12) - 8, color = c_HbT);
plot(t, Ca_IRFx2(:, 12) + NE_IRFx2(:, 12) - 8, color = [0, 0.7, 0.7]);
xlim([0, 300]);
ylim([-20, 5]);

f_corr(HbT(:, 12), Ca_IRFx2(:, 12) + NE_IRFx2(:, 12), 1);

idx = 1 : 3000;

headers = {'Time', 'Ca_norm', 'NE_norm', 'HbT_norm', ...
    'Ca_LR', 'NE_LR', 'HbT_LR', 'Ca_IRFx2', 'NE_IRFx2', 'HbT_IRFx2'};
T_data = [t(idx)', Ca(idx, 12), NE(idx), HbT(idx, 12), Ca_LR(idx, 12), ...
    NE_LR(idx, 12), Ca_LR(idx, 12) + NE_LR(idx, 12), Ca_IRFx2(idx, 12), ...
    NE_IRFx2(idx, 12), Ca_IRFx2(idx, 12) + NE_IRFx2(idx, 12)];

T = table();
for i = 1 : 10
    T.(headers{i}) = T_data(:, i);
end

writetable(T, fullfile(fig_savePath, 'Figure2_A.csv'));
saveas(f, fullfile(fig_savePath, 'Figure2_A.svg'));

%% Fig 2B

f = figure;
f_plotMap(mean(subAvg.Fig2.LR_A, 3, 'omitnan') .* plotBM, ...
    cmp = cmpbbr, ...
    clim = [-1, 1]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'Figure2_B1.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);

f = figure;
f_plotMap(mean(subAvg.Fig2.LR_B, 3, 'omitnan') .* plotBM, ...
    cmp = cmpbbr, ...
    clim = [-1, 1]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'Figure2_B2.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);

%% Fig 2C

barData = {};
barData{1} = subAvg.Fig2.LR_tA / 10;
barData{2} = subAvg.Fig2.LR_tB / 10;

f = figure;
[dataMean, dataSEM] = f_plotBar(barData, ...
    colors = [c_Ca; c_GRAB], ...
    legend = {'tA', 'tB'}, ...
    ylabel = 'r', ...
    title = 'Timing Coefficients');

T = table({NE_order.Mouse}', barData{1}', barData{2}', ...
    VariableNames = {'Mouse', 'tA', 'tB'});
writetable(T, fullfile(fig_savePath, 'Figure2_C.csv'));
saveas(f, fullfile(fig_savePath, 'Figure2_C.svg'));

%% Fig 2D

f = figure;
f_plotMap(mean(subAvg.Fig2.g_LR_perf, 3, 'omitnan') .* plotBM, ...
    cmp = cmpvir, ...
    clim = [0, 1]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'Figure2_D.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);

%% Fig 2E

f = figure;
f_plotMap((mean(subAvg.Fig2.g_LR_perf, 3, 'omitnan') - ...
    mean(subAvg.Fig1.inv_perf, 3, 'omitnan')) .* plotBM, ...
    cmp = cmpbbr, ...
    clim = [-1, 1]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'Figure2_E.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);

%% Fig 2F

f = figure;
f_plotMap(mean(subAvg.Fig2.IRFx2_A, 3, 'omitnan') .* plotBM, ...
    cmp = cmpbbr, ...
    clim = [-1, 1]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'Figure2_F1.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);

f = figure;
f_plotMap(mean(subAvg.Fig2.IRFx2_B, 3, 'omitnan') .* plotBM, ...
    cmp = cmpbbr, ...
    clim = [-1, 1]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'Figure2_F2.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);

%% Fig 2G

t = -5 : 0.1 : 10;

f = figure;

meanSig1 = mean(subAvg.Fig2.IRFx2_IRF(:, 1, :), 3);
error1 = std(subAvg.Fig2.IRFx2_IRF(:, 1, :), 0, 3) / sqrt(M_NE);
f_plotLineError(t, meanSig1, error1, ...
    color = c_Ca);

meanSig2 = mean(subAvg.Fig2.IRFx2_IRF(:, 2, :), 3);
error2 = std(subAvg.Fig2.IRFx2_IRF(:, 2, :), 0, 3) / sqrt(M_NE);
f_plotLineError(t, meanSig2, error2, ...
    color = c_GRAB, ...
    xlim = [-2, 7], ...
    xlabel = 'Time (s)', ...
    ylabel = 'a.u.', ...
    legend = {'IRF(t_0^A,\tau_A)', 'IRF(t_0^B,\tau_B)'}, ...
    title = 'Double IRF');

idx = 30:120;

T = table(t(idx)', meanSig1(idx), meanSig2(idx), error1(idx), ...
    error2(idx), VariableNames = {'Time', 'IRF_Ca', 'IRF_NE', 'SEM_Ca', ...
    'SEM_NE'});
writetable(T, fullfile(fig_savePath, 'Figure2_G.csv'));
saveas(f, fullfile(fig_savePath, 'Figure2_G.svg'));

%% Fig 2H

f = figure;
f_plotMap(mean(subAvg.Fig2.g_IRFx2_perf, 3, 'omitnan') .* plotBM, ...
    cmp = cmpvir, ...
    clim = [0, 1]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'Figure2_H.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);

%% Fig 2I

f = figure;
f_plotMap((mean(subAvg.Fig2.g_IRFx2_perf, 3, 'omitnan') - ...
    mean(subAvg.Fig1.inv_perf, 3, 'omitnan')) .* plotBM, ...
    cmp = cmpbbr, ...
    clim = [-1, 1]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath,'Figure2_I.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);

%% Fig 2J

f = figure;
f_plotMap(mean(subAvg.Fig2.s_LR_perf, 3, 'omitnan') .* plotBM, ...
    cmp = cmpvir, ...
    clim = [0, 1]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'Figure2_J1.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);

f = figure;
f_plotMap(mean(subAvg.Fig2.s_IRFx2_perf, 3, 'omitnan') .* plotBM, ...
    cmp = cmpvir, ...
    clim = [0, 1]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'Figure2_J2.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1 1 1]);

%% Fig 2K

barData = {};
barData{1} = squeeze(mean(subAvg.Fig1.inv_perf .* plotBM, [1, 2], ...
    'omitnan'));
barData{2} = squeeze(mean(subAvg.Fig1.SSp_perf .* plotBM, [1, 2], ...
    'omitnan'));
barData{3} = squeeze(mean(subAvg.Fig1.var_perf .* plotBM, [1, 2], ...
    'omitnan'));
barData{4} = squeeze(mean(subAvg.Fig2.g_LR_perf .* plotBM, [1, 2], ...
    'omitnan'));
barData{5} = squeeze(mean(subAvg.Fig2.g_IRFx2_perf .* plotBM, [1, 2], ...
    'omitnan'));
barData{6} = squeeze(mean(subAvg.Fig2.s_LR_perf .* plotBM, [1, 2], ...
    'omitnan'));
barData{7} = squeeze(mean(subAvg.Fig2.s_IRFx2_perf .* plotBM, [1, 2], ...
    'omitnan'));

f = figure;
[dataMean, dataSEM] = f_plotBar(barData, ...
    colors = [repmat(c_Yellow, 3, 1); repmat(c_darkCyan, 2, 1); ...
    repmat(c_pupil, 2, 1)], ...
    ylabel = 'r', ...
    title = 'Model Performance Comparison');
ylim([0, 1]);

saveas(f, fullfile(fig_savePath, 'Figure2_K.svg'));

T_data = NaN(numel(order), 7);
T_data(:, 1) = barData{1};
T_data(:, 2) = barData{2};
T_data(:, 3) = barData{3};
T_data(NE_Idx, 4) = barData{4};
T_data(NE_Idx, 5) = barData{5};
T_data(NE_Idx, 6) = barData{6};
T_data(NE_Idx, 7) = barData{7};

[h, p] = f_kstest(barData, 0.005);

labels = {'global', 'SSp', 'variant', 'Linear Regression', ...
    'Double IRF', 'Shuffled Linear Regression', 'Shuffled Double IRF'};

T = table;
T.Mouse = {order.Mouse}';
for i = 1 : numel(labels)
    T.(labels{i}) = T_data(:, i);
end
writetable(T, fullfile(fig_savePath, 'Figure2_K.csv'));

T = table;
for i = 1 : numel(labels)
    T.(labels{i}) = p(:, i);
end
writetable(T, fullfile(fig_savePath, 'Figure2_K_p.csv'));
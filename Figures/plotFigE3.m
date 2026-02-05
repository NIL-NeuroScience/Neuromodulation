%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              plotFigE3
% author - Brad Rauscher (created 2024)
% 
% Plots figure panels for Extended Data Figure 3. Must run 
% MAIN_plotFigures.m first.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

% calculate subject averages
NE_order = order(NE_Idx);

tmp_data = struct;
tmp_data.g_LR_perf = f_regImages(Fig2.LR_perf, refParcellation, settings.allen_masks, 0) .* BM;
tmp_data.g_IRFx2_perf = f_regImages(Fig2.IRFx2_perf, refParcellation, settings.allen_masks, 0) .* BM;
tmp_data.Hb.LR_perf = f_regImages(Hb_model.Hb_LR_perf, refParcellation, settings.allen_masks, 0) .* BM;
tmp_data.Hb.IRFx2_perf = f_regImages(Hb_model.Hb_IRFx2_perf, refParcellation, settings.allen_masks, 0) .* BM;
tmp_data.HbO.LR_perf = f_regImages(Hb_model.HbO_LR_perf, refParcellation, settings.allen_masks, 0) .* BM;
tmp_data.HbO.IRFx2_perf = f_regImages(Hb_model.HbO_IRFx2_perf, refParcellation, settings.allen_masks, 0) .* BM;
tmp_data.Hb.LR_A = f_regImages(Hb_model.Hb_LR_A, refParcellation, settings.allen_masks, 0) .* BM;
tmp_data.Hb.LR_B = f_regImages(Hb_model.Hb_LR_B, refParcellation, settings.allen_masks, 0) .* BM;
tmp_data.Hb.IRFx2_A = f_regImages(Hb_model.Hb_IRFx2_A, refParcellation, settings.allen_masks, 0) .* BM;
tmp_data.Hb.IRFx2_B = f_regImages(Hb_model.Hb_IRFx2_B, refParcellation, settings.allen_masks, 0) .* BM;
tmp_data.HbO.LR_A = f_regImages(Hb_model.HbO_LR_A, refParcellation, settings.allen_masks, 0) .* BM;
tmp_data.HbO.LR_B = f_regImages(Hb_model.HbO_LR_B, refParcellation, settings.allen_masks, 0) .* BM;
tmp_data.HbO.IRFx2_A = f_regImages(Hb_model.HbO_IRFx2_A, refParcellation, settings.allen_masks, 0) .* BM;
tmp_data.HbO.IRFx2_B = f_regImages(Hb_model.HbO_IRFx2_B, refParcellation, settings.allen_masks, 0) .* BM;

subAvg.FigE3.g_LR_perf = NaN(500, 600, M_NE);
subAvg.FigE3.g_IRFx2_perf = NaN(500, 600, M_NE);
subAvg.FigE3.Hb.LR_perf = NaN(500, 600, M_NE);
subAvg.FigE3.Hb.IRFx2_perf = NaN(500, 600, M_NE);
subAvg.FigE3.Hb.LR_A = NaN(500, 600, M_NE);
subAvg.FigE3.Hb.LR_B = NaN(500, 600, M_NE);
subAvg.FigE3.Hb.IRFx2_A = NaN(500, 600, M_NE);
subAvg.FigE3.Hb.IRFx2_B = NaN(500, 600, M_NE);
subAvg.FigE3.Hb.tA = NaN(M_NE, 1);
subAvg.FigE3.Hb.tB = NaN(M_NE, 1);
subAvg.FigE3.Hb.IRFx2_IRF = NaN(151, 2, M_NE);
subAvg.FigE3.HbO.LR_perf = NaN(500, 600, M_NE);
subAvg.FigE3.HbO.IRFx2_perf = NaN(500, 600, M_NE);
subAvg.FigE3.HbO.LR_A = NaN(500, 600, M_NE);
subAvg.FigE3.HbO.LR_B = NaN(500, 600, M_NE);
subAvg.FigE3.HbO.IRFx2_A = NaN(500, 600, M_NE);
subAvg.FigE3.HbO.IRFx2_B = NaN(500, 600, M_NE);
subAvg.FigE3.HbO.tA = NaN(M_NE, 1);
subAvg.FigE3.HbO.tB = NaN(M_NE, 1);
subAvg.FigE3.HbO.IRFx2_IRF = NaN(151, 2, M_NE);

for i = 1 : M_NE
    subAvg.FigE3.g_LR_perf(:, :, i) = mean(tmp_data.g_LR_perf(:, :, NE_order(i).Runs), 3, 'omitnan');
    subAvg.FigE3.g_IRFx2_perf(:, :, i) = mean(tmp_data.g_IRFx2_perf(:, :, NE_order(i).Runs), 3, 'omitnan');
    subAvg.FigE3.Hb.LR_perf(:, :, i) = mean(tmp_data.Hb.LR_perf(:, :, NE_order(i).Runs), 3, 'omitnan');
    subAvg.FigE3.Hb.IRFx2_perf(:, :, i) = mean(tmp_data.Hb.IRFx2_perf(:, :, NE_order(i).Runs), 3, 'omitnan');
    subAvg.FigE3.HbO.LR_perf(:, :, i) = mean(tmp_data.HbO.LR_perf(:, :, NE_order(i).Runs), 3, 'omitnan');
    subAvg.FigE3.HbO.IRFx2_perf(:, :, i) = mean(tmp_data.HbO.IRFx2_perf(:, :, NE_order(i).Runs), 3, 'omitnan');
    subAvg.FigE3.Hb.LR_A(:, :, i) = mean(tmp_data.Hb.LR_A(:, :, NE_order(i).Runs), 3, 'omitnan');
    subAvg.FigE3.Hb.LR_B(:, :, i) = mean(tmp_data.Hb.LR_B(:, :, NE_order(i).Runs), 3, 'omitnan');
    subAvg.FigE3.Hb.IRFx2_A(:, :, i) = mean(tmp_data.Hb.IRFx2_A(:, :, NE_order(i).Runs), 3, 'omitnan');
    subAvg.FigE3.Hb.IRFx2_B(:, :, i) = mean(tmp_data.Hb.IRFx2_B(:, :, NE_order(i).Runs), 3, 'omitnan');
    subAvg.FigE3.Hb.LR_tA(i) = mean([Hb_model.Hb_LR_tA{NE_order(i).Runs}]);
    subAvg.FigE3.Hb.LR_tB(i) = mean([Hb_model.Hb_LR_tB{NE_order(i).Runs}]);
    subAvg.FigE3.Hb.IRFx2_IRF(:, :, i) = mean(cat(3,Hb_model.Hb_IRFx2_IRF{NE_order(i).Runs}), 3);
    subAvg.FigE3.HbO.LR_A(:, :, i) = mean(tmp_data.HbO.LR_A(:, :, NE_order(i).Runs), 3, 'omitnan');
    subAvg.FigE3.HbO.LR_B(:, :, i) = mean(tmp_data.HbO.LR_B(:, :, NE_order(i).Runs), 3, 'omitnan');
    subAvg.FigE3.HbO.IRFx2_A(:, :, i) = mean(tmp_data.HbO.IRFx2_A(:, :, NE_order(i).Runs), 3, 'omitnan');
    subAvg.FigE3.HbO.IRFx2_B(:, :, i) = mean(tmp_data.HbO.IRFx2_B(:, :, NE_order(i).Runs), 3, 'omitnan');
    subAvg.FigE3.HbO.LR_tA(i) = mean([Hb_model.HbO_LR_tA{NE_order(i).Runs}]);
    subAvg.FigE3.HbO.LR_tB(i) = mean([Hb_model.HbO_LR_tB{NE_order(i).Runs}]);
    subAvg.FigE3.HbO.IRFx2_IRF(:, :, i) = mean(cat(3,Hb_model.HbO_IRFx2_IRF{NE_order(i).Runs}), 3);
end

fig_savePath = fullfile(savePath, 'ExtDataFig3');
[~, ~, ~] = mkdir(fig_savePath);

%% Extended Data Fig 3A

f = figure;
f_plotMap(mean(subAvg.FigE3.HbO.LR_A, 3, 'omitnan') .* plotBM, ...
    cmp = cmpbbr, ...
    clim = [-1, 1]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'ExtDataFig3_A1.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);

f = figure;
f_plotMap(mean(subAvg.FigE3.HbO.LR_B, 3, 'omitnan') .* plotBM, ...
    cmp = cmpbbr, ...
    clim = [-1, 1]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'ExtDataFig3_A2.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);

%% Extended Data Fig 3B

barData = {};
barData{1} = subAvg.FigE3.HbO.LR_tA / 10;
barData{2} = subAvg.FigE3.HbO.LR_tB / 10;

f = figure;
[meanSig, SEM] = f_plotBar(barData, ...
    colors = [c_Ca; c_GRAB], ...
    legend = {'tA', 'tB'}, ...
    ylabel = 'r', ...
    title = 'Timing Coefficients');
saveas(f, fullfile(fig_savePath, 'ExtDataFig3_B.svg'));

T = table({order(NE_Idx).Mouse}', barData{1}', barData{2}', ...
    VariableNames = {'Mouse', 'tA', 'tB'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig3_B.csv'));

%% Extended Data Fig 3C

f = figure;
f_plotMap(mean(subAvg.FigE3.HbO.LR_perf, 3, 'omitnan') .* plotBM, ...
    cmp = cmpvir, ...
    clim = [0, 1]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'ExtDataFig3_C.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);

%% Extended Data Fig 3D

f = figure;
f_plotMap((mean(subAvg.FigE3.HbO.LR_perf, 3, 'omitnan') - ...
    mean(subAvg.Fig2.g_LR_perf, 3, 'omitnan')) .* plotBM, ...
    cmp = cmpbbr, ...
    clim = [-1, 1]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'ExtDataFig3_D.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);

%% Extended Data Fig 3E

f = figure;
f_plotMap(mean(subAvg.FigE3.Hb.LR_A, 3, 'omitnan') .* plotBM, ...
    cmp = cmpbbr, ...
    clim = [-1, 1]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'ExtDataFig3_E1.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);

f = figure;
f_plotMap(mean(subAvg.FigE3.Hb.LR_B, 3, 'omitnan') .* plotBM, ...
    cmp = cmpbbr, ...
    clim = [-1, 1]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'ExtDataFig3_E2.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);

%% Extended Data Fig 3F

barData = {};
barData{1} = subAvg.FigE3.Hb.LR_tA / 10;
barData{2} = subAvg.FigE3.Hb.LR_tB / 10;

f = figure;
[meanSig, SEM] = f_plotBar(barData, ...
    colors = [c_Ca; c_GRAB], ...
    legend = {'tA', 'tB'}, ...
    ylabel = 'r', ...
    title = 'Timing Coefficients', ...
    ylim = [0, 1]);

saveas(f, fullfile(fig_savePath, 'ExtDataFig3_F.svg'));

T = table({order(NE_Idx).Mouse}', barData{1}', barData{2}', ...
    VariableNames = {'Mouse', 'tA', 'tB'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig3_F.csv'));

%% Extended Data Fig 3G

f = figure;
f_plotMap(mean(subAvg.FigE3.Hb.LR_perf, 3, 'omitnan') .* plotBM, ...
    cmp = cmpvir, ...
    clim = [0, 1]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'ExtDataFig3_G.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);

%% Extended Data Fig 3H

f = figure;
f_plotMap((mean(subAvg.FigE3.Hb.LR_perf, 3, 'omitnan') - ...
    mean(subAvg.Fig2.g_LR_perf, 3, 'omitnan')) .* plotBM, ...
    cmp = cmpbbr, ...
    clim = [-1, 1]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'ExtDataFig3_H.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);

%% Extended Data Fig 3I

f = figure;
f_plotMap(mean(subAvg.FigE3.HbO.IRFx2_A, 3, 'omitnan') .* plotBM, ...
    cmp = cmpbbr, ...
    clim = [-1, 1]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'ExtDataFig3_I1.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);

f = figure;
f_plotMap(mean(subAvg.FigE3.HbO.IRFx2_B, 3, 'omitnan') .* plotBM, ...
    cmp = cmpbbr, ...
    clim = [-1, 1]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'ExtDataFig3_I2.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);

%% Extended Data Fig 3J

t = -5 : 0.1 : 10;

f = figure;
meanSig1 = mean(subAvg.FigE3.HbO.IRFx2_IRF(:, 1, :), 3);
error1 = std(subAvg.FigE3.HbO.IRFx2_IRF(:, 1, :), 0, 3) / sqrt(M_NE);
f_plotLineError(t, meanSig1, error1, ...
    color = c_Ca);

meanSig2 = mean(subAvg.FigE3.HbO.IRFx2_IRF(:, 2, :), 3);
error2 = std(subAvg.FigE3.HbO.IRFx2_IRF(:, 2, :), 0, 3) / sqrt(M_NE);
f_plotLineError(t, meanSig2, error2, ...
    color = c_GRAB, ...
    xlim = [-2, 7], ...
    xlabel = 'Time (s)', ...
    ylabel = 'a.u.', ...
    legend = {'IRF(t_0^A,\tau_A)', 'IRF(t_0^B,\tau_B)'}, ...
    title = 'Double IRF');

idx = 30 : 120;

T = table(t(idx)', meanSig1(idx), meanSig2(idx), error1(idx), ...
    error2(idx), VariableNames = {'Time', 'IRF_Ca', 'IRF_NE', 'SEM_Ca', ...
    'SEM_NE'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig3_J.csv'));
saveas(f, fullfile(fig_savePath, 'ExtDataFig3_J.svg'));

%% Extended Data Fig 3K

f = figure;
f_plotMap(mean(subAvg.FigE3.HbO.IRFx2_perf, 3, 'omitnan') .* plotBM, ...
    cmp = cmpvir, ...
    clim = [0, 1]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'ExtDataFig3_K.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);

%% Extended Data Fig 3L

f = figure;
f_plotMap((mean(subAvg.FigE3.HbO.IRFx2_perf, 3, 'omitnan') - ...
    mean(subAvg.Fig2.g_IRFx2_perf, 3, 'omitnan')) .* plotBM, ...
    cmp = cmpbbr, ...
    clim = [-1, 1]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'ExtDataFig3_L.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);

%% Extended Data Fig 3M

f = figure;
f_plotMap(mean(subAvg.FigE3.Hb.IRFx2_A, 3, 'omitnan') .* plotBM, ...
    cmp = cmpbbr, ...
    clim = [-1, 1]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'ExtDataFig3_M1.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);

f = figure;
f_plotMap(mean(subAvg.FigE3.Hb.IRFx2_B, 3, 'omitnan') .* plotBM, ...
    cmp = cmpbbr, ...
    clim = [-1, 1]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'ExtDataFig3_M2.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);

%% Extended Data Fig 3N

t = -5 : 0.1 : 10;

f = figure;
meanSig1 = mean(subAvg.FigE3.Hb.IRFx2_IRF(:, 1, :), 3);
error1 = std(subAvg.FigE3.Hb.IRFx2_IRF(:, 1, :), 0, 3) / sqrt(M_NE);
f_plotLineError(t, meanSig1, error1, ...
    color = c_Ca);

meanSig2 = mean(subAvg.FigE3.Hb.IRFx2_IRF(:, 2, :), 3);
error2 = std(subAvg.FigE3.Hb.IRFx2_IRF(:, 2, :), 0, 3) / sqrt(M_NE);
f_plotLineError(t, meanSig2, error2, ...
    color = c_GRAB, ...
    xlim = [-2, 7], ...
    xlabel = 'Time (s)', ...
    ylabel = 'a.u.', ...
    legend = {'IRF(t_0^A,\tau_A)','IRF(t_0^B,\tau_B)'}, ...
    title = 'Double IRF');

idx = 30:120;

T = table(t(idx)', meanSig1(idx), meanSig2(idx), error1(idx), ...
    error2(idx), VariableNames = {'Time', 'IRF_Ca', 'IRF_NE', 'SEM_Ca', ...
    'SEM_NE'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig3_N.csv'));
saveas(f, fullfile(fig_savePath, 'ExtDataFig3_N.svg'));

%% Extended Data Fig 3O

f = figure;
f_plotMap(mean(subAvg.FigE3.Hb.IRFx2_perf, 3, 'omitnan') .* plotBM, ...
    cmp = cmpvir, ...
    clim = [0, 1]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'ExtDataFig3_O.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);

%% Extended Data Fig 3P

f = figure;
f_plotMap((mean(subAvg.FigE3.Hb.IRFx2_perf, 3, 'omitnan') - ...
    mean(subAvg.Fig2.g_IRFx2_perf, 3, 'omitnan')) .* plotBM, ...
    cmp = cmpbbr, ...
    clim = [-1, 1]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'ExtDataFig3_P.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);

%% Extended Data Fig 3Q

barData = {};
barData{1} = squeeze(mean(subAvg.Fig1.inv_perf .* plotBM, [1, 2], 'omitnan'));
barData{2} = squeeze(mean(subAvg.Fig1.SSp_perf .* plotBM, [1, 2], 'omitnan'));
barData{3} = squeeze(mean(subAvg.Fig1.var_perf .* plotBM, [1, 2], 'omitnan'));
barData{4} = squeeze(mean(subAvg.FigE3.g_LR_perf .* plotBM, [1, 2], 'omitnan'));
barData{5} = squeeze(mean(subAvg.FigE3.g_IRFx2_perf .* plotBM, [1, 2], 'omitnan'));
barData{6} = squeeze(mean(subAvg.FigE3.HbO.LR_perf .* plotBM, [1, 2], 'omitnan'));
barData{7} = squeeze(mean(subAvg.FigE3.HbO.IRFx2_perf .* plotBM, [1, 2], 'omitnan'));
barData{8} = squeeze(mean(subAvg.FigE3.Hb.LR_perf .* plotBM, [1, 2], 'omitnan'));
barData{9} = squeeze(mean(subAvg.FigE3.Hb.IRFx2_perf .* plotBM, [1, 2], 'omitnan'));

f = figure;
[dataMean, dataSEM] = f_plotBar(barData, ...
    colors = [repmat(c_Yellow, 3, 1); repmat(c_darkCyan, 2, 1); ...
    repmat(c_Ca, 2, 1); repmat(c_pupil, 2, 1)], ...
    ylabel = 'r', ...
    title = 'Model Performance Comparison', ...
    ylim = [0, 1]);

[h, p] = f_kstest(barData, 0.01);

T_data = NaN(numel(order), 9);
T_data(:, 1) = barData{1};
T_data(:, 2) = barData{2};
T_data(:, 3) = barData{3};
T_data(NE_Idx, 4) = barData{4};
T_data(NE_Idx, 5) = barData{5};
T_data(NE_Idx, 6) = barData{6};
T_data(NE_Idx, 7) = barData{7};
T_data(NE_Idx, 8) = barData{8};
T_data(NE_Idx, 9) = barData{9};

saveas(f, fullfile(fig_savePath, 'ExtDataFig3_Q.svg'));

labels = {'global', 'SSp', 'variant', 'Linear Regression', 'Double IRF', ...
    'HbO Linear Regression', 'HbO Double IRF', 'HbR Linear Regression', ...
    'HbR Double IRF'};

T = table;
T.Mouse = {order.Mouse}';
for i = 1 : numel(labels)
    T.(labels{i}) = T_data(:, i);
end
writetable(T, fullfile(fig_savePath, 'ExtDataFig3_Q.csv'));

T = table;
for i = 1 : numel(labels)
    T.(labels{i}) = p(:, i);
end
writetable(T, fullfile(fig_savePath, 'ExtDataFig3_Q_p.csv'));
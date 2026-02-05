%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              plotFigE8
% author - Brad Rauscher (created 2024)
% 
% Plots figure panels for Extended Data Figure 8. Must run 
% MAIN_plotFigures.m first.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all

% calculate subject averages
NE_order = order(NE_Idx);

subAvg.FigE8.low_NE_Ca = NaN(300, M_NE);
subAvg.FigE8.high_NE_Ca = NaN(300, M_NE);
subAvg.FigE8.low_NE_HbT = NaN(300, M_NE);
subAvg.FigE8.high_NE_HbT = NaN(300, M_NE);
subAvg.FigE8.lowF_lowNE_Ca = NaN(12, M_NE);
subAvg.FigE8.lowF_highNE_Ca = NaN(12, M_NE);
subAvg.FigE8.medF_lowNE_Ca = NaN(12, M_NE);
subAvg.FigE8.medF_highNE_Ca = NaN(12, M_NE);
subAvg.FigE8.highF_lowNE_Ca = NaN(12, M_NE);
subAvg.FigE8.highF_highNE_Ca = NaN(12, M_NE);
subAvg.FigE8.lowF_lowNE_HbT = NaN(12, M_NE);
subAvg.FigE8.lowF_highNE_HbT = NaN(12, M_NE);
subAvg.FigE8.medF_lowNE_HbT = NaN(12, M_NE);
subAvg.FigE8.medF_highNE_HbT = NaN(12, M_NE);

for i = 1 : M_NE
    subAvg.FigE8.low_NE_Ca(:, i) = mean(cat(2, ...
        spectra.low_NE_Ca{NE_order(i).Runs}), 2);
    subAvg.FigE8.high_NE_Ca(:, i) = mean(cat(2, ...
        spectra.high_NE_Ca{NE_order(i).Runs}), 2);
    subAvg.FigE8.low_NE_HbT(:, i) = mean(cat(2, ...
        spectra.low_NE_HbT{NE_order(i).Runs}), 2);
    subAvg.FigE8.high_NE_HbT(:, i) = mean(cat(2, ...
        spectra.high_NE_HbT{NE_order(i).Runs}), 2);
    subAvg.FigE8.lowF_lowNE_Ca(:, i) = mean(cat(2, ...
        spectra.lowF_lowNE_Ca{NE_order(i).Runs}), 2);
    subAvg.FigE8.lowF_highNE_Ca(:, i) = mean(cat(2, ...
        spectra.lowF_highNE_Ca{NE_order(i).Runs}), 2);
    subAvg.FigE8.medF_lowNE_Ca(:, i) = mean(cat(2, ...
        spectra.medF_lowNE_Ca{NE_order(i).Runs}), 2);
    subAvg.FigE8.medF_highNE_Ca(:, i) = mean(cat(2, ...
        spectra.medF_highNE_Ca{NE_order(i).Runs}), 2);
    subAvg.FigE8.highF_lowNE_Ca(:, i) = mean(cat(2, ...
        spectra.highF_lowNE_Ca{NE_order(i).Runs}), 2);
    subAvg.FigE8.highF_highNE_Ca(:, i) = mean(cat(2, ...
        spectra.highF_highNE_Ca{NE_order(i).Runs}), 2);
    subAvg.FigE8.lowF_lowNE_HbT(:, i) = mean(cat(2, ...
        spectra.lowF_lowNE_HbT{NE_order(i).Runs}), 2);
    subAvg.FigE8.lowF_highNE_HbT(:, i) = mean(cat(2, ...
        spectra.lowF_highNE_HbT{NE_order(i).Runs}), 2);
    subAvg.FigE8.medF_lowNE_HbT(:, i) = mean(cat(2, ...
        spectra.medF_lowNE_HbT{NE_order(i).Runs}), 2);
    subAvg.FigE8.medF_highNE_HbT(:, i) = mean(cat(2, ...
        spectra.medF_highNE_HbT{NE_order(i).Runs}), 2);
end

fig_savePath = fullfile(savePath, 'ExtDataFig8');
[~, ~, ~] = mkdir(fig_savePath);

% extract NE peaks

fr = spectra.fr;

tBef = 50;
tAft = 50;

[peaks, trials] = f_findPeaks(spectra.NE, [0, 0.02], 10, 0.005, ...
    [tBef, tAft]);

%% Extended Data Fig 8A

tmp_trials = cell(numel(spectra.NE), 1);
tmp_subAvg = NaN(tBef * 10 + tAft * 10 + 1, M_NE);

for i = 1 : numel(spectra.NE)
    for idx = 1 : numel(peaks{i})
        tmp_trials{i}(:, idx) = spectra.NE{i}(peaks{i}(idx) - tBef * 10 ...
            : peaks{i}(idx) + tAft * 10);
    end
end
for i = 1 : M_NE
    tmp_subAvg(:, i) = mean(cat(2, tmp_trials{NE_order(i).Runs}), 2);
end

t = -50 : 0.1 : 50;

f = figure(Position = [100, 100, 1000, 200]);
meanSig = mean(tmp_subAvg, 2);
error = std(tmp_subAvg, 0, 2) / sqrt(M_NE);
f_plotLineError(t, meanSig, error, ...
    color = c_GRAB, ...
    ylim = 1.5 * [-1, 1]);

clear tmp_trials tmp_subAvg;

saveas(f, fullfile(fig_savePath, 'ExtDataFig8_A.svg'));

T = table(t', meanSig, error, ...
    VariableNames = {'Time', 'mean_NE', 'SEM_NE'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig8_A.csv'));

%% Extended Data Fig 8B

tmp_subAvg = NaN(tBef * 10 + tAft * 10 + 1, 300, M_NE);

for m = 1 : M_NE
    trials_mouse = NE_order(m).Runs;
    tmp_trials = NaN(tBef * 10 + tAft * 10 + 1, 300, ...
        sum(trials(trials_mouse)));
    running_idx = 1;
    
    for i = 1 : numel(trials_mouse)
        trial_idx = trials_mouse(i);
        for p = 1 : numel(peaks{trial_idx})
            tmp_trials(:, :, running_idx) = ...
                spectra.SPG_Ca{trial_idx}(peaks{trial_idx}(p) - ...
                tBef * 10 : peaks{trial_idx}(p) + tAft * 10, :);
            running_idx = running_idx + 1;
        end
    end    

    tmp_subAvg(:, :, m) = mean(tmp_trials, 3);
end

f = figure(Position = [100, 100, 1000, 300]);
f_plotMap(flipud(mean(tmp_subAvg, 3)'), ...
    cmp = cmpinf, ...
    clim = [0, 1.25]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'ExtDataFig8_B.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);

clear tmp_trials tmp_subAvg;

%% Extended Data Fig 8C

tmp_subAvg = NaN(tBef * 10 + tAft * 10 + 1, 300, M_NE);

for m = 1 : M_NE
    trials_mouse = NE_order(m).Runs;
    tmp_trials = NaN(tBef * 10 + tAft * 10 + 1, 300, ...
        sum(trials(trials_mouse)));
    running_idx = 1;
    
    for i = 1 : numel(trials_mouse)
        trial_idx = trials_mouse(i);
        for p = 1 : numel(peaks{trial_idx})
            tmp_trials(:, :, running_idx) = ...
                spectra.SPG_HbT{trial_idx}(peaks{trial_idx}(p) - ...
                tBef * 10 : peaks{trial_idx}(p) + tAft * 10, :, :);
            running_idx = running_idx + 1;
        end
    end    

    tmp_subAvg(:, :, m) = mean(tmp_trials, 3);
end

f = figure(Position = [100, 100, 1000, 300]);
f_plotMap(flipud(mean(tmp_subAvg, 3)'), ...
    cmp = cmpinf, ...
    clim = [0, 0.2]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'ExtDataFig8_C.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);

clear tmp_trials tmp_subAvg;

%% Extended Data Fig 8C

tmp_trials = cell(numel(spectra.NE), 1);
tmp_subAvg = NaN(tBef * 10 + tAft * 10 + 1, 300, M_NE);

for i = 1 : numel(spectra.NE)
    for idx = 1 : numel(peaks{i})
        tmp_trials{i}(:, :, idx) = spectra.SPG_HbT{i}(peaks{i}(idx) ...
            - tBef * 10 : peaks{i}(idx) + tAft * 10, :);
    end
end
for i = 1 : M_NE
    tmp_subAvg(:,:,i) = mean(cat(3, tmp_trials{NE_order(i).Runs}), 3);
end

f = figure(Position = [100, 100, 1000, 300]);
f_plotMap(flipud(mean(tmp_subAvg, 3)'), ...
    cmp = cmpinf, ...
    clim = [0, 0.2]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'ExtDataFig8_C.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);

clear tmp_trials tmp_subAvg;

%% Extended Data Fig 8D

f = figure;
meanSig1 = mean(subAvg.FigE8.low_NE_Ca, 2);
error1 = std(subAvg.FigE8.low_NE_Ca, 0, 2) / sqrt(M_NE);
f_plotLineError(fr, fr .* meanSig1, fr .* error1, ...
    color = c_darkCyan, ...
    log = 1);

meanSig2 = mean(subAvg.FigE8.high_NE_Ca, 2);
error2 = std(subAvg.FigE8.high_NE_Ca, 0, 2) / sqrt(M_NE);
f_plotLineError(fr, fr .* meanSig2, fr .* error2, ...
    color = c_Orange, ...
    log = 1, ...
    xlim = [0.01, 5], ...
    xlabel = 'F (Hz', ...
    ylabel = 'Power (normalized by F)', ...
    legend = {'low NE', 'high NE'});
set(gca, YScale = 'linear');

saveas(f, fullfile(fig_savePath, 'ExtDataFig8_D1.svg'));

T = table(fr, fr .* meanSig1, fr .* meanSig2, fr .* error1, ...
    fr .* error2, VariableNames = {'F (Hz)', 'mean_low_NE', ...
    'mean_high_NE', 'SEM_low_NE', 'SEM_high_NE'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig8_D1.csv'));

f = figure;
meanSig1 = mean(subAvg.FigE8.low_NE_HbT, 2);
error1 = std(subAvg.FigE8.low_NE_HbT, 0, 2) / sqrt(M_NE);
f_plotLineError(fr, fr .* meanSig1, fr .* error1, ...
    color = c_darkCyan, ...
    log = 1);

meanSig2 = mean(subAvg.FigE8.high_NE_HbT, 2);
error2 = std(subAvg.FigE8.high_NE_HbT, 0, 2) / sqrt(M_NE);
f_plotLineError(fr, fr .* meanSig2, fr .* error2, ...
    color = c_Orange, ...
    log = 1, ...
    xlim = [0.01, 5], ...
    xlabel = 'F (Hz)', ...
    ylabel = 'Power (normalized by F)', ...
    legend = {'low NE', 'high NE'});
set(gca, YScale  = 'linear');

saveas(f, fullfile(fig_savePath, 'ExtDataFig8_D2.svg'));

T = table(fr, fr .* meanSig1, fr .* meanSig2, fr .* error1, ...
    fr .* error2, VariableNames = {'F (Hz)', 'mean_low_NE', ...
    'mean_high_NE', 'SEM_low_NE', 'SEM_high_NE'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig5_D2.csv'));

%% Extended Data Fig 8E

f = figure;
f_plotAllenMap(mean(subAvg.FigE8.lowF_lowNE_Ca, 2), ...
    cmp = cmpinf, ...
    mask = plotBM, ...
    clim=[0, 1]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'ExtDataFig8_E1.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);

f = figure;
f_plotAllenMap(mean(subAvg.FigE8.lowF_highNE_Ca, 2), ...
    cmp = cmpinf, ...
    mask = plotBM, ...
    clim = [0, 1]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'ExtDataFig8_E2.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);

f = figure;
f_plotAllenMap(mean(subAvg.FigE8.lowF_highNE_Ca, 2) - ...
    mean(subAvg.FigE8.lowF_lowNE_Ca, 2), ...
    cmp = cmpbbr, ...
    mask = plotBM, ...
    clim = 0.3 * [-1, 1]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'ExtDataFig8_E3.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);


f = figure;
f_plotAllenMap(mean(subAvg.FigE8.medF_lowNE_Ca, 2), ...
    cmp = cmpinf, ...
    mask = plotBM, ...
    clim = [0, 1]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'ExtDataFig8_E4.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);

f = figure;
f_plotAllenMap(mean(subAvg.FigE8.medF_highNE_Ca, 2), ...
    cmp = cmpinf, ...
    mask = plotBM, ...
    clim = [0, 1]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'ExtDataFig8_E5.jpg'), ...
    Resolution = 300, BackgroundColor = [1, 1, 1]);

f = figure;
f_plotAllenMap(mean(subAvg.FigE8.medF_highNE_Ca, 2) - ...
    mean(subAvg.FigE8.medF_lowNE_Ca, 2), ...
    cmp = cmpbbr, ...
    mask = plotBM, ...
    clim = 0.3 * [-1, 1]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'ExtDataFig8_E6.jpg'), ...
    Resolution = 300, BackgroundColor = [1, 1, 1]);


f = figure;
f_plotAllenMap(mean(subAvg.FigE8.highF_lowNE_Ca, 2), ...
    cmp = cmpinf, ...
    mask = plotBM, ...
    clim = [0, 0.4]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'ExtDataFig8_E7.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);

f = figure;
f_plotAllenMap(mean(subAvg.FigE8.highF_highNE_Ca, 2), ...
    cmp = cmpinf, ...
    mask = plotBM, ...
    clim = [0, 0.4]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'ExtDataFig8_E8.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);

f = figure;
f_plotAllenMap(mean(subAvg.FigE8.highF_highNE_Ca, 2) - ...
    mean(subAvg.FigE8.highF_lowNE_Ca, 2), ...
    cmp = cmpbbr, ...
    mask = plotBM, ...
    clim = 0.1 * [-1, 1]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'ExtDataFig8_E9.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);

%% Extended Data Fig 8F

f = figure;
f_plotAllenMap(mean(subAvg.FigE8.lowF_lowNE_HbT, 2), ...
    cmp = cmpinf, ...
    mask = plotBM, ...
    clim = [0, 0.15]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'ExtDataFig8_F1.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);

f = figure;
f_plotAllenMap(mean(subAvg.FigE8.lowF_highNE_HbT, 2), ...
    cmp = cmpinf, ...
    mask = plotBM, ...
    clim = [0, 0.15]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'ExtDataFig8_F2.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);

f = figure;
f_plotAllenMap(mean(subAvg.FigE8.lowF_highNE_HbT, 2) - ...
    mean(subAvg.FigE8.lowF_lowNE_HbT, 2), ...
    cmp = cmpbbr, ...
    mask = plotBM, ...
    clim = 0.02 * [-1, 1]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'ExtDataFig8_F3.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);


f = figure;
f_plotAllenMap(mean(subAvg.FigE8.medF_lowNE_HbT, 2), ...
    cmp = cmpinf, ...
    mask = plotBM, ...
    clim = [0, 0.15]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'ExtDataFig8_F4.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);

f = figure;
f_plotAllenMap(mean(subAvg.FigE8.medF_highNE_HbT, 2), ...
    cmp = cmpinf, ...
    mask = plotBM, ...
    clim = [0, 0.15]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'ExtDataFig8_F5.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);

f = figure;
f_plotAllenMap(mean(subAvg.FigE8.medF_highNE_HbT, 2) - ...
    mean(subAvg.FigE8.medF_lowNE_HbT, 2), ...
    cmp = cmpbbr, ...
    mask = plotBM, ...
    clim = 0.03 * [-1, 1]);
colorbar off;

exportgraphics(f, fullfile(fig_savePath, 'ExtDataFig8_F6.jpg'), ...
    Resolution = 300, ...
    BackgroundColor = [1, 1, 1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EXTRA FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find NE peaks
function [peaks, trials] = f_findPeaks(NE, fwin, fr, th, t)

N = numel(NE);
peaks = cell(N, 1);
trials = zeros(N, 1);

for i = 1 : N
    
    if isempty(NE{i})
        continue
    end

    sig = f_bpf(NE{i}, fwin, fr);
    diffNE = diff(sig);
    transitions = diffNE > th;

    if sum(transitions)
        tIdx = find(diff(transitions) == 1);
        if isempty(tIdx)
            continue
        end
        tDown = find(diff(transitions) == -1);
        tDown(tDown < tIdx(1)) = [];
        if numel(tIdx) > numel(tDown)
            tIdx(numel(tDown) + 1 : end) = [];
        end

        transitions = round(mean([tIdx, tDown], 2));
        transitions(transitions < t(1) * 10 + 1) = [];
        transitions(transitions > numel(sig) - t(2) * 10) = [];
        peaks{i} = transitions;
        trials(i) = numel(transitions);
    end
end

end


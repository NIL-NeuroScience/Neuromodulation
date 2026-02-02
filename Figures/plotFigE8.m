close all

% calculate subject averages

NE_order = order(NE_Idx);

subAvg.FigE8.low_NE_Ca = NaN(300,numel(NE_order));
subAvg.FigE8.high_NE_Ca = NaN(300,numel(NE_order));
subAvg.FigE8.low_NE_HbT = NaN(300,numel(NE_order));
subAvg.FigE8.high_NE_HbT = NaN(300,numel(NE_order));
subAvg.FigE8.lowF_lowNE_Ca = NaN(12,numel(NE_order));
subAvg.FigE8.lowF_highNE_Ca = NaN(12,numel(NE_order));
subAvg.FigE8.medF_lowNE_Ca = NaN(12,numel(NE_order));
subAvg.FigE8.medF_highNE_Ca = NaN(12,numel(NE_order));
subAvg.FigE8.highF_lowNE_Ca = NaN(12,numel(NE_order));
subAvg.FigE8.highF_highNE_Ca = NaN(12,numel(NE_order));
subAvg.FigE8.lowF_lowNE_HbT = NaN(12,numel(NE_order));
subAvg.FigE8.lowF_highNE_HbT = NaN(12,numel(NE_order));
subAvg.FigE8.medF_lowNE_HbT = NaN(12,numel(NE_order));
subAvg.FigE8.medF_highNE_HbT = NaN(12,numel(NE_order));

for i = 1:numel(NE_order)
    subAvg.FigE8.low_NE_Ca(:,i) = mean(cat(2,spectra.low_NE_Ca{NE_order(i).Runs}),2);
    subAvg.FigE8.high_NE_Ca(:,i) = mean(cat(2,spectra.high_NE_Ca{NE_order(i).Runs}),2);
    subAvg.FigE8.low_NE_HbT(:,i) = mean(cat(2,spectra.low_NE_HbT{NE_order(i).Runs}),2);
    subAvg.FigE8.high_NE_HbT(:,i) = mean(cat(2,spectra.high_NE_HbT{NE_order(i).Runs}),2);
    subAvg.FigE8.lowF_lowNE_Ca(:,i) = mean(cat(2,spectra.lowF_lowNE_Ca{NE_order(i).Runs}),2);
    subAvg.FigE8.lowF_highNE_Ca(:,i) = mean(cat(2,spectra.lowF_highNE_Ca{NE_order(i).Runs}),2);
    subAvg.FigE8.medF_lowNE_Ca(:,i) = mean(cat(2,spectra.medF_lowNE_Ca{NE_order(i).Runs}),2);
    subAvg.FigE8.medF_highNE_Ca(:,i) = mean(cat(2,spectra.medF_highNE_Ca{NE_order(i).Runs}),2);
    subAvg.FigE8.highF_lowNE_Ca(:,i) = mean(cat(2,spectra.highF_lowNE_Ca{NE_order(i).Runs}),2);
    subAvg.FigE8.highF_highNE_Ca(:,i) = mean(cat(2,spectra.highF_highNE_Ca{NE_order(i).Runs}),2);
    subAvg.FigE8.lowF_lowNE_HbT(:,i) = mean(cat(2,spectra.lowF_lowNE_HbT{NE_order(i).Runs}),2);
    subAvg.FigE8.lowF_highNE_HbT(:,i) = mean(cat(2,spectra.lowF_highNE_HbT{NE_order(i).Runs}),2);
    subAvg.FigE8.medF_lowNE_HbT(:,i) = mean(cat(2,spectra.medF_lowNE_HbT{NE_order(i).Runs}),2);
    subAvg.FigE8.medF_highNE_HbT(:,i) = mean(cat(2,spectra.medF_highNE_HbT{NE_order(i).Runs}),2);
end

plotBM = refBM;
plotBM(:,1:300) = NaN;

fig_savePath = fullfile(savePath,'ExtDataFig8');
[~, ~, ~] = mkdir(fig_savePath);

%% extract NE peaks

tBef = 50;
tAft = 50;

[peaks,trials] = f_findPeaks(spectra.NE,[0,0.02],10,0.005,[tBef,tAft]);

%% Fig Spectra A

tmp_trials = cell(numel(spectra.NE),1);
tmp_subAvg = NaN(tBef*10+tAft*10+1,numel(NE_order));
for i = 1:numel(spectra.NE)
    for idx = 1:numel(peaks{i})
        tmp_trials{i}(:,idx) = spectra.NE{i}(peaks{i}(idx)-tBef*10:peaks{i}(idx)+tAft*10);
    end
end
for i = 1:numel(NE_order)
    tmp_subAvg(:,i) = mean(cat(2,tmp_trials{NE_order(i).Runs}),2);
end

f = figure(Position=[100,100,1000,200]);
meanSig = mean(tmp_subAvg,2);
error = std(tmp_subAvg,0,2)/sqrt(M_NE);
f_plotLineError(-50:0.1:50,meanSig,error,color=c_GRAB,log=0,lineWidth=3);
% axis off;
ylim([1.5*[-1,1]]);

saveas(f, fullfile(fig_savePath, 'ExtDataFig8_A.svg'));
t = (-50:0.1:50)';

T = table(t,meanSig,error, ...
    VariableNames={'Time','mean_NE','SEM_NE'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig8_A.csv'));

%% Fig Spectra B
tmp_trials = cell(numel(spectra.NE),1);
tmp_subAvg = NaN(tBef*10+tAft*10+1,300,12,numel(NE_order));

for i = 1:numel(spectra.NE)
    for idx = 1:numel(peaks{i})
        tmp_trials{i}(:,:,:,idx) = spectra.SPG_Ca{i}(peaks{i}(idx)-tBef*10:peaks{i}(idx)+tAft*10,:,:);
    end
end

for i = 1:numel(NE_order)
    tmp_subAvg = mean(cat(4,tmp_trials{NE_order(i).Runs}),4);
end

f = figure(Position=[100,100,1000,300]);
imagesc(flipud(mean(tmp_subAvg,[3,4])'));
axis off;
colormap cmpinf;
clim([0 1.25]);
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig8_B.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

%% Fig Spectra C
tmp_trials = cell(numel(spectra.NE),1);
tmp_subAvg = NaN(tBef*10+tAft*10+1,300,12,numel(NE_order));

for i = 1:numel(spectra.NE)
    for idx = 1:numel(peaks{i})
        tmp_trials{i}(:,:,:,idx) = spectra.SPG_HbT{i}(peaks{i}(idx)-tBef*10:peaks{i}(idx)+tAft*10,:,:);
    end
end

for i = 1:numel(NE_order)
    tmp_subAvg = mean(cat(4,tmp_trials{NE_order(i).Runs}),4);
end

f = figure(Position=[100,100,1000,300]);
imagesc(flipud(mean(tmp_subAvg,[3,4])'));
axis off;
colormap cmpinf;
clim([0 0.2]);
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig8_C.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

%% Fig Spectra D

f = figure;
meanSig1 = mean(subAvg.FigE8.low_NE_Ca,2);
error1 = std(subAvg.FigE8.low_NE_Ca,0,2)/sqrt(M_NE);
f_plotLineError(spectra.fr,spectra.fr.*meanSig1,spectra.fr.*error1,color=c_darkCyan,log=1,lineWidth=3);
meanSig2 = mean(subAvg.FigE8.high_NE_Ca,2);
error2 = std(subAvg.FigE8.high_NE_Ca,0,2)/sqrt(M_NE);
f_plotLineError(spectra.fr,spectra.fr.*meanSig2,spectra.fr.*error2,color=c_Orange,log=1,lineWidth=3);
set(gca,'YScale','linear','FontSize',14);
xlim([0.01 5]);
xlabel('F (Hz)');
ylabel('Power (normalized by F)');
legend('','low NE','','high NE');

saveas(f, fullfile(fig_savePath, 'ExtDataFig8_D1.svg'));

T = table(spectra.fr,spectra.fr.*meanSig1,spectra.fr.*meanSig2,spectra.fr.*error1,spectra.fr.*error2, ...
    VariableNames={'F (Hz)','mean_low_NE','mean_high_NE','SEM_low_NE','SEM_high_NE'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig8_D1.csv'));

f = figure;
meanSig1 = mean(subAvg.FigE8.low_NE_HbT,2);
error1 = std(subAvg.FigE8.low_NE_HbT,0,2)/sqrt(M_NE);
f_plotLineError(spectra.fr,spectra.fr.*meanSig1,spectra.fr.*error1,color=c_darkCyan,log=1,lineWidth=3);
meanSig2 = mean(subAvg.FigE8.high_NE_HbT,2);
error2 = std(subAvg.FigE8.high_NE_HbT,0,2)/sqrt(M_NE);
f_plotLineError(spectra.fr,spectra.fr.*meanSig2,spectra.fr.*error2,color=c_Orange,log=1,lineWidth=3);
set(gca,'YScale','linear','FontSize',14);
xlim([0.01 5]);
xlabel('F (Hz)');
ylabel('Power (normalized by F)');
legend('','low NE','','high NE');

saveas(f, fullfile(fig_savePath, 'ExtDataFig8_D2.svg'));

T = table(spectra.fr,spectra.fr.*meanSig1,spectra.fr.*meanSig2,spectra.fr.*error1,spectra.fr.*error2, ...
    VariableNames={'F (Hz)','mean_low_NE','mean_high_NE','SEM_low_NE','SEM_high_NE'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig5_D2.csv'));

%% Fig Spectra E

f = figure;f_plotAllenMap(mean(subAvg.FigE8.lowF_lowNE_Ca,2),cmp=cmpinf,mask=plotBM,clim=[0,1]);
colorbar off;
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig8_E1.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

f = figure;f_plotAllenMap(mean(subAvg.FigE8.lowF_highNE_Ca,2),cmp=cmpinf,mask=plotBM,clim=[0,1]);
colorbar off;
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig8_E2.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

f = figure;f_plotAllenMap(mean(subAvg.FigE8.lowF_highNE_Ca,2)-mean(subAvg.FigE8.lowF_lowNE_Ca,2),cmp=cmpbbr,mask=plotBM,clim=0.3*[-1,1]);
colorbar off;
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig8_E3.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

%% Fig Spectra E 2

f = figure;f_plotAllenMap(mean(subAvg.FigE8.medF_lowNE_Ca,2),cmp=cmpinf,mask=plotBM,clim=[0,1]);
colorbar off;
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig8_E4.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

f = figure;f_plotAllenMap(mean(subAvg.FigE8.medF_highNE_Ca,2),cmp=cmpinf,mask=plotBM,clim=[0,1]);
colorbar off;
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig8_E5.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

f = figure;f_plotAllenMap(mean(subAvg.FigE8.medF_highNE_Ca,2)-mean(subAvg.FigE8.medF_lowNE_Ca,2),cmp=cmpbbr,mask=plotBM,clim=0.3*[-1,1]);
colorbar off;
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig8_E6.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

%%

f = figure;f_plotAllenMap(mean(subAvg.FigE8.highF_lowNE_Ca,2),cmp=cmpinf,mask=plotBM,clim=[0,0.4]);
colorbar off;
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig8_E7.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

f = figure;f_plotAllenMap(mean(subAvg.FigE8.highF_highNE_Ca,2),cmp=cmpinf,mask=plotBM,clim=[0,0.4]);
colorbar off;
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig8_E8.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

f = figure;f_plotAllenMap(mean(subAvg.FigE8.highF_highNE_Ca,2)-mean(subAvg.FigE8.highF_lowNE_Ca,2),cmp=cmpbbr,mask=plotBM,clim=0.1*[-1,1]);
colorbar off;
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig8_E9.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

%% Fig Spectra F

f = figure;f_plotAllenMap(mean(subAvg.FigE8.lowF_lowNE_HbT,2),cmp=cmpinf,mask=plotBM,clim=[0,0.15]);
colorbar off;
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig8_F1.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

f = figure;f_plotAllenMap(mean(subAvg.FigE8.lowF_highNE_HbT,2),cmp=cmpinf,mask=plotBM,clim=[0,0.15]);
colorbar off;
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig8_F2.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

f = figure;f_plotAllenMap(mean(subAvg.FigE8.lowF_highNE_HbT,2)-mean(subAvg.FigE8.lowF_lowNE_HbT,2),cmp=cmpbbr,mask=plotBM,clim=0.02*[-1,1]);
colorbar off;
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig8_F3.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

%%

f = figure;f_plotAllenMap(mean(subAvg.FigE8.medF_lowNE_HbT,2),cmp=cmpinf,mask=plotBM,clim=[0,0.15]);
colorbar off;
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig8_F4.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

f = figure;f_plotAllenMap(mean(subAvg.FigE8.medF_highNE_HbT,2),cmp=cmpinf,mask=plotBM,clim=[0,0.15]);
colorbar off;
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig8_F5.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

f = figure;f_plotAllenMap(mean(subAvg.FigE8.medF_highNE_HbT,2)-mean(subAvg.FigE8.medF_lowNE_HbT,2),cmp=cmpbbr,mask=plotBM,clim=0.03*[-1,1]);
colorbar off;
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig8_F6.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

%%
function [peaks,trials] = f_findPeaks(NE,fwin,fr,th,t)

N = numel(NE);
peaks = cell(N,1);
trials = zeros(N,1);

for i = 1:N
    
    if isempty(NE{i})
        continue
    end

    sig = f_bpf(NE{i},fwin,fr);
    diffNE = diff(sig);
    transitions = diffNE>th;

    if sum(transitions)
        tIdx = find(diff(transitions)==1);
        if isempty(tIdx)
            continue
        end
        tDown = find(diff(transitions)==-1);
        tDown(tDown<tIdx(1)) = [];
        if numel(tIdx)>numel(tDown)
            tIdx(numel(tDown)+1:end) = [];
        end

        transitions = round(mean([tIdx tDown],2));
        transitions(transitions<t(1)*10+1) = [];
        transitions(transitions>numel(sig)-t(2)*10) = [];
        peaks{i} = transitions;
        trials(i) = numel(transitions);
    end
end

end


close all

% load data
parentDir = f_path();
filePath = fullfile(parentDir,'Figures/results/Blocker/blocker_metadata.mat');
load(filePath);

load(fullfile(parentDir,'Figures/plot_types/refAllen.mat'));
savePath = fullfile(parentDir,'results');

fig_savePath = fullfile(savePath,'ExtDataFig6');
[~, ~, ~] = mkdir(fig_savePath);

plotBM = refBM;
plotBM(:,1:300) = NaN;

%% order metadata

[E6_sessions,E6_ImagingOrder,E6_log] = f_orderMeta(metadata);

%% create placeholder matrices

n_fr = numel(metadata.Thy1_302.d24_11_18.Run01.fr);

IRF_width = 101;

comb = struct;

comb.BM = cell(numel(E6_sessions),1);
comb.allen = cell(numel(E6_sessions),1);

comb.IRFx1.perf = cell(numel(E6_sessions),1);
comb.IRFx1_var.perf = cell(numel(E6_sessions),1);

comb.IRFx1.IRF = NaN(IRF_width,numel(E6_sessions));
comb.IRFx1_var.IRF = NaN(IRF_width,12,numel(E6_sessions));

comb.SPC.HbT = NaN(n_fr,numel(E6_sessions));
comb.COH.Ca_HbT = NaN(n_fr,numel(E6_sessions));
comb.XC.Ca_HbT = NaN(101,numel(E6_sessions));

comb.COH.Ca_HbT_allen = NaN(n_fr,12,numel(E6_sessions));

mice = fieldnames(metadata);

m = numel(mice);
runningIdx = 1;

for mIdx = 1:m
    date = fieldnames(metadata.(mice{mIdx}));
    d = numel(date);

    for dIdx = 1:d
        run = fieldnames(metadata.(mice{mIdx}).(date{dIdx}));
        
        tmpDate = date{dIdx}(2:end);
        tmpDate = strrep(tmpDate,'_','-');
        % parcellation = load(files.allen);parcellation = parcellation.parcellation;

        r = numel(run);

        for rIdx = 1:r
            interRun = numel(metadata.(mice{mIdx}).(date{dIdx}).(run{rIdx}));

            for irIdx = 1:interRun

                runData = metadata.(mice{mIdx}).(date{dIdx}).(run{rIdx})(irIdx);
                fr = runData.fr;

                comb.BM{runningIdx} = runData.brain_mask;

                comb.allen{runningIdx} = runData.parcellation;

                comb.IRFx1.perf{runningIdx} = runData.IRFx1.inv.perf;
                comb.IRFx1_var.perf{runningIdx} = runData.IRFx1.var.perf;

                comb.IRFx1.IRF(:,runningIdx) = runData.IRFx1.inv.IRF;

                comb.SPC.HbT(:,runningIdx) = runData.SPC.overall.HbT/sum(runData.SPC.overall.HbT);
                comb.COH.Ca_HbT(:,runningIdx) = runData.COH.overall.C.Ca_HbT;
                comb.XC.Ca_HbT(:,runningIdx) = runData.XC.overall.Ca_HbT;
                comb.COH.Ca_HbT_allen(:,:,runningIdx) = runData.COH.allen.C.Ca_HbT';

                for i = 1:12
                    comb.IRFx1_var.IRF(:,i,runningIdx) = squeeze(mean(runData.IRFx1.var.IRF.*runData.allen(:,:,i),[1 2],'omitnan'));
                end

                runningIdx = runningIdx+1;
            end
        end
    end
end

%% register to allen atlas

comb.BM = f_regImages(comb.BM,refParcellation,comb.allen,1);
E6_BM = comb.BM;

comb.IRFx1.perf = f_regImages(comb.IRFx1.perf,refParcellation,comb.allen,0).*E6_BM;
comb.IRFx1_var.perf = f_regImages(comb.IRFx1_var.perf,refParcellation,comb.allen,0).*E6_BM;

%% combine mice

subAvg.FigE6 = struct;
subAvg.FigE6.IRFx1.perf = NaN([size(refBM) m]);
subAvg.FigE6.IRFx1_var.perf = NaN([size(refBM) m]);

subAvg.FigE6.IRFx1.IRF = NaN(IRF_width,m);
subAvg.FigE6.IRFx1_var.IRF = NaN(IRF_width,12,m);

subAvg.FigE6.SPC.HbT = NaN(n_fr,m);
subAvg.FigE6.COH.Ca_HbT = NaN(n_fr,m);
subAvg.FigE6.COH.Ca_HbT_allen = NaN(n_fr,12,m);
subAvg.FigE6.XC.Ca_HbT = NaN(101,m);

remove = 1:2;
for idx = 1:m
    tIdx = zeros(numel(E6_log),1);
    tIdx(E6_ImagingOrder(idx).runs) = 1;
    tIdx = logical(tIdx.*(~ismember([E6_log.Run],remove))');
    subAvg.FigE6.IRFx1.perf(:,:,idx) = mean(comb.IRFx1.perf(:,:,tIdx).*comb.BM(:,:,tIdx),3,'omitnan');
    subAvg.FigE6.IRFx1_var.perf(:,:,idx) = mean(comb.IRFx1_var.perf(:,:,tIdx).*comb.BM(:,:,tIdx),3,'omitnan');
    
    subAvg.FigE6.IRFx1.IRF(:,idx) = mean(comb.IRFx1.IRF(:,tIdx),2);
    subAvg.FigE6.IRFx1_var.IRF(:,:,idx) = mean(comb.IRFx1_var.IRF(:,:,tIdx),3);

    subAvg.FigE6.SPC.HbT(:,idx) = mean(comb.SPC.HbT(:,tIdx),2);
    subAvg.FigE6.COH.Ca_HbT(:,idx) = mean(comb.COH.Ca_HbT(:,tIdx),2);
    subAvg.FigE6.COH.Ca_HbT_allen(:,:,idx) = mean(comb.COH.Ca_HbT_allen(:,:,tIdx),3);
    subAvg.FigE6.XC.Ca_HbT(:,idx) = mean(comb.XC.Ca_HbT(:,tIdx),2);
end

subAvg.post = subAvg.FigE6;

remove = 3:8;
for idx = 1:m
    tIdx = zeros(numel(E6_log),1);
    tIdx(E6_ImagingOrder(idx).runs) = 1;
    tIdx = logical(tIdx.*(~ismember([E6_log.Run],remove))');
    subAvg.FigE6.IRFx1.perf(:,:,idx) = mean(comb.IRFx1.perf(:,:,tIdx).*comb.BM(:,:,tIdx),3,'omitnan');
    subAvg.FigE6.IRFx1_var.perf(:,:,idx) = mean(comb.IRFx1_var.perf(:,:,tIdx).*comb.BM(:,:,tIdx),3,'omitnan');
    
    subAvg.FigE6.IRFx1.IRF(:,idx) = mean(comb.IRFx1.IRF(:,tIdx),2);
    subAvg.FigE6.IRFx1_var.IRF(:,:,idx) = mean(comb.IRFx1_var.IRF(:,:,tIdx),3);

    subAvg.FigE6.SPC.HbT(:,idx) = mean(comb.SPC.HbT(:,tIdx),2);
    subAvg.FigE6.COH.Ca_HbT(:,idx) = mean(comb.COH.Ca_HbT(:,tIdx),2);
    subAvg.FigE6.COH.Ca_HbT_allen(:,:,idx) = mean(comb.COH.Ca_HbT_allen(:,:,tIdx),3);
    subAvg.FigE6.XC.Ca_HbT(:,idx) = mean(comb.XC.Ca_HbT(:,tIdx),2);
end

subAvg.baseline = subAvg.FigE6;

%% plot FigBlocker B

signalPath = fullfile(parentDir,'Figures/results/Blocker/rawData.mat');
signals = load(signalPath);

Ca = signals.data1.rfp_HD(:,5);
HbT = signals.data1.HbT(:,5);
Pupil = signals.data1.Pupil;
Whisking = signals.data1.Whisking;
Movement = signals.data1.Accelerometer;

t = 0.1:0.1:600;

f = figure;
tiledlayout(4,1);

nexttile;hold on;
plot(t,Ca,color=c_Ca);
plot(t,HbT-10,color=c_HbT);
plot(60*[1,1],[0,10],'-k',lineWidth=2);
xlim([50 250]);
box off;

nexttile;hold on;
plot(t,Pupil,color=c_pupil);
plot(60*[1,1],[0,0.4],'-k',lineWidth=2);
xlim([50 250]);
box off;

nexttile;
plot(t,Whisking,color=[0,0.7,0.7]);
xlim([50 250]);
box off;

nexttile;
plot(t,Movement,color=[0,0,0]);
xlim([50 250]);
box off;

idx = 500:2500;

saveas(f, fullfile(fig_savePath, 'ExtDataFig6_B.svg'));

T = table(t(idx)',Ca(idx),HbT(idx),Pupil(idx),Whisking(idx),Movement(idx), ...
    VariableNames={'Time','Ca','HbT','Pupil','Whisking','Accelerometer'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig6_B.csv'));

%% plot FigBlocker C

signalPath = fullfile(parentDir,'Figures/results/Blocker/rawData.mat');
signals = load(signalPath);

Ca = signals.data2.rfp_HD(:,5);
HbT = signals.data2.HbT(:,5);
Pupil = signals.data2.Pupil;
Whisking = signals.data2.Whisking;
Movement = signals.data2.Accelerometer;

t = 0.1:0.1:600;

f = figure;
tiledlayout(4,1);

nexttile;hold on;
plot(t,Ca,color=c_Ca);
plot(t,HbT-10,color=c_HbT);
plot(260*[1,1],[0,10],'-k',lineWidth=2);
xlim([250 450]);
box off;

nexttile;hold on;
plot(t,Pupil,color=c_pupil);
plot(260*[1,1],[0,0.4],'-k',lineWidth=2);
xlim([250 450]);
box off;

nexttile;
plot(t,Whisking,color=[0,0.7,0.7]);
xlim([250 450]);
box off;

nexttile;
plot(t,Movement,color=[0,0,0]);
xlim([250 450]);
box off;

idx = 2500:4500;

saveas(f, fullfile(fig_savePath, 'ExtDataFig6_C.svg'));

T = table(t(idx)',Ca(idx),HbT(idx),Pupil(idx),Whisking(idx),Movement(idx), ...
    VariableNames={'Time','Ca','HbT','Pupil','Whisking','Accelerometer'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig6_C.csv'));

%% plot FigBlocker D

sig = subAvg.baseline.COH.Ca_HbT;
meanSig1 = mean(sig(2:end,:),2);
error1 = std(sig(2:end,:),0,2)/sqrt(m);

sig = subAvg.post.COH.Ca_HbT;
meanSig2 = mean(sig(2:end,:),2);
error2 = std(sig(2:end,:),0,2)/sqrt(m);

tmpfr = fr(2:end);

f = figure;
f_plotLineError(tmpfr,meanSig1,error1,color=c_Orange,log=1);
f_plotLineError(tmpfr,meanSig2,error2,color=c_darkCyan,log=1);
xlim([0.05, 5]);
xlabel('Frequency (Hz)');
ylabel('C');
set(gca,'YScale','linear','FontSize',14);
ylim([0,0.8]);

saveas(f, fullfile(fig_savePath, 'ExtDataFig6_D.svg'));
T = table(tmpfr',meanSig1,meanSig2,error1,error2, ...
    VariableNames={'F (Hz)','mean_baseline','mean_PPA','SEM_baseline','SEM_PPA'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig6_D.csv'));

%% plot FigBlocker E

fr_Range = [0 0.1;0.1 0.5];
fr_Idx = 1;
fr_Range = fr_Range(fr_Idx,:);

[~,fr_Range] = min(abs(fr'-fr_Range));

sig = squeeze(mean(subAvg.baseline.COH.Ca_HbT_allen(fr_Range(1):fr_Range(2),:,:),1));
f = figure;f_plotAllenMap(mean(sig,2),cmp=cmpinf,mask=plotBM,clim=[0,1]);
colorbar off;
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig6_E_baseLow.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

sig = squeeze(mean(subAvg.post.COH.Ca_HbT_allen(fr_Range(1):fr_Range(2),:,:),1));
f = figure;f_plotAllenMap(mean(sig,2),cmp=cmpinf,mask=plotBM,clim=[0,1]);
colorbar off;
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig6_E_postLow.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

sig = squeeze(mean(subAvg.post.COH.Ca_HbT_allen(fr_Range(1):fr_Range(2),:,:),1))-squeeze(mean(subAvg.baseline.COH.Ca_HbT_allen(fr_Range(1):fr_Range(2),:,:),1));
f = figure;f_plotAllenMap(mean(sig,2),cmp=cmpbbr,mask=plotBM,clim=0.3*[-1,1]);
colorbar off;
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig6_E_diffLow.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);


fr_Range = [0 0.1;0.1 0.5];
fr_Idx = 2;
fr_Range = fr_Range(fr_Idx,:);

[~,fr_Range] = min(abs(fr'-fr_Range));

sig = squeeze(mean(subAvg.baseline.COH.Ca_HbT_allen(fr_Range(1):fr_Range(2),:,:),1));
f = figure;f_plotAllenMap(mean(sig,2),cmp=cmpinf,mask=plotBM,clim=[0,1]);
colorbar off;
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig6_E_baseHigh.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

sig = squeeze(mean(subAvg.post.COH.Ca_HbT_allen(fr_Range(1):fr_Range(2),:,:),1));
f = figure;f_plotAllenMap(mean(sig,2),cmp=cmpinf,mask=plotBM,clim=[0,1]);
colorbar off;
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig6_E_postHigh.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

sig = squeeze(mean(subAvg.post.COH.Ca_HbT_allen(fr_Range(1):fr_Range(2),:,:),1))-squeeze(mean(subAvg.baseline.COH.Ca_HbT_allen(fr_Range(1):fr_Range(2),:,:),1));
f = figure;f_plotAllenMap(mean(sig,2),cmp=cmpbbr,mask=plotBM,clim=0.3*[-1,1]);
colorbar off;
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig6_E_diffHigh.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

%% plot FigBlocker F

sig = subAvg.baseline.XC.Ca_HbT;
meanSig1 = mean(sig,2);
error1 = std(sig,0,2)/sqrt(m);

sig = subAvg.post.XC.Ca_HbT;
meanSig2 = mean(sig,2);
error2 = std(sig,0,2)/sqrt(m);

c = c_darkCyan;
% c = c_Orange;

t = 5:-0.1:-5;

f = figure;
f_plotLineError(t,meanSig1,error1,color=c_Orange,lineWidth=3);
f_plotLineError(t,meanSig2,error2,color=c_darkCyan,lineWidth=3);
xlim([-5 5]);
ylim([-0.2, 0.4]);
xlabel('Time (s)');
ylabel('r');
set(gca,'FontSize',14);

xlim([-5 5]);

saveas(f, fullfile(fig_savePath, 'ExtDataFig6_F.svg'));
T = table(t',meanSig1,meanSig2,error1,error2, ...
    VariableNames={'Time','mean_baseline','mean_PPA','SEM_baseline','SEM_PPA'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig6_F.csv'));

%% plot FigBlocker G

sig = squeeze(mean(subAvg.baseline.IRFx1_var.IRF(:,[4,5],:),2));
meanSig1 = mean(sig,2);
error1 = std(sig,0,2)/sqrt(m);

t = 0:0.1:10;
f = figure;
f_plotLineError(t,meanSig1,error1,color=c_Orange);
xlim([0 7]);
xlabel('Time (s)');
ylabel('a.u.');

sig = squeeze(mean(subAvg.post.IRFx1_var.IRF(:,[4,5],:),2));
meanSig2 = mean(sig,2);
error2 = std(sig,0,2)/sqrt(m);

f_plotLineError(t,meanSig2,error2,color=c_darkCyan);
xlim([0 7]);
xlabel('Time (s)');
ylabel('a.u.');

plot([0,7],[0,0],'-k',lineWidth=2);
ylim([-0.04 0.1]);

saveas(f, fullfile(fig_savePath, 'ExtDataFig6_G1.svg'));
T = table(t',meanSig1,meanSig2,error1,error2, ...
    VariableNames={'Time','mean_baseline','mean_PPA','SEM_baseline','SEM_PPA'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig6_G1.csv'));


sig = squeeze(mean(subAvg.baseline.IRFx1_var.IRF(:,2,:),2));
meanSig1 = mean(sig,2);
error1 = std(sig,0,2)/sqrt(m);

t = 0:0.1:10;
f = figure;
f_plotLineError(t,meanSig1,error1,color=c_Orange);
xlim([0 7]);
xlabel('Time (s)');
ylabel('a.u.');

sig = squeeze(mean(subAvg.post.IRFx1_var.IRF(:,2,:),2));
meanSig2 = mean(sig,2);
error2 = std(sig,0,2)/sqrt(m);

f_plotLineError(t,meanSig2,error2,color=c_darkCyan);
xlim([0 7]);
xlabel('Time (s)');
ylabel('a.u.');

plot([0,7],[0,0],'-k',lineWidth=2);
ylim([-0.04 0.1]);

saveas(f, fullfile(fig_savePath, 'ExtDataFig6_G2.svg'));
T = table(t',meanSig1,meanSig2,error1,error2, ...
    VariableNames={'Time','mean_baseline','mean_PPA','SEM_baseline','SEM_PPA'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig6_G2.csv'));

%% plot FigBlocker H

sig = mean(subAvg.baseline.IRFx1_var.perf,3,'omitnan');
f = figure(Position = [100, 100, 500, 400]);
imagesc(sig.*plotBM,AlphaData=plotBM);
colormap cmpvir;
clim([0, 1]);
axis image off;
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig6_H_base.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

sig = mean(subAvg.post.IRFx1_var.perf,3,'omitnan');
f = figure(Position = [100, 100, 500, 400]);
imagesc(sig.*plotBM,AlphaData=plotBM);
colormap cmpvir;
clim([0, 1]);
axis image off;
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig6_H_post.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

sig = mean(subAvg.post.IRFx1_var.perf,3,'omitnan')-mean(subAvg.baseline.IRFx1_var.perf,3,'omitnan');
f = figure(Position = [100, 100, 500, 400]);
imagesc(sig.*plotBM,AlphaData=plotBM);
colormap cmpbbr;
clim(0.5*[-1, 1]);
axis image off;
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig6_H_diff.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

%% plot FigBlocker I

tmpMask = refParcellation.Masks(:,:,3,2).*refBM;
tmpMask(tmpMask==0) = NaN;

post = squeeze(mean(subAvg.post.IRFx1_var.perf.*tmpMask,[1 2],'omitnan'));
baseline = squeeze(mean(subAvg.baseline.IRFx1_var.perf.*tmpMask,[1 2],'omitnan'));

barData = {};
barData{1} = baseline;
barData{2} = post;

f = figure;
[dataMean, dataSEM] = f_plotBar(barData,colors=[c_Orange;c_darkCyan],ylabel='r')
ylim([0, 1]);

[h,p] = f_kstest(barData,0.05);

saveas(f, fullfile(fig_savePath, 'ExtDataFig6_I.svg'));
T = table(mice,baseline,post, ...
    VariableNames={'Mouse','baseline','post'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig6_I.csv'));

%%
% %% plot spectra
% 
% sig = subAvg.post.SPC.HbT;
% sig(1,:) = [];
% sig = log10(sig);
% 
% meanSig = mean(sig,2);
% 
% tmp_fr = log10(fr(2:end));
% err = std(sig,0,2)/sqrt(sum(idx));
% 
% c = c_darkCyan;
% % c = c_Orange;
% 
% % f = figure;hold on;
% fill([tmp_fr fliplr(tmp_fr)],[meanSig+err;flipud(meanSig-err)],c,'EdgeColor','none',FaceAlpha=0.3);
% plot(tmp_fr,meanSig,Color=c,LineWidth=3);
% xlim(log10([0.05 5]));
% axis off;
% 
% %%
% scatter(ones(sum(idx),1)*1,sig,'filled','XJitter','randn','XJitterWidth',0.1,MarkerFaceColor=c_darkCyan,MarkerFaceAlpha=0.3);
% z = sum(sig<0);
% scatter(ones(z,1)*1,zeros(z,1),60,'o','XJitter','randn','XJitterWidth',0.4,'MarkerEdgeColor',c_darkCyan);
% 
% er.Color = [0 0 0];
% er.LineStyle = 'none';
% er.LineWidth = 2;
% b.FaceColor = 'flat';
% b.CData = [1 1 1];
% b.ShowBaseLine = 'off';
% b.EdgeColor = c_darkCyan;
% b.LineWidth = 2;
% axis off;
% ylim([0 1]);
% 
% %% plot A and B
% 
% imAlpha = f_nan2zero(refBM);
% imAlpha(:,1:300) = 0;
% 
% f = figure;
% imagesc(mean(subAvg.LR_unf.A,3,'omitnan'),AlphaData=imAlpha);
% axis image off;
% clim(1*[-1 1]);
% colormap cmpbwr;
% 
% %% plot LR lag
% 
% sig = subAvg.LR.lag;
% 
% err = std(sig,0,2)/sqrt(m);
% meanSig = mean(sig,2);
% 
% f = figure(Position=[100 100 375 400]);hold on;
% b(1) = bar(1,meanSig(1),0.6);
% b(2) = bar(2,meanSig(2),0.6);
% scatter(1*ones(m,1),sig(1,:),'filled','XJitter','randn','XJitterWidth',0.4,MarkerFaceColor=c_Ca,MarkerFaceAlpha=0.3);
% scatter(2*ones(m,1),sig(2,:),'filled','XJitter','randn','XJitterWidth',0.4,MarkerFaceColor=c_GRAB,MarkerFaceAlpha=0.3);
% 
% er = errorbar(1:2,meanSig,err);
% 
% er.Color = [0 0 0];
% er.LineStyle = 'none';
% er.LineWidth = 2;
% b(1).FaceColor = 'flat';
% b(1).CData = [1 1 1];
% b(1).ShowBaseLine = 'off';
% b(1).EdgeColor = c_Ca;
% b(1).LineWidth = 2;
% b(2).FaceColor = 'flat';
% b(2).CData = [1 1 1];
% b(2).ShowBaseLine = 'off';
% b(2).EdgeColor = c_GRAB;
% b(2).LineWidth = 2;
% ylim([-0.2 1]);
% 
% axis off;
% %% adjust transparency
% 
% filename = 'NE_reg_perf.jpg';
% exportgraphics(f,filename,'Resolution',300,'BackgroundColor',[1 1 1]);
% 
% img = imread(filename);
% 
% imAlpha = sum(img,3);
% imAlpha(imAlpha==765) = 0;
% imAlpha(find(imAlpha)) = 1;
% 
% imwrite(img,filename,'Alpha',imAlpha);
% 
% %% compare exposures
% 
% f = figure;
% imagesc(sum(comb.BM,3,'omitnan')/size(comb.BM,3),'AlphaData',refBM)
% axis off image;
% colormap cmpwvir;
% 
% %% create Fig2 Struct
% 
% Fig2 = struct;
% Fig2.trials = log;
% Fig2.data.B.A = comb.LR.A;
% Fig2.data.B.B = comb.LR.B;
% Fig2.data.C.lag = comb.LR.lag;
% Fig2.data.D.perf = comb.LR.perf;
% Fig2.data.F.A = comb.IRFx2.A;
% Fig2.data.F.B = comb.IRFx2.B;
% Fig2.data.G.IRF = comb.IRFx2.IRF;
% Fig2.data.H.perf = comb.IRFx2.perf;
% Fig2.data.J.LR = comb.s.LR.perf;
% Fig2.data.J.IRFx2 = comb.s.IRFx2.perf;
% 
% %% plot all mice
% 
% for mIdx = 6
%     % mIdx = 6;
%     plotImg = subAvg.IRFx2.perf(:,:,mIdx);
%     % plotImg = tmp(:,:,mIdx);
% 
% 
%     imAlpha = f_nan2zero(refBM);
%     imAlpha(:,1:300) = 0;
%     B = bwboundaries(imAlpha,'noholes');
% 
%     imAlpha = ~isnan(plotImg);
% 
%     f = figure;hold on;
%     imagesc(comb.IRFx2.perf(:,:,92),AlphaData=imAlpha);
%     clim([0 1]);
%     colormap cmpvir;
%     plot(B{1}(:,2),B{1}(:,1),'-k','LineWidth',3,Color=0*[1 1 1]); %[0.9 0.85 0.5]
%     axis image off;
%     set(gca,'YDir','reverse');
% 
%     % filename = [mice2{mIdx} '_inv_perf.jpg'];
%     % f_adjTrans(f,filename);
% end
% % close all;

%% utility functions

function allenMap = f_allenMap(masks,values)
    allenMap = NaN(size(masks(:,:,1)));
    for idx = 1:12
        tmp = logical(f_nan2zero(masks(:,:,idx)));
        allenMap(tmp(:)) = values(idx);
    end
end

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

function [idx] = f_sortGRAB(metadata,GRAB)
    mice = fieldnames(metadata);
    m = numel(mice);
    idx = zeros(m,1);


    for mIdx = 1:m
        date = fieldnames(metadata.(mice{mIdx}));
        date = date{1};
        run = fieldnames(metadata.(mice{mIdx}).(date));
        run = run{1};
        idx(mIdx) = strcmp(metadata.(mice{mIdx}).(date).(run)(1).GRAB,GRAB);
    end

end

function ImagingOrder = f_GRABidx(ImagingOrder,metadata,GRAB)
    n = numel(ImagingOrder);
    for i = 1:n
        dates = fieldnames(metadata.(ImagingOrder(i).mouse));
        runs = fieldnames(metadata.(ImagingOrder(i).mouse).(dates{1}));
        idx(i) = strcmp(metadata.(ImagingOrder(i).mouse).(dates{1}).(runs{1})(1).GRAB,GRAB);
    end
    ImagingOrder = ImagingOrder(idx);
end

function [subjectAvg] = f_subAvg(data,sessions)
    mice = unique({sessions.mouse});
    m = numel(mice);
    
    for i = 1:m
        order = strcmp({sessions.mouse},mice{i});
        subjectAvg(i) = mean(data(order,:),1);
    end
end

function [reg] = f_regAllen(data,allen,refAllen,isBM)
    
    refDim = size(refAllen.Masks);
    reg = zeros(refDim(1:3));

    n = numel(allen);

    for i = 1:n
        if nargin == 4
            reg(:,:,i) = f_ImgReg_allen(refAllen,allen{i},data{i},isBM);
        else
            reg(:,:,i) = f_ImgReg_allen(refAllen,allen{i},data{i});
        end
    end
end
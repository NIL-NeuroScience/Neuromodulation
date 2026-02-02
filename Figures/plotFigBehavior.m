% calculate subject averages

NE_order = order(NE_Idx);

tmp.Beh.R_rfp_HD_HbT = f_regImages(Behavior.R.rfp_HD_low_HbT_low,refParcellation,settings.allen_masks,0).*BM;
tmp.Beh.R_gfp_HD_HbT = f_regImages(Behavior.R.gfp_HD_low_HbT_low,refParcellation,settings.allen_masks,0).*BM;
tmp.Beh.R_rfp_HD_gfp_HD = f_regImages(Behavior.R.rfp_HD_low_gfp_HD_low,refParcellation,settings.allen_masks,0).*BM;
tmp.Beh.NE_IRF_perf = f_regImages(Behavior.NE_IRF_perf,refParcellation,settings.allen_masks,0).*BM;
tmp.Beh.SPG_rfp_HD = f_adjust_SPG(Behavior.SPG.rfp_HD,4097,1);
tmp.Beh.SPG_gfp_HD = f_adjust_SPG(Behavior.SPG.gfp_HD,4097,1);
tmp.Beh.SPG_HbT = f_adjust_SPG(Behavior.SPG.HbT,4097,1);
tmp.Beh.COH_rfp_HD_gfp_HD = f_adjust_SPG(Behavior.COH.rfp_HD_gfp_HD,4097,0);
tmp.Beh.COH_rfp_HD_HbT = f_adjust_SPG(Behavior.COH.rfp_HD_HbT,4097,0);
tmp.Beh.COH_gfp_HD_HbT = f_adjust_SPG(Behavior.COH.gfp_HD_HbT,4097,0);
tmp.Beh.PHI_rfp_HD_gfp_HD = f_adjust_SPG(Behavior.PHI.rfp_HD_gfp_HD,4097,0);
tmp.Beh.PHI_rfp_HD_HbT = f_adjust_SPG(Behavior.PHI.rfp_HD_HbT,4097,0);
tmp.Beh.PHI_gfp_HD_HbT = f_adjust_SPG(Behavior.PHI.gfp_HD_HbT,4097,0);

subAvg.Beh.XC_gfp_HD_HbT = NaN(201,M);
subAvg.Beh.XC_rfp_HD_HbT = NaN(201,M);
subAvg.Beh.XC_rfp_HD_gfp_HD = NaN(201,M);
subAvg.Beh.R_rfp_HD_HbT = NaN(500,600,M);
subAvg.Beh.R_gfp_HD_HbT = NaN(500,600,M);
subAvg.Beh.R_rfp_HD_gfp_HD = NaN(500,600,M);
subAvg.Beh.R_behavior = NaN(6,6,M);
subAvg.Beh.SPG_rfp_HD = NaN(4097,M);
subAvg.Beh.SPG_gfp_HD = NaN(4097,M);
subAvg.Beh.SPG_HbT = NaN(4097,M);
subAvg.Beh.COH_rfp_HD_gfp_HD = NaN(4097,M);
subAvg.Beh.COH_rfp_HD_HbT = NaN(4097,M);
subAvg.Beh.COH_gfp_HD_HbT = NaN(4097,M);
subAvg.Beh.PHI_rfp_HD_gfp_HD = NaN(4097,M);
subAvg.Beh.PHI_rfp_HD_HbT = NaN(4097,M);
subAvg.Beh.PHI_gfp_HD_HbT = NaN(4097,M);
subAvg.Beh.NE_IRF_perf = NaN(500,600,M);
subAvg.Beh.NE_IRF = NaN(151,M);
subAvg.Beh.GRAB_conn = NaN(12,12,M);
subAvg.Beh.XC_gfp_HD_HbT_allen = NaN(201,12,M);

for i = 1:M
    subAvg.Beh.XC_gfp_HD_HbT(:,i) = mean(cat(2,Behavior.XC.gfp_HD_HbT{order(i).Runs}),2);
    subAvg.Beh.XC_rfp_HD_HbT(:,i) = mean(cat(2,Behavior.XC.rfp_HD_HbT{order(i).Runs}),2);
    subAvg.Beh.XC_rfp_HD_gfp_HD(:,i) = mean(cat(2,Behavior.XC.rfp_HD_gfp_HD{order(i).Runs}),2);
    subAvg.Beh.R_rfp_HD_HbT(:,:,i) = mean(tmp.Beh.R_rfp_HD_HbT(:,:,order(i).Runs),3,'omitnan');
    subAvg.Beh.R_gfp_HD_HbT(:,:,i) = mean(tmp.Beh.R_gfp_HD_HbT(:,:,order(i).Runs),3,'omitnan');
    subAvg.Beh.R_rfp_HD_gfp_HD(:,:,i) = mean(tmp.Beh.R_rfp_HD_gfp_HD(:,:,order(i).Runs),3,'omitnan');
    subAvg.Beh.R_behavior(:,:,i) = mean(cat(3,Behavior.R.signals{order(i).Runs}),3,'omitnan');
    subAvg.Beh.SPG_rfp_HD(:,i) = mean(cat(2,tmp.Beh.SPG_rfp_HD{order(i).Runs}),2);
    subAvg.Beh.SPG_gfp_HD(:,i) = mean(cat(2,tmp.Beh.SPG_gfp_HD{order(i).Runs}),2);
    subAvg.Beh.SPG_HbT(:,i) = mean(cat(2,tmp.Beh.SPG_HbT{order(i).Runs}),2);
    subAvg.Beh.COH_rfp_HD_gfp_HD(:,i) = mean(cat(2,tmp.Beh.COH_rfp_HD_gfp_HD{order(i).Runs}),2);
    subAvg.Beh.COH_rfp_HD_HbT(:,i) = mean(cat(2,tmp.Beh.COH_rfp_HD_HbT{order(i).Runs}),2);
    subAvg.Beh.COH_gfp_HD_HbT(:,i) = mean(cat(2,tmp.Beh.COH_gfp_HD_HbT{order(i).Runs}),2);
    subAvg.Beh.PHI_rfp_HD_gfp_HD(:,i) = mean(cat(2,tmp.Beh.PHI_rfp_HD_gfp_HD{order(i).Runs}),2);
    subAvg.Beh.PHI_rfp_HD_HbT(:,i) = mean(cat(2,tmp.Beh.PHI_rfp_HD_HbT{order(i).Runs}),2);
    subAvg.Beh.PHI_gfp_HD_HbT(:,i) = mean(cat(2,tmp.Beh.PHI_gfp_HD_HbT{order(i).Runs}),2);
    subAvg.Beh.NE_IRF_perf(:,:,i) = mean(tmp.Beh.NE_IRF_perf(:,:,order(i).Runs),3,'omitnan');
    try 
        subAvg.Beh.NE_IRF(:,i) = mean(cat(2,Behavior.NE_IRF_IRF{order(i).Runs}),2);
        subAvg.Beh.XC_gfp_HD_HbT_allen(:,:,i) = mean(cat(3,GRAB_FC.gfp_HD_vs_HbT_low{order(i).Runs}),3);
    end
    subAvg.Beh.GRAB_conn(:,:,i) = mean(cat(3,GRAB_FC.FC_detrend{order(i).Runs}),3);
end

fr = Behavior.SPG.fr;

plotBM = refBM;
plotBM(:,1:300) = NaN;

fig_savePath = fullfile(savePath,'ExtDataFig2');
[~, ~, ~] = mkdir(fig_savePath);

%% Fig Behavior A

f = figure;
meanSig1 = mean(subAvg.Beh.SPG_rfp_HD,2);
error1 = std(subAvg.Beh.SPG_rfp_HD,0,2)/sqrt(M);
f_plotLineError(fr,meanSig1,error1,color=c_Ca,log=1,lineWidth=3);
meanSig2 = mean(subAvg.Beh.SPG_gfp_HD(:,ACh_Idx),2);
error2 = std(subAvg.Beh.SPG_gfp_HD(:,ACh_Idx),0,2)/sqrt(M_ACh);
f_plotLineError(fr,meanSig2,error2,color=c_Orange,log=1,lineWidth=3);
meanSig3 = mean(subAvg.Beh.SPG_gfp_HD(:,NE_Idx),2);
error3 = std(subAvg.Beh.SPG_gfp_HD(:,NE_Idx),0,2)/sqrt(M_NE);
f_plotLineError(fr,meanSig3,error3,color=c_GRAB,log=1,lineWidth=3);
meanSig4 = mean(subAvg.Beh.SPG_HbT,2);
error4 = std(subAvg.Beh.SPG_HbT,0,2)/sqrt(M);
f_plotLineError(fr,meanSig4,error4,color=c_HbT,log=1,lineWidth=3);

xlim([0.05,5]);
xlabel('F (Hz)');
ylabel('Normalized Power');
% legend('','Ca^2^+','','ACh','','NE','','HbT');
set(gca,'FontSize',14);
set(gcf,'Renderer','painters');

saveas(f, fullfile(fig_savePath, 'ExtDataFig2_A.svg'));
T = table(fr,meanSig1,meanSig2,meanSig3,meanSig4,error1,error2,error3,error4, ...
    VariableNames={'F (Hz)','mean_Ca','mean_ACh','mean_NE','mean_HbT', ...
    'SEM_Ca','SEM_ACh','SEM_NE','SEM_HbT'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig2_A.csv'));

%% Fig Behavior B

t = 10:-0.1:-10;

f = figure;
meanSig1 = mean(subAvg.Beh.XC_rfp_HD_gfp_HD(:,NE_Idx),2);
error1 = std(subAvg.Beh.XC_rfp_HD_gfp_HD(:,NE_Idx),0,2)/sqrt(M_NE);
f_plotLineError(t,meanSig1,error1,color=c_GRAB,lineWidth=3);
meanSig2 = mean(subAvg.Beh.XC_rfp_HD_gfp_HD(:,ACh_Idx),2);
error2 = std(subAvg.Beh.XC_rfp_HD_gfp_HD(:,ACh_Idx),0,2)/sqrt(M_ACh);
f_plotLineError(t,meanSig2,error2,color=c_Orange,lineWidth=3);
xlim([-5 5]);
xlabel('Time (s)');
ylabel('r');
set(gca,'FontSize',14);
title('x vs. Ca^2^+');
legend('','NE','','ACh');

saveas(f, fullfile(fig_savePath, 'ExtDataFig2_B.svg'));
T = table(t',meanSig1,meanSig2,error1,error2, ...
    VariableNames={'Lag','mean_NE','mean_ACh', ...
    'SEM_NE','SEM_ACh'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig2_B.csv'));

%% Fig Behavior C

R_Beh = mean(subAvg.Beh.R_behavior,3);
gfp_R = mean(subAvg.Beh.R_behavior(:,:,ACh_Idx),3);
R_Beh(2,:) = gfp_R(2,:);
R_Beh(:,2) = gfp_R(2,:);
gfp_R = mean(subAvg.Beh.R_behavior(:,:,NE_Idx),3);
R_Beh(4:7,:) = R_Beh(3:6,:);
R_Beh(:,4:7) = R_Beh(:,3:6);
R_Beh(3,[1,3,4,5,6,7]) = gfp_R(2,:);
R_Beh([1,3,4,5,6,7],3) = gfp_R(2,:);
R_Beh(2,3) = 0;R_Beh(3,2) = 0;

barData = {};
barData{1} = squeeze(subAvg.Beh.R_behavior(2,4,ACh_Idx));
barData{2} = squeeze(subAvg.Beh.R_behavior(2,5,ACh_Idx));
barData{3} = squeeze(subAvg.Beh.R_behavior(2,6,ACh_Idx));
barData{4} = squeeze(subAvg.Beh.R_behavior(2,4,NE_Idx));
barData{5} = squeeze(subAvg.Beh.R_behavior(2,5,NE_Idx));
barData{6} = squeeze(subAvg.Beh.R_behavior(2,6,NE_Idx));
barData{7} = squeeze(subAvg.Beh.R_behavior(1,4,:));
barData{8} = squeeze(subAvg.Beh.R_behavior(1,5,:));
barData{9} = squeeze(subAvg.Beh.R_behavior(1,6,:));
barData{10} = squeeze(subAvg.Beh.R_behavior(3,4,:));
barData{11} = squeeze(subAvg.Beh.R_behavior(3,5,:));
barData{12} = squeeze(subAvg.Beh.R_behavior(3,6,:));

f = figure;
[dataMean, dataSEM] = f_plotBar(barData,colors=repmat([c_pupil;0,0.7,0.7;0,0,0],4,1),legend={'Pupil diameter','Whisking','Movement'},ylabel='r',title='Model Performance Comparison',ylim=[0,0.8])
saveas(f, fullfile(fig_savePath, 'ExtDataFig2_C1.svg'));

labels = {'ACh_pupil','ACh_whisking','ACh_movement', ...
    'NE_pupil','NE_whisking','NE_movement', ...
    'Ca_pupil','Ca_whisking','Ca_movement', ...
    'HbT_pupil','HbT_whisking','HbT_movement'};
T = table({order.Mouse}',VariableNames={'Mouse'});

for i = 1:3
    tmp = NaN(numel({order.Mouse}),1);
    tmp(ACh_Idx) = barData{i};
    T.(labels{i}) = tmp;
end
for i = 4:6
    tmp = NaN(numel({order.Mouse}),1);
    tmp(NE_Idx) = barData{i};
    T.(labels{i}) = tmp;
end
for i = 7:12
    T.(labels{i}) = barData{i};
end
writetable(T, fullfile(fig_savePath, 'ExtDataFig2_C1.csv'));

f = figure;
imagesc(R_Beh);
colormap cmpbbr;
clim([-1 1]);
axis image off;
c = colorbar;
c.Label.String = 'r';
set(gca,'FontSize',14);
colorbar off;
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig2_C2.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);
writetable(table(R_Beh), fullfile(fig_savePath, 'ExtDataFig2_C2.csv'));


%% Fig Behavior D

f = figure;
f_plotMap(mean(subAvg.Beh.R_rfp_HD_gfp_HD(:,:,ACh_Idx),3,'omitnan').*plotBM,cmp=cmpvir,bounds=[0 1],title='ACh vs. Ca^2^+',clabel='r');
% for i = 1:12
%     f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
% end
colorbar off;
title '';
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig2_D1.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

barData = {};
barData{1} = squeeze(mean(subAvg.Beh.R_rfp_HD_gfp_HD(:,:,ACh_Idx).*plotBM,[1,2],'omitnan'));

f = figure(Position=[100,100,200,330]);
[dataMean, dataSEM] = f_plotBar(barData,colors=c_darkCyan,ylabel='r')
ylim([0 1]);

saveas(f, fullfile(fig_savePath, 'ExtDataFig2_D2.svg'));
T = table({order(ACh_Idx).Mouse}',barData{1}, ...
    VariableNames={'Mouse','r'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig2_D2.csv'));

%% Fig Behavior E

f = figure(Position=[100,100,450,450]);
meanSig1 = mean(subAvg.Beh.COH_rfp_HD_gfp_HD(:,ACh_Idx),2);
error1 = std(subAvg.Beh.COH_rfp_HD_gfp_HD(:,ACh_Idx),0,2)/sqrt(M_ACh);
f_plotLineError(fr,meanSig1,error1,color=c_Orange,log=1,lineWidth=3);
ylim([0,1]);
set(gca,'YScale','linear','FontSize',14);
ylabel('Coherence');

yyaxis right;
meanSig2 = mean(subAvg.Beh.PHI_rfp_HD_gfp_HD(:,ACh_Idx),2);
error2 = std(subAvg.Beh.PHI_rfp_HD_gfp_HD(:,ACh_Idx),0,2)/sqrt(M_ACh);
f_plotLineError(fr,meanSig2,error2,color=[0,0,0],log=1,lineWidth=3);
ylim(pi*[-1 1]);
xlim([0.05,5]);
set(gca,'YScale','linear','FontSize',14);
xlabel('F (Hz)');
ylabel('Phi (rad)');

saveas(f, fullfile(fig_savePath, 'ExtDataFig2_E.svg'));
T = table(fr,meanSig1,meanSig2,error1,error2, ...
    VariableNames={'F (Hz)','mean_coherence','mean_phase','SEM_coherence','SEM_phase'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig2_E.csv'));

%% Fig Behavior F

run = 67;

runData = load('sub-Thy1-296_ses-24-04-18_run-01_irun-01_behavior+ophys.mat');

Ca = runData.Fig1.Ca_allen;
ACh = GRAB_FC.GRAB{run};

t = 0.1:0.1:600;

f = figure;hold on;
plot(t,2*Ca(:,5),color=c_Ca);
plot(t,ACh(:,5)-10,color=c_Orange);
plot(t,2*Ca(:,3)-20,color=c_Ca);
plot(t,ACh(:,3)-30,color=c_Orange);
plot([110, 110],[0, 10],'-k');
xlim([100, 220]);

idx = 1000:2200;

saveas(f, fullfile(fig_savePath, 'ExtDataFig2_F.svg'));
T = table(t(idx)',Ca(idx,5),ACh(idx,5),Ca(idx,3),ACh(idx,3), ...
    VariableNames={'Time','Ca_SSpll','ACh_SSpll','Ca_SSpbfd','ACh_SSpbfd'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig2_F.csv'));

%% Fig Behavior G

run = 67;
runData = load('sub-Thy1-296_ses-24-04-18_run-01_irun-01_behavior+ophys.mat');

Ca = runData.Fig1.Ca_allen;
ACh = GRAB_FC.GRAB{run};

lm = fitlm(ACh(:,5),Ca(:,5));
lm = table2array(lm.Coefficients);
f_corr(ACh(:,5),Ca(:,5),1)

f = figure;hold on;
scatter(ACh(:,5),Ca(:,5),75,'filled',MarkerFaceAlpha=0.1,MarkerFaceColor=c_Orange);
plot([-30,30],[-30,30]*lm(2,1)+lm(1,1),'-k');
ylim([-10,15]);
xlim([-30,30]);
saveas(f, fullfile(fig_savePath, 'ExtDataFig2_G1.svg'));
T = table(Ca(idx,5),ACh(idx,5), ...
    VariableNames={'Ca_SSpll','ACh_SSpll'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig2_G1.csv'));

run = 129;
runData = load('sub-Thy1-332_ses-25-02-04_run-01_irun-01_behavior+ophys.mat');

Ca = runData.Fig1.Ca_allen;
NE = GRAB_FC.GRAB{run};

lm = fitlm(NE(:,5),Ca(:,5));
lm = table2array(lm.Coefficients);
f_corr(NE(:,5),Ca(:,5),1)

f = figure;hold on;
scatter(NE(:,5),Ca(:,5),75,'filled',MarkerFaceAlpha=0.1,MarkerFaceColor=c_GRAB);
plot([-6,8],[-6,8]*lm(2,1)+lm(1,1),'-k');
ylim([-10,15]);
xlim([-6,8]);
saveas(f, fullfile(fig_savePath, 'ExtDataFig2_G2.svg'));
T = table(Ca(idx,5),NE(idx,5), ...
    VariableNames={'Ca_SSpll','NE_SSpll'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig2_G2.csv'));

%% Fig Behavior G
% 
% mIdx = 3;
% runIdx = order(mIdx).Runs;
% 
% cmp = cmpinf;
% cmp = cmp(1:end-20,:);
% 
% X = Fig1.GRAB_global(runIdx);
% Y = Fig1.SSp_perf_dt(runIdx);
% f = figure;
% f_multiScatter(X,Y,cmp=cmp,alpha=0.5,lineWidth=2,xlim=[-6 8],ylim=[-1 1],ylabel='r',xlabel='ACh');
% 
% barData = {};
% barData{1} = squeeze(subAvg.Fig1.SSp_perf_vs_GRAB_global(2,ACh_Idx));
% 
% f = figure(Position=[100,100,200,330]);
% [dataMean, dataSEM] = f_plotBar(barData,colors=c_darkCyan,ylabel='r');
% ylim([-1 1]);

%% Fig Behavior H

f = figure;
f_plotMap(mean(subAvg.Beh.R_rfp_HD_gfp_HD(:,:,NE_Idx),3,'omitnan').*plotBM,cmp=cmpvir,bounds=[0 1],title='NE vs. Ca^2^+',clabel='r');
% for i = 1:12
%     f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
% end
colorbar off;
title '';
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig2_H1.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

barData = {};
barData{1} = squeeze(mean(subAvg.Beh.R_rfp_HD_gfp_HD(:,:,NE_Idx).*plotBM,[1,2],'omitnan'));

f = figure(Position=[100,100,200,330]);
[dataMean, dataSEM] = f_plotBar(barData,colors=c_darkCyan,ylabel='r')
ylim([0 1]);

saveas(f, fullfile(fig_savePath, 'ExtDataFig2_H2.svg'));
T = table({order(NE_Idx).Mouse}',barData{1}, ...
    VariableNames={'Mouse','r'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig2_H2.csv'));

%% Fig Behavior I

f = figure(Position=[100,100,450,450]);
meanSig1 = mean(subAvg.Beh.COH_rfp_HD_gfp_HD(:,NE_Idx),2);
error1 = std(subAvg.Beh.COH_rfp_HD_gfp_HD(:,NE_Idx),0,2)/sqrt(M_NE);
f_plotLineError(fr,meanSig1,error1,color=c_GRAB,log=1,lineWidth=3);
ylim([0,1]);
set(gca,'YScale','linear','FontSize',14);
ylabel('Coherence');

yyaxis right;
meanSig2 = mean(subAvg.Beh.PHI_rfp_HD_gfp_HD(:,NE_Idx),2);
error2 = std(subAvg.Beh.PHI_rfp_HD_gfp_HD(:,NE_Idx),0,2)/sqrt(M_NE);
f_plotLineError(fr,meanSig2,error2,color=[0,0,0],log=1,lineWidth=3);
ylim(pi*[-1 1]);
xlim([0.05,5]);
set(gca,'YScale','linear','FontSize',14);
xlabel('F (Hz)');
ylabel('Phi (rad)');

saveas(f, fullfile(fig_savePath, 'ExtDataFig2_I.svg'));
T = table(fr,meanSig1,meanSig2,error1,error2, ...
    VariableNames={'F (Hz)','mean_coherence','mean_phase','SEM_coherence','SEM_phase'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig2_I.csv'));

%% Fig Behavior J

meanSig = mean(subAvg.Beh.NE_IRF(:,NE_Idx),2);
error = std(subAvg.Beh.NE_IRF(:,NE_Idx),0,2)/sqrt(M_NE);
t = -5:0.1:10;

f = figure(Position=[100,100,300,200]);
f_plotLineError(t,meanSig,error,color=c_darkCyan);
xlim([-3 7]);
xlabel('Time (s)');
ylabel('a.u.');
set(gca,'FontSize',14);
title('NE IRF');

saveas(f, fullfile(fig_savePath, 'ExtDataFig2_J1.svg'));
T = table(round(t',2),meanSig,error, ...
    VariableNames={'Time','mean','SEM'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig2_J1.csv'));

f = figure;
f_plotMap(mean(subAvg.Beh.NE_IRF_perf(:,:,NE_Idx),3,'omitnan').*plotBM,cmp=cmpvir,bounds=[0 1],title='NE IRF vs. Ca^2^+',clabel='r');
% for i = 1:12
%     f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
% end
colorbar off;
title '';
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig2_J2.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

barData = {};
barData{1} = squeeze(mean(subAvg.Beh.NE_IRF_perf(:,:,NE_Idx).*plotBM,[1,2],'omitnan'));

f = figure(Position=[100,100,200,330]);
[dataMean, dataSEM] = f_plotBar(barData,colors=c_darkCyan,ylabel='r')
ylim([0 1]);

saveas(f, fullfile(fig_savePath, 'ExtDataFig2_J3.svg'));
T = table({order(NE_Idx).Mouse}',barData{1}, ...
    VariableNames={'Mouse','r'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig2_J3.csv'));


%% Fig Behavior K

f = figure;
f_plotFC(mean(subAvg.Beh.GRAB_conn(:,:,NE_Idx),3),1,cmp=cmpvir,bounds=[0 1],title='NE Connectivity',clabel='r');
colorbar off;
title '';
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig2_K.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);
writetable(table(mean(subAvg.Beh.GRAB_conn(:,:,NE_Idx),3)), fullfile(fig_savePath, 'ExtDataFig2_K.csv'));

%% Fig Behavior L

f = figure(Position=[100,100,450,450]);
meanSig1 = mean(subAvg.Beh.COH_rfp_HD_HbT,2);
error1 = std(subAvg.Beh.COH_rfp_HD_HbT,0,2)/sqrt(M);
f_plotLineError(fr,meanSig1,error1,color=c_Ca,log=1,lineWidth=3);
ylim([0,1]);
set(gca,'YScale','linear','FontSize',14);
ylabel('Coherence');

yyaxis right;
meanSig2 = mean(subAvg.Beh.PHI_rfp_HD_HbT,2);
error2 = std(subAvg.Beh.PHI_rfp_HD_HbT,0,2)/sqrt(M);
f_plotLineError(fr,meanSig2,error2,color=[0,0,0],log=1,lineWidth=3);
ylim(pi*[-1 1]);
xlim([0.05,5]);
set(gca,'YScale','linear','FontSize',14);
xlabel('F (Hz)');
ylabel('Phi (rad)');

saveas(f, fullfile(fig_savePath, 'ExtDataFig2_L1.svg'));
T = table(fr,meanSig1,meanSig2,error1,error2, ...
    VariableNames={'F (Hz)','mean_coherence','mean_phase','SEM_coherence','SEM_phase'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig2_L1.csv'));

f = figure(Position=[100,100,450,450]);
meanSig1 = mean(subAvg.Beh.COH_gfp_HD_HbT(:,NE_Idx),2);
error1 = std(subAvg.Beh.COH_gfp_HD_HbT(:,NE_Idx),0,2)/sqrt(M_NE);
f_plotLineError(fr,meanSig1,error1,color=c_GRAB,log=1,lineWidth=3);
ylim([0,1]);
set(gca,'YScale','linear','FontSize',14);
ylabel('Coherence');

yyaxis right;
meanSig2 = mean(subAvg.Beh.PHI_gfp_HD_HbT(:,NE_Idx),2);
error2 = std(subAvg.Beh.PHI_gfp_HD_HbT(:,NE_Idx),0,2)/sqrt(M_NE);
f_plotLineError(fr,meanSig2,error2,color=[0,0,0],log=1,lineWidth=3);
ylim(pi*[-1 1]);
xlim([0.05,5]);
set(gca,'YScale','linear','FontSize',14);
xlabel('F (Hz)');
ylabel('Phi (rad)');

saveas(f, fullfile(fig_savePath, 'ExtDataFig2_L2.svg'));
T = table(fr,meanSig1,meanSig2,error1,error2, ...
    VariableNames={'F (Hz)','mean_coherence','mean_phase','SEM_coherence','SEM_phase'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig2_L2.csv'));

f = figure(Position=[100,100,450,450]);
meanSig1 = mean(subAvg.Beh.COH_gfp_HD_HbT(:,ACh_Idx),2);
error1 = std(subAvg.Beh.COH_gfp_HD_HbT(:,ACh_Idx),0,2)/sqrt(M_ACh);
f_plotLineError(fr,meanSig1,error1,color=c_Orange,log=1,lineWidth=3);
ylim([0,1]);
set(gca,'YScale','linear','FontSize',14);
ylabel('Coherence');

yyaxis right;
meanSig2 = mean(subAvg.Beh.PHI_gfp_HD_HbT(:,ACh_Idx),2);
error2 = std(subAvg.Beh.PHI_gfp_HD_HbT(:,ACh_Idx),0,2)/sqrt(M_ACh);
f_plotLineError(fr,meanSig2,error2,color=[0,0,0],log=1,lineWidth=3);
ylim(pi*[-1 1]);
xlim([0.05,5]);
set(gca,'YScale','linear','FontSize',14);
xlabel('F (Hz)');
ylabel('Phi (rad)');

saveas(f, fullfile(fig_savePath, 'ExtDataFig2_L3.svg'));
T = table(fr,meanSig1,meanSig2,error1,error2, ...
    VariableNames={'F (Hz)','mean_coherence','mean_phase','SEM_coherence','SEM_phase'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig2_L3.csv'));

%% Fig Behavior M

t = 10:-0.1:-10;

f = figure;
meanSig1 = mean(subAvg.Beh.XC_gfp_HD_HbT(:,NE_Idx),2);
error1 = std(subAvg.Beh.XC_gfp_HD_HbT(:,NE_Idx),0,2)/sqrt(M_NE);
f_plotLineError(t,meanSig1,error1,color=c_GRAB);
meanSig2 = mean(subAvg.Beh.XC_gfp_HD_HbT(:,ACh_Idx),2);
error2 = std(subAvg.Beh.XC_gfp_HD_HbT(:,ACh_Idx),0,2)/sqrt(M_ACh);
f_plotLineError(t,meanSig2,error2,color=c_Orange);
meanSig3 = mean(subAvg.Beh.XC_rfp_HD_HbT,2);
error3 = std(subAvg.Beh.XC_rfp_HD_HbT,0,2)/sqrt(M);
f_plotLineError(t,meanSig3,error3,color=c_Ca);
xlim([-5 5]);
xlabel('Time (s)');
ylabel('r');
set(gca,'FontSize',14);
title('x vs. HbT');
legend('','NE','','ACh','','Ca^2^+');

saveas(f, fullfile(fig_savePath, 'ExtDataFig2_M.svg'));
T = table(t',meanSig1,meanSig2,meanSig3,error1,error2,error3, ...
    VariableNames={'Lag','mean_NE','mean_ACh','mean_Ca', ...
    'SEM_NE','SEM_ACh','SEM_Ca'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig2_M.csv'));

%% Extended data fig 4 E

fig_savePath = fullfile(savePath,'ExtDataFig4');
[~, ~, ~] = mkdir(fig_savePath);

t = 10:-0.1:-10;

f = figure;
meanSig = mean(subAvg.Beh.XC_gfp_HD_HbT_allen(:,2,NE_Idx),3);
error = std(subAvg.Beh.XC_gfp_HD_HbT_allen(:,2,NE_Idx),0,3)/sqrt(M_NE);
f_plotLineError(t,meanSig,error,color=c_GRAB,lineWidth=3);
xlim([-5 5]);
ylim(0.7*[-1 1]);
xlabel('Time (s)');
ylabel('r');
set(gca,'FontSize',14);
title('NE(MOs) vs. HbT');
saveas(f, fullfile(fig_savePath, 'ExtDataFig4_E1.svg'));
T = table(t',meanSig,error, ...
    VariableNames={'Lag','mean_r','SEM_r'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig4_E1.csv'));

f = figure;
meanSig = mean(subAvg.Beh.XC_gfp_HD_HbT_allen(:,5,NE_Idx),3);
error = std(subAvg.Beh.XC_gfp_HD_HbT_allen(:,5,NE_Idx),0,3)/sqrt(M_NE);
f_plotLineError(t,meanSig,error,color=c_GRAB,lineWidth=3);
xlim([-5 5]);
ylim(0.7*[-1 1]);
xlabel('Time (s)');
ylabel('r');
set(gca,'FontSize',14);
title('NE(SSp-ll) vs. HbT');
saveas(f, fullfile(fig_savePath, 'ExtDataFig4_E2.svg'));
T = table(t',meanSig,error, ...
    VariableNames={'Lag','mean_r','SEM_r'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig4_E2.csv'));

f = figure;
meanSig = mean(subAvg.Beh.XC_gfp_HD_HbT_allen(:,12,NE_Idx),3);
error = std(subAvg.Beh.XC_gfp_HD_HbT_allen(:,12,NE_Idx),0,3)/sqrt(M_NE);
f_plotLineError(t,meanSig,error,color=c_GRAB,lineWidth=3);
xlim([-5 5]);
ylim(0.7*[-1 1]);
xlabel('Time (s)');
ylabel('r');
set(gca,'FontSize',14);
title('NE(VISp) vs. HbT');
saveas(f, fullfile(fig_savePath, 'ExtDataFig4_E3.svg'));
T = table(t',meanSig,error, ...
    VariableNames={'Lag','mean_r','SEM_r'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig4_E3.csv'));

%%
function spg = f_adjust_SPG(data,length,norm)
    N = numel(data);
    spg = data;
    for i = 1:N
        if numel(spg{i}) ~= length
            spg{i} = movmean(spg{i},2);
            spg{i} = spg{i}(1:2:end);
        end
        if norm
            spg{i} = spg{i}/mean(spg{i})/5;
        end
    end
end
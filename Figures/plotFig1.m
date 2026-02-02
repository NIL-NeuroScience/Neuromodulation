close all

% calculate subject averages

tmp = struct;
tmp.inv_perf = f_regImages(Fig1.inv_perf,refParcellation,settings.allen_masks,0).*BM;
tmp.SSp_perf = f_regImages(Fig1.SSp_perf,refParcellation,settings.allen_masks,0).*BM;
tmp.var_perf = f_regImages(Fig1.var_perf,refParcellation,settings.allen_masks,0).*BM;
tmp.inv_IRF = cat(2,Fig1.inv_IRF{:});
tmp.SSp_IRF = cat(2,Fig1.SSp_IRF{:});
tmp.var_IRF = cat(3,Fig1.var_IRF{:});
tmp.SSp_perf_vs_GRAB = cat(2,Fig1.SSp_perf_vs_GRAB{:});

subAvg.Fig1.inv_perf = NaN(500,600,M);
subAvg.Fig1.SSp_perf = NaN(500,600,M);
subAvg.Fig1.var_perf = NaN(500,600,M);
subAvg.Fig1.inv_IRF = NaN(101,M);
subAvg.Fig1.SSp_IRF = NaN(101,M);
subAvg.Fig1.var_IRF = NaN(101,12,M);
subAvg.Fig1.SSp_perf_vs_GRAB = NaN(12,M);
subAvg.Fig1.SSp_perf_vs_GRAB_global = NaN(12,M);

for i = 1:M
    subAvg.Fig1.inv_perf(:,:,i) = mean(tmp.inv_perf(:,:,order(i).Runs),3,'omitnan');
    subAvg.Fig1.SSp_perf(:,:,i) = mean(tmp.SSp_perf(:,:,order(i).Runs),3,'omitnan');
    subAvg.Fig1.var_perf(:,:,i) = mean(tmp.var_perf(:,:,order(i).Runs),3,'omitnan');
    subAvg.Fig1.inv_IRF(:,i) = mean(tmp.inv_IRF(:,order(i).Runs),2,'omitnan');
    subAvg.Fig1.SSp_IRF(:,i) = mean(tmp.SSp_IRF(:,order(i).Runs),2,'omitnan');
    subAvg.Fig1.var_IRF(:,:,i) = mean(tmp.var_IRF(:,:,order(i).Runs),3,'omitnan');
    subAvg.Fig1.SSp_perf_vs_GRAB(:,i) = mean(tmp.SSp_perf_vs_GRAB(:,order(i).Runs),2,'omitnan');
    subAvg.Fig1.SSp_perf_vs_GRAB_global(:,i) = mean(cat(1,Fig1.SSp_perf_vs_GRAB_global{order(i).Runs}),1,'omitnan');
end

plotBM = refBM;
plotBM(:,1:300) = NaN;

fig_savePath = fullfile(savePath,'Figure1');
[~, ~, ~] = mkdir(fig_savePath);

%% Fig 1B

run = 129;

Ca = spectra.Ca{run};
HbT = spectra.HbT{run};
NE = GRAB_FC.GRAB{run};
Pupil = Behavior.signals{run}(:,4);
Whisking = Behavior.signals{run}(:,5);
Accelerometer = Behavior.signals{run}(:,6);

t = 0.1:0.1:600;

f = figure;

tiledlayout(5,1);

nexttile;hold on;
plot(t,Ca(:,5),color=c_Ca);
plot(t,NE(:,5),color=c_GRAB);
plot(t,2*HbT(:,5),color=c_HbT);
plot([50,50],[0,10],'-k',LineWidth=2);
axis off;
xlim([50, 350]);
ylim([-10, 20]);

nexttile;hold on;
plot(t,Ca(:,3),color=c_Ca);
plot(t,NE(:,3),color=c_GRAB);
plot(t,2*HbT(:,3),color=c_HbT);
axis off;
xlim([50, 350]);
ylim([-10, 20]);

nexttile;hold on;
plot(t,Pupil,color=c_pupil);
plot([50, 50],[0.2, 0.7],'-k',LineWidth=2);
axis off;
xlim([50, 350]);

nexttile;hold on;
plot(t,Whisking,color=[0,0.7,0.7]);
axis off;
xlim([50, 350]);

nexttile;hold on;
plot(t,Accelerometer,color=[0,0,0]);
xlim([50, 350]);
saveas(f, fullfile(fig_savePath, 'Figure1_B.svg'));

idx = 500:3500;

T = table(t(idx)',Ca(idx,5),NE(idx,5),HbT(idx,5),Ca(idx,3),NE(idx,3), ...
    HbT(idx,3),Pupil(idx),Whisking(idx),Accelerometer(idx), ...
    VariableNames={'Time','Ca_SSpll','NE_SSpll','HbT_SSpll', ...
    'Ca_SSpbfd','NE_SSpbfd','HbT_SSpbfd','Pupil','Whisking','Accelerometer'});
writetable(T, fullfile(fig_savePath, 'Figure1_B.csv'));

%% Fig 1C

f = figure(Position = [100, 100, 500, 400]);
% f_plotMap(mean(subAvg.Fig1.inv_perf,3,'omitnan').*plotBM,cmp=cmpvir,bounds=[0 1],title='Global IRF',clabel='r');
% % for i = 1:12
% %     f_plotAllenRegion(i,2,linewidth=2,color=[1 0 0]);
% % end
% set(gca,'YDir','reverse');
imagesc(mean(subAvg.Fig1.inv_perf,3,'omitnan').*plotBM,AlphaData=plotBM);
colormap cmpvir;
clim([0, 1]);
axis image off;
exportgraphics(f, fullfile(fig_savePath,'Figure1_C1.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

meanSig = mean(subAvg.Fig1.inv_IRF,2);
error = std(subAvg.Fig1.inv_IRF,0,2)/sqrt(M);

f = figure(Position=[100,100,300,200]);
f_plotLineError(0:0.1:10,meanSig,error,color=c_darkCyan);
xlim([0 7]);
xlabel('Time (s)');
ylabel('a.u.');
set(gca,'FontSize',14);
title('Global IRF');

T = table((0:0.1:7)', meanSig(1:71), error(1:71), ...
    VariableNames={'Time','mean','SEM'});
writetable(T, fullfile(fig_savePath, 'Figure1_C2.csv'));
saveas(f, fullfile(fig_savePath, 'Figure1_C2.svg'));

%% Fig 1D

f = figure(Position = [100, 100, 500, 400]);
% f_plotMap(mean(subAvg.Fig1.SSp_perf,3,'omitnan').*plotBM,cmp=cmpvir,bounds=[0 1],title='Global IRF',clabel='r');
% for i = 1:12
%     f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
% end
% f_plotAllenRegion(4,2,linewidth=3,color=[1 0 0]);
% f_plotAllenRegion(5,2,linewidth=3,color=[1 0 0]);
imagesc(mean(subAvg.Fig1.SSp_perf,3,'omitnan').*plotBM,AlphaData=plotBM);
colormap cmpvir;
clim([0, 1]);
axis image off;
exportgraphics(f, fullfile(fig_savePath,'Figure1_D1.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

meanSig = mean(subAvg.Fig1.SSp_IRF,2);
error = std(subAvg.Fig1.SSp_IRF,0,2)/sqrt(M);

f = figure(Position=[100,100,300,200]);
f_plotLineError(0:0.1:10,meanSig,error,color=c_darkCyan);
xlim([0 7]);
xlabel('Time (s)');
ylabel('a.u.');
set(gca,'FontSize',14);
title('SSp IRF');

T = table((0:0.1:7)', meanSig(1:71), error(1:71), ...
    VariableNames={'Time','mean','SEM'});
writetable(T, fullfile(fig_savePath, 'Figure1_D2.csv'));
saveas(f, fullfile(fig_savePath, 'Figure1_D2.svg'));

%% Fig 1E

f = figure(Position = [100, 100, 500, 400]);
% f_plotMap(mean(subAvg.Fig1.var_perf,3,'omitnan').*plotBM,cmp=cmpvir,bounds=[0 1],title='Global IRF',clabel='r');
% for i = 1:12
%     f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
% end
% f_plotAllenRegion(2,2,linewidth=3,color=c_darkCyan);
% f_plotAllenRegion(4,2,linewidth=3,color=c_Orange);
% f_plotAllenRegion(5,2,linewidth=3,color=c_Orange);
imagesc(mean(subAvg.Fig1.var_perf,3,'omitnan').*plotBM,AlphaData=plotBM);
colormap cmpvir;
clim([0, 1]);
axis image off;
exportgraphics(f, fullfile(fig_savePath,'Figure1_E1.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

f = figure(Position=[100,100,300,200]);
meanSig1 = mean(subAvg.Fig1.var_IRF(:,2),3);
error1 = std(subAvg.Fig1.var_IRF(:,2,:),0,3)/sqrt(M);
f_plotLineError(0:0.1:10,meanSig1,error1,color=c_darkCyan);
meanSig2 = mean(subAvg.Fig1.var_IRF(:,[4,5],:),2);
error2 = std(meanSig2,0,3)/sqrt(M);
meanSig2 = mean(meanSig2,3);
f_plotLineError(0:0.1:10,meanSig2,error2,color=c_Orange);
xlim([0 7]);
xlabel('Time (s)');
ylabel('a.u.');
legend('','MOs','','SSp-tr/ll');
title('Variant IRF');
set(gca,'FontSize',14);

T = table((0:0.1:7)', meanSig1(1:71), meanSig1(1:71), error1(1:71), error1(1:71), ...
    VariableNames={'Time','mean_MOs','mean_SSp_tr_ll','SEM_MOs','SEM_SSp_tr_ll'});
writetable(T, fullfile(fig_savePath, 'Figure1_E2.csv'));
saveas(f, fullfile(fig_savePath, 'Figure1_E2.svg'));

%% plot 1F

barData = {};
barData{1} = squeeze(mean(subAvg.Fig1.inv_perf.*plotBM,[1,2],'omitnan'));
barData{2} = squeeze(mean(subAvg.Fig1.SSp_perf.*plotBM,[1,2],'omitnan'));
barData{3} = squeeze(mean(subAvg.Fig1.var_perf.*plotBM,[1,2],'omitnan'));

f = figure;
[dataMean, dataSEM] = f_plotBar(barData,colors=c_darkCyan,ylabel='r');
saveas(f, fullfile(fig_savePath, 'Figure1_F.svg'));
T = table({order.Mouse}', barData{1}, barData{2}, barData{3}, ...
    VariableNames={'Mouse','global','SSp','variant'});
writetable(T, fullfile(fig_savePath, 'Figure1_F.csv'));

[h, p] = f_SRtest(barData,0.05);
T = table(p(:,1),p(:,2),p(:,3), ...
    VariableNames={'global','SSp','variant'});
writetable(T, fullfile(fig_savePath, 'Figure1_F_p.csv'));
% fprintf('%.2e', p);

%% plot 1G

mIdx = 15;
runIdx = order(mIdx).Runs;

cmp = cmpinf;
cmp = cmp(1:end-20,:);

X = Fig1.GRAB_global(runIdx);
Y = Fig1.SSp_perf_dt(runIdx);
f = figure;
f_multiScatter(X,Y,cmp=cmp,alpha=0.5,lineWidth=2,xlim=[-6 8],ylim=[-1 1],ylabel='r',xlabel='NE');
saveas(f, fullfile(fig_savePath, 'Figure1_G1.svg'));
T = table();
for i = 1 : numel(X)
    T.(sprintf('Run%02i_NE',i)) = X{i};
    T.(sprintf('Run%02i_r',i)) = Y{i};
end
writetable(T, fullfile(fig_savePath, 'Figure1_G1.csv'));

barData = {};
barData{1} = squeeze(subAvg.Fig1.SSp_perf_vs_GRAB_global(2,NE_Idx));

f = figure(Position=[100,100,200,330]);
[dataMean, dataSEM] = f_plotBar(barData,colors=c_darkCyan,ylabel='r')
ylim(0.5*[-1 1]);
saveas(f, fullfile(fig_savePath, 'Figure1_G2.svg'));
T = table({order(NE_Idx).Mouse}', barData{1}', ...
    VariableNames={'Mouse','r'});
writetable(T, fullfile(fig_savePath, 'Figure1_G2.csv'));

%% plot 1H

f = figure(Position = [100, 100, 500, 400]);
f_plotAllenMap(mean(subAvg.Fig1.SSp_perf_vs_GRAB_global(:,NE_Idx),2),cmp=cmpbbr,cLabel='r',mask=plotBM,clim=[-0.3, 0.3]);
colorbar off;
exportgraphics(f, fullfile(fig_savePath,'Figure1_H.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

%% plot 1J

f = figure(Position = [100, 100, 500, 400]);
f_plotAllenMap(std(subAvg.Fig1.SSp_perf_vs_GRAB_global(:,NE_Idx),0,2)/sqrt(sum(NE_Idx)),cmp=cmpinf,cLabel='SEM',mask=plotBM,clim=[0, 0.1]);
colorbar off;
exportgraphics(f, fullfile(fig_savePath,'Figure1_J.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);
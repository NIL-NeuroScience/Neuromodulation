close all

% calculate subject averages

NE_order = order(NE_Idx);

subAvg.Fig3.lowNE_FC_Ca = NaN(12,12,numel(NE_order));
subAvg.Fig3.lowNE_FC_HbT = NaN(12,12,numel(NE_order));
subAvg.Fig3.highNE_FC_Ca = NaN(12,12,numel(NE_order));
subAvg.Fig3.highNE_FC_HbT = NaN(12,12,numel(NE_order));
subAvg.Fig3.FC_Ca_HbT_vs_GRAB = NaN(numel(NE_order),1);
subAvg.Fig3.FC_Ca_vs_GRAB = NaN(12,12,numel(NE_order));
subAvg.Fig3.FC_HbT_vs_GRAB = NaN(12,12,numel(NE_order));

for i = 1:numel(NE_order)
    subAvg.Fig3.lowNE_FC_Ca(:,:,i) = mean(cat(3,Fig3.lowNE_Ca{NE_order(i).Runs}),3,'omitnan');
    subAvg.Fig3.lowNE_FC_HbT(:,:,i) = mean(cat(3,Fig3.lowNE_HbT{NE_order(i).Runs}),3,'omitnan');
    subAvg.Fig3.highNE_FC_Ca(:,:,i) = mean(cat(3,Fig3.highNE_Ca{NE_order(i).Runs}),3,'omitnan');
    subAvg.Fig3.highNE_FC_HbT(:,:,i) = mean(cat(3,Fig3.highNE_HbT{NE_order(i).Runs}),3,'omitnan');
    subAvg.Fig3.FC_Ca_HbT_vs_GRAB(i) = mean(cat(1,Fig3.FC_Ca_HbT_vs_GRAB{NE_order(i).Runs}));
    subAvg.Fig3.FC_Ca_vs_GRAB(:,:,i) = mean(cat(3,Fig3.FC_Ca_vs_GRAB{NE_order(i).Runs}),3);
    subAvg.Fig3.FC_HbT_vs_GRAB(:,:,i) = mean(cat(3,Fig3.FC_HbT_vs_GRAB{NE_order(i).Runs}),3);
end

plotBM = refBM;
plotBM(:,1:300) = NaN;

fig_savePath = fullfile(savePath,'Figure3');
[~, ~, ~] = mkdir(fig_savePath);

%% Fig 3 A

run = 59;

Ca = spectra.Ca{run};
HbT = spectra.HbT{run};
NE = GRAB_FC.GRAB{run};
Pupil = Behavior.signals{run}(:,4);
Whisking = Behavior.signals{run}(:,5);
Accelerometer = Behavior.signals{run}(:,6);

t = 0.1:0.1:600;

f = figure;
tiledlayout(6,1);

nexttile;hold on;
plot(t,HbT(:,2),color=0.8*c_HbT);
plot(t,HbT(:,5)-5,color=1.2*c_HbT);
plot([50, 50],[0, 10],'-k',LineWidth=2);
xlim([50,550]);
axis off;

nexttile;hold on;
plot(t,Ca(:,2),color=0.8*c_Ca);
plot(t,Ca(:,5)-5,color=1.2*c_Ca);
plot([50, 50],[0, 10],'-k',LineWidth=2);
xlim([50,550]);
axis off;

nexttile;hold on;
plot(t,NE(:,2),color=0.8*c_GRAB);
plot(t,NE(:,5)-5,color=1.2*c_GRAB);
plot([50, 50],[0, 10],'-k',LineWidth=2);
xlim([50,550]);
axis off;

nexttile;hold on;
plot(t,Pupil,color=c_pupil);
plot([50, 50],[0.2, 0.5],'-k',LineWidth=2);
axis off;
xlim([50, 550]);

nexttile;hold on;
plot(t,Whisking,color=[0,0.7,0.7]);
axis off;
xlim([50, 550]);

nexttile;hold on;
plot(t,Accelerometer,color=[0,0,0]);
xlim([50, 550]);
saveas(f, fullfile(fig_savePath, 'Figure3_A.svg'));

idx = 500:5500;

T = table(t(idx)',Ca(idx,5),NE(idx,5),HbT(idx,5),Ca(idx,2),NE(idx,2), ...
    HbT(idx,2),Pupil(idx),Whisking(idx),Accelerometer(idx), ...
    VariableNames={'Time','Ca_SSpll','NE_SSpll','HbT_SSpll', ...
    'Ca_MOs','NE_MOs','HbT_MOs','Pupil','Whisking','Accelerometer'});
writetable(T, fullfile(fig_savePath, 'Figure3_A.csv'));

%% Fig 3 B

Ca_FC = Fig3.FC_Ca_MOs_SSpll{run};
HbT_FC = Fig3.FC_HbT_MOs_SSpll{run};

t = 15:6:585;

f = figure;hold on;
plot(t,Ca_FC,color=c_Ca);
plot([55,55],[0,1],'-k');
plot(t,HbT_FC-1.5,color=c_HbT);
plot([55,55],[-1.5,-0.5],'-k');
xlim([50,550]);
saveas(f, fullfile(fig_savePath, 'Figure3_B.svg'));

T = table(t', Ca_FC, HbT_FC, ...
    VariableNames={'Time','FC_Ca','FC_HbT'});
writetable(T, fullfile(fig_savePath, 'Figure3_B.csv'));

%% Fig 3 C

FC = Fig3.FC_Ca_vs_HbT{run};

t = 15:6:585;

f = figure;hold on;
plot(t,FC,color=[0,0.7,0.7]);
plot([55,55],[0,1],'-k');
xlim([50,550]);
saveas(f, fullfile(fig_savePath, 'Figure3_C.svg'));

T = table(t', FC, ...
    VariableNames={'Time','FC_Ca_vs_HbT'});
writetable(T, fullfile(fig_savePath, 'Figure3_C.csv'));

%% Fig 3 D

mIdx = 8;
runIdx = order(mIdx).Runs;

cmp = cmpinf;
cmp = cmp(1:220,:);

X = Fig3.NE(runIdx);
Y = Fig3.FC_Ca_MOs_SSpll(runIdx);
f = figure;
f_multiScatter(X,Y,cmp=cmp,alpha=0.5,lineWidth=2,xlim=[-4 5],ylim=[0 1],ylabel='r',xlabel='NE');
saveas(f, fullfile(fig_savePath, 'Figure3_D1.svg'));

T = table();
for i = 1 : numel(X)
    T.(sprintf('Run%02i_NE',i)) = X{i};
    T.(sprintf('Run%02i_FC_Ca',i)) = Y{i};
end
writetable(T, fullfile(fig_savePath, 'Figure3_D1.csv'));

Y = Fig3.FC_HbT_MOs_SSpll(runIdx);
f = figure;
f_multiScatter(X,Y,cmp=cmp,alpha=0.5,lineWidth=2,xlim=[-4 5],ylim=[0 1],ylabel='r',xlabel='NE');
saveas(f, fullfile(fig_savePath, 'Figure3_D2.svg'));

T = table();
for i = 1 : numel(X)
    T.(sprintf('Run%02i_NE',i)) = X{i};
    T.(sprintf('Run%02i_FC_HbT',i)) = Y{i};
end
writetable(T, fullfile(fig_savePath, 'Figure3_D2.csv'));

X = cat(1,X{:});
f = figure;
h = histfit(X,30,'kernel');
h(1).FaceAlpha = 0.2;
h(1).FaceColor = [0 0 0];
h(1).EdgeColor = 'none';
h(2).Color = [0 0 0];
h(2).LineWidth = 1;
xlim([-4 5]);
axis off;
saveas(f, fullfile(fig_savePath, 'Figure3_D3.svg'));

%% Fig 3 E
X = Fig3.NE(runIdx);
Y = Fig3.FC_Ca_vs_HbT(runIdx);
f = figure;
f_multiScatter(X,Y,cmp=cmp,alpha=0.5,lineWidth=2,xlim=[-4 5],ylim=[0 1],ylabel='r',xlabel='NE');
saveas(f, fullfile(fig_savePath, 'Figure3_E.svg'));

T = table();
for i = 1 : numel(X)
    T.(sprintf('Run%02i_NE',i)) = X{i};
    T.(sprintf('Run%02i_FC',i)) = Y{i};
end
writetable(T, fullfile(fig_savePath, 'Figure3_E.csv'));


%% Fig 3 F

barData = {};
barData{1} = squeeze(subAvg.Fig3.FC_Ca_vs_GRAB(2,5,:));
barData{2} = squeeze(subAvg.Fig3.FC_HbT_vs_GRAB(2,5,:));
barData{3} = subAvg.Fig3.FC_Ca_HbT_vs_GRAB;

f = figure;
[dataMean, dataSEM] = f_plotBar(barData,colors=[c_Ca;c_HbT;c_darkCyan],legend={'Ca^2^+','HbT','Ca^2^+ vs. HbT'},ylabel='r')
saveas(f, fullfile(fig_savePath, 'Figure3_F.svg'));

T = table({order(NE_Idx).Mouse}',barData{1},barData{2},barData{3}, ...
    VariableNames={'Mouse','FC_Ca','FC_HbT','FC_Ca_vs_HbT'});
writetable(T, fullfile(fig_savePath, 'Figure3_F.csv'));

%% Fig 3 G

f = figure;
f_plotFC(mean(subAvg.Fig3.lowNE_FC_Ca,3),1,cmp=cmpvir,clim=[0 1],title='Low NE Ca++ Connectivity',clabel='r');
colorbar off;
title '';
exportgraphics(f, fullfile(fig_savePath,'Figure3_G1.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);
writetable(table(mean(subAvg.Fig3.lowNE_FC_Ca,3)), fullfile(fig_savePath, 'Figure3_G1.csv'));

f = figure;
f_plotFC(mean(subAvg.Fig3.highNE_FC_Ca,3),1,cmp=cmpvir,clim=[0 1],title='High NE Ca++ Connectivity',clabel='r');
colorbar off;
title '';
exportgraphics(f, fullfile(fig_savePath,'Figure3_G2.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);
writetable(table(mean(subAvg.Fig3.highNE_FC_Ca,3)), fullfile(fig_savePath, 'Figure3_G2.csv'));

MouseID = {log(strcmp({log.GRAB},'GRAB_NE')).Mouse};
g1 = cat(3,Fig3.lowNE_Ca{strcmp({log.GRAB},'GRAB_NE')});
g1 = reshape(g1,144,[])';
g2 = cat(3,Fig3.highNE_Ca{strcmp({log.GRAB},'GRAB_NE')});
g2 = reshape(g2,144,[])';
[h,p] = f_lme(MouseID,g1,g2,0.05);

idx = logical(eye(12));
h(idx(:)) = 0;

f = figure;
f_plotFC(mean(subAvg.Fig3.highNE_FC_Ca,3)-mean(subAvg.Fig3.lowNE_FC_Ca,3),0,cmp=cmpbbr,clim=0.25*[-1 1],title='High-Low NE',clabel='\Deltar');
f_overlayStats_FC(reshape(h,12,12));
saveas(f, fullfile(fig_savePath, 'Figure3_G3.svg'));
writetable(table(mean(subAvg.Fig3.highNE_FC_Ca,3)-mean(subAvg.Fig3.lowNE_FC_Ca,3)), fullfile(fig_savePath, 'Figure3_G3.csv'));
writetable(table(reshape(p,12,12)), fullfile(fig_savePath, 'Figure3_G3_P.csv'));

f = figure;
f_plotFC(mean(subAvg.Fig3.highNE_FC_Ca,3)-mean(subAvg.Fig3.lowNE_FC_Ca,3),0,cmp=cmpbbr,clim=0.25*[-1 1],title='High-Low NE',clabel='\Deltar');
colorbar off;
title '';
exportgraphics(f, fullfile(fig_savePath,'Figure3_G3.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

%% Fig 3 H

f = figure;
f_plotFC(mean(subAvg.Fig3.lowNE_FC_HbT,3),1,cmp=cmpvir,clim=[0 1],title='Low NE HbT Connectivity',clabel='r');
colorbar off;
title '';
exportgraphics(f, fullfile(fig_savePath,'Figure3_H1.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);
writetable(table(mean(subAvg.Fig3.lowNE_FC_Ca,3)), fullfile(fig_savePath, 'Figure3_H1.csv'));

f = figure;
f_plotFC(mean(subAvg.Fig3.highNE_FC_HbT,3),1,cmp=cmpvir,clim=[0 1],title='High NE HbT Connectivity',clabel='r');
colorbar off;
title '';
exportgraphics(f, fullfile(fig_savePath,'Figure3_H2.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);
writetable(table(mean(subAvg.Fig3.highNE_FC_Ca,3)), fullfile(fig_savePath, 'Figure3_H2.csv'));

MouseID = {log(strcmp({log.GRAB},'GRAB_NE')).Mouse};
g1 = cat(3,Fig3.lowNE_HbT{strcmp({log.GRAB},'GRAB_NE')});
g1 = reshape(g1,144,[])';
g2 = cat(3,Fig3.highNE_HbT{strcmp({log.GRAB},'GRAB_NE')});
g2 = reshape(g2,144,[])';
[h,p] = f_lme(MouseID,g1,g2,0.05);

idx = logical(eye(12));
h(idx(:)) = 0;

f = figure;
f_plotFC(mean(subAvg.Fig3.highNE_FC_HbT,3)-mean(subAvg.Fig3.lowNE_FC_HbT,3),0,cmp=cmpbbr,clim=0.1*[-1 1],title='High-Low NE',clabel='\Deltar');
f_overlayStats_FC(reshape(h,12,12));
saveas(f, fullfile(fig_savePath, 'Figure3_H3.svg'));
writetable(table(mean(subAvg.Fig3.highNE_FC_Ca,3)-mean(subAvg.Fig3.lowNE_FC_Ca,3)), fullfile(fig_savePath, 'Figure3_H3.csv'));
writetable(table(reshape(p,12,12)), fullfile(fig_savePath, 'Figure3_H3_P.csv'));

f = figure;
f_plotFC(mean(subAvg.Fig3.highNE_FC_HbT,3)-mean(subAvg.Fig3.lowNE_FC_HbT,3),0,cmp=cmpbbr,clim=0.1*[-1 1],title='High-Low NE',clabel='\Deltar');
colorbar off;
title '';
exportgraphics(f, fullfile(fig_savePath,'Figure3_H3.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

%% Fig 3 I

regions = [2,5,12];

allenCaLow = mean(subAvg.Fig3.lowNE_FC_Ca,3);
allenCaLow = allenCaLow(regions,:);
allenCaHigh = mean(subAvg.Fig3.highNE_FC_Ca,3);
allenCaHigh = allenCaHigh(regions,:);

for i = 1:3
    f = figure;f_plotAllenMap(allenCaLow(i,:),cmp=cmpvir,mask=plotBM,clim=[0,1]);
    f_plotAllenRegion(regions(i),2,linewidth=3,color=[0 0 0]);
    colorbar off;
    exportgraphics(f, fullfile(fig_savePath,sprintf('Figure3_I%01i_low.jpg',i)),'Resolution',300,'BackgroundColor',[1 1 1]);
    
    f = figure;f_plotAllenMap(allenCaHigh(i,:),cmp=cmpvir,mask=plotBM,clim=[0,1]);
    f_plotAllenRegion(regions(i),2,linewidth=3,color=[0 0 0]);
    colorbar off;
    exportgraphics(f, fullfile(fig_savePath,sprintf('Figure3_I%01i_high.jpg',i)),'Resolution',300,'BackgroundColor',[1 1 1]);
    
    f = figure;f_plotAllenMap(allenCaHigh(i,:)-allenCaLow(i,:),cmp=cmpbbr,mask=plotBM,clim=0.25*[-1,1]);
    f_plotAllenRegion(regions(i),2,linewidth=3,color=0.7*[1 1 1]);
    colorbar off;
    exportgraphics(f, fullfile(fig_savePath,sprintf('Figure3_I%01i_diff.jpg',i)),'Resolution',300,'BackgroundColor',[1 1 1]);
end

%% Fig 3 J

regions = [2,5,12];

allenHbTLow = mean(subAvg.Fig3.lowNE_FC_HbT,3);
allenHbTLow = allenHbTLow(regions,:);
allenHbTHigh = mean(subAvg.Fig3.highNE_FC_HbT,3);
allenHbTHigh = allenHbTHigh(regions,:);

for i = 1:3
    f = figure;f_plotAllenMap(allenHbTLow(i,:),cmp=cmpvir,mask=plotBM,clim=[0,1]);
    f_plotAllenRegion(regions(i),2,linewidth=3,color=[0 0 0]);
    colorbar off;
    exportgraphics(f, fullfile(fig_savePath,sprintf('Figure3_J%01i_low.jpg',i)),'Resolution',300,'BackgroundColor',[1 1 1]);

    f = figure;f_plotAllenMap(allenHbTHigh(i,:),cmp=cmpvir,mask=plotBM,clim=[0,1]);
    f_plotAllenRegion(regions(i),2,linewidth=3,color=[0 0 0]);
    colorbar off;
    exportgraphics(f, fullfile(fig_savePath,sprintf('Figure3_J%01i_high.jpg',i)),'Resolution',300,'BackgroundColor',[1 1 1]);

    f = figure;f_plotAllenMap(allenHbTHigh(i,:)-allenHbTLow(i,:),cmp=cmpbbr,mask=plotBM,clim=0.1*[-1,1]);
    f_plotAllenRegion(regions(i),2,linewidth=3,color=0.7*[1 1 1]);
    colorbar off;
    exportgraphics(f, fullfile(fig_savePath,sprintf('Figure3_J%01i_diff.jpg',i)),'Resolution',300,'BackgroundColor',[1 1 1]);
end

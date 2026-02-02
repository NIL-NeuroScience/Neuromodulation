close all

% calculate subject averages

NE_order = order(NE_Idx);

tmp = struct;
tmp.inv_perf = f_regImages(Fig1.inv_perf,refParcellation,settings.allen_masks,0).*BM;
tmp.SSp_perf = f_regImages(Fig1.SSp_perf,refParcellation,settings.allen_masks,0).*BM;
tmp.var_perf = f_regImages(Fig1.var_perf,refParcellation,settings.allen_masks,0).*BM;
tmp.LR_perf = f_regImages(Fig2.LR_perf,refParcellation,settings.allen_masks,0).*BM;
tmp.IRFx2_perf = f_regImages(Fig2.IRFx2_perf,refParcellation,settings.allen_masks,0).*BM;
tmp.LR_A = f_regImages(Fig2.LR_A,refParcellation,settings.allen_masks,0).*BM;
tmp.LR_B = f_regImages(Fig2.LR_B,refParcellation,settings.allen_masks,0).*BM;
tmp.IRFx2_A = f_regImages(Fig2.IRFx2_A,refParcellation,settings.allen_masks,0).*BM;
tmp.IRFx2_B = f_regImages(Fig2.IRFx2_B,refParcellation,settings.allen_masks,0).*BM;

subAvg.FigE1.inv_perf = NaN(500,600,numel(NE_order));
subAvg.FigE1.SSp_perf = NaN(500,600,numel(NE_order));
subAvg.FigE1.var_perf = NaN(500,600,numel(NE_order));
subAvg.FigE1.LR_perf = NaN(500,600,numel(NE_order));
subAvg.FigE1.IRFx2_perf = NaN(500,600,numel(NE_order));
subAvg.FigE1.LR_A = NaN(500,600,numel(NE_order));
subAvg.FigE1.LR_B = NaN(500,600,numel(NE_order));
subAvg.FigE1.IRFx2_A = NaN(500,600,numel(NE_order));
subAvg.FigE1.IRFx2_B = NaN(500,600,numel(NE_order));

for i = 1:numel(NE_order)
    subAvg.FigE1.inv_perf(:,:,i) = mean(tmp.inv_perf(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.FigE1.SSp_perf(:,:,i) = mean(tmp.SSp_perf(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.FigE1.var_perf(:,:,i) = mean(tmp.var_perf(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.FigE1.LR_perf(:,:,i) = mean(tmp.LR_perf(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.FigE1.IRFx2_perf(:,:,i) = mean(tmp.IRFx2_perf(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.FigE1.LR_A(:,:,i) = mean(tmp.LR_A(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.FigE1.LR_B(:,:,i) = mean(tmp.LR_B(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.FigE1.IRFx2_A(:,:,i) = mean(tmp.IRFx2_A(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.FigE1.IRFx2_B(:,:,i) = mean(tmp.IRFx2_B(:,:,NE_order(i).Runs),3,'omitnan');
end

plotBM = refBM;
plotBM(:,1:300) = NaN;

BM_outline = plotBM;
BM_outline(isnan(BM_outline)) = 0;
BM_outline = bwboundaries(BM_outline);
BM_outline = BM_outline{1};

plotOrder = [1 2 7 8 3 4 5 6];
Mice = {order(NE_Idx).Mouse};

fig_savePath = fullfile(savePath,'ExtDataFig1');
[~, ~, ~] = mkdir(fig_savePath);

%% Fig AllMice A

for i = 1:8
    f = figure;
    f_plotMap(mean(subAvg.FigE1.inv_perf(:,:,plotOrder(i)),3,'omitnan'),cmp=cmpvir,clim=[0 1],title='LR Performance',clabel='r');
    plot(BM_outline(:,2),BM_outline(:,1),'-k',lineWidth=3);
    set(gca,'YDir','reverse');
    title('');colorbar off;
    filename = fullfile(fig_savePath,['ExtDataFig1_A_' Mice{plotOrder(i)} '.jpg']);
    exportgraphics(f,filename,'Resolution',300,'BackgroundColor',[1 1 1]);
end

%% Fig AllMice B

for i = 1:8
    f = figure;
    f_plotMap(mean(subAvg.FigE1.SSp_perf(:,:,plotOrder(i)),3,'omitnan'),cmp=cmpvir,clim=[0 1],title='LR Performance',clabel='r');
    plot(BM_outline(:,2),BM_outline(:,1),'-k',lineWidth=3);
    set(gca,'YDir','reverse');
    title('');colorbar off;
    filename = fullfile(fig_savePath,['ExtDataFig1_B_' Mice{plotOrder(i)} '.jpg']);
    exportgraphics(f,filename,'Resolution',300,'BackgroundColor',[1 1 1]);
end

%% Fig AllMice C

for i = 1:8
    f = figure;
    f_plotMap(mean(subAvg.FigE1.var_perf(:,:,plotOrder(i)),3,'omitnan'),cmp=cmpvir,clim=[0 1],title='LR Performance',clabel='r');
    plot(BM_outline(:,2),BM_outline(:,1),'-k',lineWidth=3);
    set(gca,'YDir','reverse');
    title('');colorbar off;
    filename = fullfile(fig_savePath,['ExtDataFig1_C_' Mice{plotOrder(i)} '.jpg']);
    exportgraphics(f,filename,'Resolution',300,'BackgroundColor',[1 1 1]);
end

%% Fig AllMice D

for i = 1:8
    f = figure;
    f_plotMap(mean(subAvg.FigE1.LR_perf(:,:,plotOrder(i)),3,'omitnan'),cmp=cmpvir,clim=[0 1],title='LR Performance',clabel='r');
    plot(BM_outline(:,2),BM_outline(:,1),'-k',lineWidth=3);
    set(gca,'YDir','reverse');
    title('');colorbar off;
    filename = fullfile(fig_savePath,['ExtDataFig1_D_' Mice{plotOrder(i)} '.jpg']);
    exportgraphics(f,filename,'Resolution',300,'BackgroundColor',[1 1 1]);
end

%% Fig AllMice E

for i = 1:8
    f = figure;
    f_plotMap(mean(subAvg.FigE1.IRFx2_perf(:,:,plotOrder(i)),3,'omitnan'),cmp=cmpvir,clim=[0 1],title='LR Performance',clabel='r');
    plot(BM_outline(:,2),BM_outline(:,1),'-k',lineWidth=3);
    set(gca,'YDir','reverse');
    title('');colorbar off;
    filename = fullfile(fig_savePath,['ExtDataFig1_E_' Mice{plotOrder(i)} '.jpg']);
    exportgraphics(f,filename,'Resolution',300,'BackgroundColor',[1 1 1]);
end

%% Fig AllMice F

for i = 1:8
    f = figure;
    f_plotMap(mean(subAvg.FigE1.LR_A(:,:,plotOrder(i)),3,'omitnan'),cmp=cmpbbr,clim=[-1 1],title='LR Performance',clabel='r');
    plot(BM_outline(:,2),BM_outline(:,1),'-k',lineWidth=3,color=0.7*[1,1,1]);
    set(gca,'YDir','reverse');
    title('');colorbar off;
    filename = fullfile(fig_savePath,['ExtDataFig1_F_' Mice{plotOrder(i)} '.jpg']);
    exportgraphics(f,filename,'Resolution',300,'BackgroundColor',[1 1 1]);
end

%% Fig AllMice G

for i = 1:8
    f = figure;
    f_plotMap(mean(subAvg.FigE1.LR_B(:,:,plotOrder(i)),3,'omitnan'),cmp=cmpbbr,clim=[-1 1],title='LR Performance',clabel='r');
    plot(BM_outline(:,2),BM_outline(:,1),'-k',lineWidth=3,color=0.7*[1,1,1]);
    set(gca,'YDir','reverse');
    title('');colorbar off;
    filename = fullfile(fig_savePath,['ExtDataFig1_G_' Mice{plotOrder(i)} '.jpg']);
    exportgraphics(f,filename,'Resolution',300,'BackgroundColor',[1 1 1]);
end


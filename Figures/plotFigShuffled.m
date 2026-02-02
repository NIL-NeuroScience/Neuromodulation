% calculate subject averages

NE_order = order(NE_Idx);

tmp = struct;
tmp.Unf.LR_perf = f_regImages(unfiltered.LR_perf,refParcellation,settings.allen_masks,0).*BM;
tmp.Unf.IRFx2_perf = f_regImages(unfiltered.IRFx2_perf,refParcellation,settings.allen_masks,0).*BM;
tmp.Shu.LR_perf = f_regImages(shuffled.LR_perf,refParcellation,settings.allen_masks,0).*BM;
tmp.Shu.IRFx2_perf = f_regImages(shuffled.IRFx2_perf,refParcellation,settings.allen_masks,0).*BM;
tmp.Unf.LR_A = f_regImages(unfiltered.LR_A,refParcellation,settings.allen_masks,0).*BM;
tmp.Unf.LR_B = f_regImages(unfiltered.LR_B,refParcellation,settings.allen_masks,0).*BM;
tmp.Unf.IRFx2_A = f_regImages(unfiltered.IRFx2_A,refParcellation,settings.allen_masks,0).*BM;
tmp.Unf.IRFx2_B = f_regImages(unfiltered.IRFx2_B,refParcellation,settings.allen_masks,0).*BM;
tmp.Shu.LR_A = f_regImages(shuffled.LR_A,refParcellation,settings.allen_masks,0).*BM;
tmp.Shu.LR_B = f_regImages(shuffled.LR_B,refParcellation,settings.allen_masks,0).*BM;
tmp.Shu.IRFx2_A = f_regImages(shuffled.IRFx2_A,refParcellation,settings.allen_masks,0).*BM;
tmp.Shu.IRFx2_B = f_regImages(shuffled.IRFx2_B,refParcellation,settings.allen_masks,0).*BM;

subAvg.Unf.LR_perf = NaN(500,600,numel(NE_order));
subAvg.Unf.IRFx2_perf = NaN(500,600,numel(NE_order));
subAvg.Fig2.LR_perf = NaN(500,600,numel(NE_order));
subAvg.Fig2.IRFx2_perf = NaN(500,600,numel(NE_order));
subAvg.Unf.LR_A = NaN(500,600,numel(NE_order));
subAvg.Unf.LR_B = NaN(500,600,numel(NE_order));
subAvg.Unf.IRFx2_A = NaN(500,600,numel(NE_order));
subAvg.Unf.IRFx2_B = NaN(500,600,numel(NE_order));
subAvg.Unf.tA = NaN(numel(NE_order),1);
subAvg.Unf.tB = NaN(numel(NE_order),1);
subAvg.Unf.IRFx2_IRF = NaN(151,2,numel(NE_order));
subAvg.Shu.LR_A = NaN(500,600,numel(NE_order));
subAvg.Shu.LR_B = NaN(500,600,numel(NE_order));
subAvg.Shu.IRFx2_A = NaN(500,600,numel(NE_order));
subAvg.Shu.IRFx2_B = NaN(500,600,numel(NE_order));
subAvg.Shu.tA = NaN(numel(NE_order),1);
subAvg.Shu.tB = NaN(numel(NE_order),1);
subAvg.Shu.IRFx2_IRF = NaN(151,2,numel(NE_order));

for i = 1:numel(NE_order)
    subAvg.Unf.LR_perf(:,:,i) = mean(tmp.Unf.LR_perf(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.Unf.IRFx2_perf(:,:,i) = mean(tmp.Unf.IRFx2_perf(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.Shu.LR_perf(:,:,i) = mean(tmp.Shu.LR_perf(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.Shu.IRFx2_perf(:,:,i) = mean(tmp.Shu.IRFx2_perf(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.Unf.LR_A(:,:,i) = mean(tmp.Unf.LR_A(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.Unf.LR_B(:,:,i) = mean(tmp.Unf.LR_B(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.Unf.IRFx2_A(:,:,i) = mean(tmp.Unf.IRFx2_A(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.Unf.IRFx2_B(:,:,i) = mean(tmp.Unf.IRFx2_B(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.Unf.LR_tA(i) = mean([unfiltered.LR_tA{NE_order(i).Runs}]);
    subAvg.Unf.LR_tB(i) = mean([unfiltered.LR_tB{NE_order(i).Runs}]);
    subAvg.Unf.IRFx2_IRF(:,:,i) = mean(cat(3,unfiltered.IRFx2_IRF{NE_order(i).Runs}),3);
    subAvg.Shu.LR_A(:,:,i) = mean(tmp.Shu.LR_A(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.Shu.LR_B(:,:,i) = mean(tmp.Shu.LR_B(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.Shu.IRFx2_A(:,:,i) = mean(tmp.Shu.IRFx2_A(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.Shu.IRFx2_B(:,:,i) = mean(tmp.Shu.IRFx2_B(:,:,NE_order(i).Runs),3,'omitnan');
    subAvg.Shu.LR_tA(i) = mean([shuffled.LR_tA{NE_order(i).Runs}]);
    subAvg.Shu.LR_tB(i) = mean([shuffled.LR_tB{NE_order(i).Runs}]);
    subAvg.Shu.IRFx2_IRF(:,:,i) = mean(cat(3,shuffled.IRFx2_IRF{NE_order(i).Runs}),3);
end

plotBM = refBM;
plotBM(:,1:300) = NaN;

fig_savePath = fullfile(savePath,'ExtDataFig5');
[~, ~, ~] = mkdir(fig_savePath);

%% Fig unfiltered A

f = figure;
f_plotMap(mean(subAvg.Unf.LR_A,3,'omitnan').*plotBM,cmp=cmpbbr,bounds=[-1 1],title='LR A',clabel='');
% for i = 1:12
%     f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
% end
colorbar off;
title '';
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig5_A1.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

f = figure;
f_plotMap(mean(subAvg.Unf.LR_B,3,'omitnan').*plotBM,cmp=cmpbbr,bounds=[-1 1],title='LR B',clabel='');
% for i = 1:12
%     f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
% end
colorbar off;
title '';
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig5_A2.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

%% Fig unfiltered B

barData = {};
barData{1} = subAvg.Unf.LR_tA / 10;
barData{2} = subAvg.Unf.LR_tB / 10;

f = figure;
[dataMean, dataSEM] = f_plotBar(barData,colors=[c_Ca;c_GRAB],legend={'tA','tB'},ylabel='r',title='Timing Coefficients')
ylim([-0.2, 1]);
saveas(f, fullfile(fig_savePath, 'ExtDataFig5_B.svg'));

T = table({order(NE_Idx).Mouse}',barData{1}',barData{2}', ...
    VariableNames={'Mouse','tA','tB'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig5_B.csv'));

%% Fig unfiltered C

f = figure;
f_plotMap(mean(subAvg.Unf.LR_perf,3,'omitnan').*plotBM,cmp=cmpvir,bounds=[0 1],title='LR Performance',clabel='r');
% for i = 1:12
%     f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
% end
colorbar off;
title '';
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig5_C.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

%% Fig unfiltered D

f = figure;
f_plotMap((mean(subAvg.Unf.LR_perf,3,'omitnan')-mean(subAvg.Fig2.g_LR_perf,3,'omitnan')).*plotBM,cmp=cmpbbr,bounds=[-1 1],title='LR vs. Global IRF',clabel='\Deltar');
% for i = 1:12
%     f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
% end
colorbar off;
title '';
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig5_D.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

%% Fig unfiltered E

f = figure;
f_plotMap(mean(subAvg.Unf.IRFx2_A,3,'omitnan').*plotBM,cmp=cmpbbr,bounds=[-1 1],title='IRFx2 A',clabel='');
% for i = 1:12
%     f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
% end
colorbar off;
title '';
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig5_E1.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

f = figure;
f_plotMap(mean(subAvg.Unf.IRFx2_B,3,'omitnan').*plotBM,cmp=cmpbbr,bounds=[-1 1],title='IRFx2 B',clabel='');
% for i = 1:12
%     f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
% end
colorbar off;
title '';
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig5_E2.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

%% Fig unfiltered F

f = figure;
meanSig1 = mean(subAvg.Unf.IRFx2_IRF(:,1,:),3);
error1 = std(subAvg.Unf.IRFx2_IRF(:,1,:),0,3)/sqrt(numel(NE_order));
f_plotLineError(-5:0.1:10,meanSig1,error1,color=c_Ca);
meanSig2 = mean(subAvg.Unf.IRFx2_IRF(:,2,:),3);
error2 = std(subAvg.Unf.IRFx2_IRF(:,2,:),0,3)/sqrt(numel(NE_order));
f_plotLineError(-5:0.1:10,meanSig2,error2,color=c_GRAB);
xlim([-2 7]);
xlabel('Time (s)');
ylabel('a.u.');
legend('','IRF(t_0^A,\tau_A)','','IRF(t_0^B,\tau_B)');
title('Double IRF');
set(gca,'FontSize',14);

t = round(-5:0.1:10,2);
idx = 30:120;

T = table(t(idx)', meanSig1(idx), meanSig2(idx), error1(idx), error2(idx), ...
    VariableNames={'Time','IRF_Ca','IRF_NE','SEM_Ca','SEM_NE'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig5_F.csv'));
saveas(f, fullfile(fig_savePath, 'ExtDataFig5_F.svg'));

%% Fig unfiltered G

f = figure;
f_plotMap(mean(subAvg.Unf.IRFx2_perf,3,'omitnan').*plotBM,cmp=cmpvir,bounds=[0 1],title='IRFx2 Performance',clabel='r');
% for i = 1:12
%     f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
% end
colorbar off;
title '';
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig5_G.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

%% Fig unfiltered H

f = figure;
f_plotMap((mean(subAvg.Unf.IRFx2_perf,3,'omitnan')-mean(subAvg.Fig2.g_IRFx2_perf,3,'omitnan')).*plotBM,cmp=cmpbbr,bounds=[-1 1],title='IRFx2 vs. Global IRF',clabel='\Deltar');
% for i = 1:12
%     f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
% end
colorbar off;
title '';
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig5_H.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

%% Fig shuffled I

f = figure;
f_plotMap(mean(subAvg.Shu.LR_A,3,'omitnan').*plotBM,cmp=cmpbbr,bounds=[-1 1],title='LR A',clabel='');
% for i = 1:12
%     f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
% end
colorbar off;
title '';
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig5_I1.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

f = figure;
f_plotMap(mean(subAvg.Shu.LR_B,3,'omitnan').*plotBM,cmp=cmpbbr,bounds=[-1 1],title='LR B',clabel='');
% for i = 1:12
%     f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
% end
colorbar off;
title '';
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig5_I2.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

%% Fig shuffled J

barData = {};
barData{1} = subAvg.Shu.LR_tA / 10;
barData{2} = subAvg.Shu.LR_tB / 10;

f = figure;
[dataMean, dataSEM] = f_plotBar(barData,colors=[c_Ca;c_GRAB],legend={'tA','tB'},ylabel='r',title='Timing Coefficients')
ylim([-0.2, 1]);
saveas(f, fullfile(fig_savePath, 'ExtDataFig5_J.svg'));

T = table({order(NE_Idx).Mouse}',barData{1}',barData{2}', ...
    VariableNames={'Mouse','tA','tB'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig5_J.csv'));

%% Fig shuffled K

f = figure;
f_plotMap((mean(subAvg.Shu.LR_perf,3,'omitnan')-mean(subAvg.Fig1.inv_perf,3,'omitnan')).*plotBM,cmp=cmpbbr,bounds=[-1 1],title='LR vs. Global IRF',clabel='\Deltar');
% for i = 1:12
%     f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
% end
colorbar off;
title '';
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig5_K.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

%% Fig shuffled L

f = figure;
f_plotMap(mean(subAvg.Shu.IRFx2_A,3,'omitnan').*plotBM,cmp=cmpbbr,bounds=[-1 1],title='IRFx2 A',clabel='');
% for i = 1:12
%     f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
% end
colorbar off;
title '';
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig5_L1.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

f = figure;
f_plotMap(mean(subAvg.Shu.IRFx2_B,3,'omitnan').*plotBM,cmp=cmpbbr,bounds=[-1 1],title='IRFx2 B',clabel='');
% for i = 1:12
%     f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
% end
colorbar off;
title '';
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig5_L2.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

%% Fig shuffled M

f = figure;
meanSig1 = mean(subAvg.Shu.IRFx2_IRF(:,1,:),3);
error1 = std(subAvg.Shu.IRFx2_IRF(:,1,:),0,3)/sqrt(numel(NE_order));
f_plotLineError(-5:0.1:10,meanSig1,error1,color=c_Ca);
meanSig2 = mean(subAvg.Shu.IRFx2_IRF(:,2,:),3);
error2 = std(subAvg.Shu.IRFx2_IRF(:,2,:),0,3)/sqrt(numel(NE_order));
f_plotLineError(-5:0.1:10,meanSig2,error2,color=c_GRAB);
xlim([-2 7]);
xlabel('Time (s)');
ylabel('a.u.');
legend('','IRF(t_0^A,\tau_A)','','IRF(t_0^B,\tau_B)');
title('Double IRF');
set(gca,'FontSize',14);

t = round(-5:0.1:10,2);
idx = 30:120;

T = table(t(idx)', meanSig1(idx), meanSig2(idx), error1(idx), error2(idx), ...
    VariableNames={'Time','IRF_Ca','IRF_NE','SEM_Ca','SEM_NE'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig5_M.csv'));
saveas(f, fullfile(fig_savePath, 'ExtDataFig5_M.svg'));

%% Fig shuffled N

f = figure;
f_plotMap((mean(subAvg.Shu.IRFx2_perf,3,'omitnan')-mean(subAvg.Fig1.inv_perf,3,'omitnan')).*plotBM,cmp=cmpbbr,bounds=[-1 1],title='IRFx2 vs. Global IRF',clabel='\Deltar');
% for i = 1:12
%     f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
% end
colorbar off;
title '';
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig5_N.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

%% Fig shuffled O

barData = {};
barData{1} = squeeze(mean(subAvg.Fig1.inv_perf.*plotBM,[1,2],'omitnan'));
barData{2} = squeeze(mean(subAvg.Fig1.SSp_perf.*plotBM,[1,2],'omitnan'));
barData{3} = squeeze(mean(subAvg.Fig1.var_perf.*plotBM,[1,2],'omitnan'));
barData{4} = squeeze(mean(subAvg.Fig2.g_LR_perf.*plotBM,[1,2],'omitnan'));
barData{5} = squeeze(mean(subAvg.Fig2.g_IRFx2_perf.*plotBM,[1,2],'omitnan'));
barData{6} = squeeze(mean(subAvg.Unf.LR_perf.*plotBM,[1,2],'omitnan'));
barData{7} = squeeze(mean(subAvg.Unf.IRFx2_perf.*plotBM,[1,2],'omitnan'));
barData{8} = squeeze(mean(subAvg.Shu.LR_perf.*plotBM,[1,2],'omitnan'));
barData{9} = squeeze(mean(subAvg.Shu.IRFx2_perf.*plotBM,[1,2],'omitnan'));

f = figure;
[dataMean, dataSEM] = f_plotBar(barData,colors=[repmat(c_Yellow,3,1);repmat(c_darkCyan,2,1);repmat(c_Ca,2,1);repmat(c_pupil,2,1)],legend={'Invariant','SSp','Variant','LR','IRFx2','unfiltered LR','unfiltered IRFx2','shuffled LR','shuffled IRFx2'},ylabel='r',title='Model Performance Comparison')
ylim([0, 1]);
legend('off');

[h,p] = f_kstest(barData,0.01);

T = NaN(numel(order),9);
T(:,1) = barData{1};
T(:,2) = barData{2};
T(:,3) = barData{3};
T(NE_Idx,4) = barData{4};
T(NE_Idx,5) = barData{5};
T(NE_Idx,6) = barData{6};
T(NE_Idx,7) = barData{7};
T(NE_Idx,8) = barData{8};
T(NE_Idx,9) = barData{9};

saveas(f, fullfile(fig_savePath, 'ExtDataFig5_O.svg'));
T = table({order.Mouse}',T(:,1),T(:,2),T(:,3),T(:,4),T(:,5),T(:,6),T(:,7),T(:,8),T(:,9), ...
    VariableNames={'Mouse','global','SSp','variant','Linear Regression', ...
    'Double IRF','Unf Linear Regression','Unf Double IRF', ...
    'Shuffled Linear Regression','Shuffled Double IRF'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig5_O.csv'));

labels = {'global','SSp','variant','Linear Regression','Double IRF', ...
    'Unfiltered Linear Regression','Unfiltered Double IRF','Shuffled Linear Regression','Shuffled Double IRF'};
T = table;
for i = 1:numel(labels); T.(labels{i}) = p(:,i); end
writetable(T, fullfile(fig_savePath, 'ExtDataFig5_O_p.csv'));
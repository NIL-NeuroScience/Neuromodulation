close all

% calculate subject averages

NE_order = order(NE_Idx);

subAvg.FigE9.Ca_vs_HbT = NaN(12,12,numel(NE_order));
subAvg.FigE9.Ca_vs_HbT_reg = NaN(12,12,numel(NE_order));
subAvg.FigE9.global_Ca_vs_HbT_reg = NaN(12,12,numel(NE_order));
subAvg.FigE9.lowNE_HbT_reg = NaN(12,12,numel(NE_order));
subAvg.FigE9.highNE_HbT_reg = NaN(12,12,numel(NE_order));
subAvg.FigE9.global_lowNE_HbT_reg = NaN(12,12,numel(NE_order));
subAvg.FigE9.global_highNE_HbT_reg = NaN(12,12,numel(NE_order));

for i = 1:numel(NE_order)
    subAvg.FigE9.Ca_vs_HbT(:,:,i) = mean(cat(3,NE_reg.FC_R_HbT{NE_order(i).Runs}),3);
    subAvg.FigE9.Ca_vs_HbT_reg(:,:,i) = mean(cat(3,NE_reg.FC_R_HbT_reg{NE_order(i).Runs}),3);
    subAvg.FigE9.global_Ca_vs_HbT_reg(:,:,i) = mean(cat(3,NE_reg.global_FC_R_HbT_reg{NE_order(i).Runs}),3);
    subAvg.FigE9.lowNE_HbT_reg(:,:,i) = mean(cat(3,NE_reg.lowNE_FC{NE_order(i).Runs}),3);
    subAvg.FigE9.highNE_HbT_reg(:,:,i) = mean(cat(3,NE_reg.highNE_FC{NE_order(i).Runs}),3);
    subAvg.FigE9.global_lowNE_HbT_reg(:,:,i) = mean(cat(3,NE_reg.global_lowNE_FC{NE_order(i).Runs}),3);
    subAvg.FigE9.global_highNE_HbT_reg(:,:,i) = mean(cat(3,NE_reg.global_highNE_FC{NE_order(i).Runs}),3);
end

plotBM = refBM;
plotBM(:,1:300) = NaN;

fig_savePath = fullfile(savePath,'ExtDataFig9');
[~, ~, ~] = mkdir(fig_savePath);

%% Fig NE_reg B

runData = load('sub-Thy1-288_ses-24-05-06_run-01_irun-01_behavior+ophys.mat');

t = 15:6:585;

f = figure;hold on;
plot(t,squeeze(runData.NE_reg.FC.HbT_reg(2,5,:)),color=[0,0.7,0.7]);
plot([55,55],[0,1],'-k');
ylim([0,1]);
xlim([50,550]);
saveas(f, fullfile(fig_savePath, 'ExtDataFig9_B.svg'));

T = table(t', squeeze(runData.NE_reg.FC.HbT_reg(2,5,:)), ...
    VariableNames={'Time','FC'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig9_B.csv'));

%% Fig NE_reg C

f = figure;
f_plotFC(mean(subAvg.FigE9.Ca_vs_HbT,3),0,cmp=cmpbbr,clim=[-1 1],title='HbT',clabel='r');
colorbar off;title '';
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig9_C1.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);
writetable(table(mean(subAvg.FigE9.Ca_vs_HbT,3)), fullfile(fig_savePath, 'ExtDataFig9_C1.csv'));

f = figure;
f_plotFC(mean(subAvg.FigE9.Ca_vs_HbT_reg,3),0,cmp=cmpbbr,clim=[-1 1],title='NE \\ HbT',clabel='r');
colorbar off;title '';
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig9_C2.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);
writetable(table(mean(subAvg.FigE9.Ca_vs_HbT_reg,3)), fullfile(fig_savePath, 'ExtDataFig9_C2.csv'));

f = figure;
f_plotFC(mean(subAvg.FigE9.global_Ca_vs_HbT_reg,3),0,cmp=cmpbbr,clim=[-1 1],title='Global NE \\ HbT',clabel='r');
colorbar off;title '';
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig9_C3.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);
writetable(table(mean(subAvg.FigE9.global_Ca_vs_HbT_reg,3)), fullfile(fig_savePath, 'ExtDataFig9_C3.csv'));

%% Fig NE_reg D

f = figure;
f_plotFC(mean(subAvg.FigE9.lowNE_HbT_reg,3),1,cmp=cmpvir,clim=[0 1],title='Low NE Ca++ Connectivity',clabel='r');
colorbar off;title '';
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig9_D1.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);
writetable(table(mean(subAvg.FigE9.global_Ca_vs_HbT_reg,3)), fullfile(fig_savePath, 'ExtDataFig9_D1.csv'));

f = figure;
f_plotFC(mean(subAvg.FigE9.highNE_HbT_reg,3),1,cmp=cmpvir,clim=[0 1],title='High NE Ca++ Connectivity',clabel='r');
colorbar off;title '';
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig9_D2.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);
writetable(table(mean(subAvg.FigE9.global_Ca_vs_HbT_reg,3)), fullfile(fig_savePath, 'ExtDataFig9_D2.csv'));

MouseID = {log(strcmp({log.GRAB},'GRAB_NE')).Mouse};
g1 = cat(3,NE_reg.lowNE_FC{strcmp({log.GRAB},'GRAB_NE')});
g1 = reshape(g1,144,[])';
g2 = cat(3,NE_reg.highNE_FC{strcmp({log.GRAB},'GRAB_NE')});
g2 = reshape(g2,144,[])';
[h,p] = f_lme(MouseID,g1,g2,0.05);

idx = logical(eye(12));
h(idx(:)) = 0;

f = figure;
f_plotFC(mean(subAvg.FigE9.highNE_HbT_reg,3)-mean(subAvg.FigE9.lowNE_HbT_reg,3),0,cmp=cmpbbr,clim=0.1*[-1 1],title='High-Low NE',clabel='\Deltar');
f_overlayStats_FC(reshape(h,12,12));
title('');colorbar off;
saveas(f, fullfile(fig_savePath, 'ExtDataFig9_D3.svg'));
writetable(table(mean(subAvg.FigE9.highNE_HbT_reg,3)-mean(subAvg.FigE9.lowNE_HbT_reg,3)), fullfile(fig_savePath, 'ExtDataFig9_D3.csv'));
writetable(table(reshape(p,12,12)), fullfile(fig_savePath, 'ExtDataFig9_D3_P.csv'));

f = figure;
f_plotFC(mean(subAvg.FigE9.highNE_HbT_reg,3)-mean(subAvg.FigE9.lowNE_HbT_reg,3),0,cmp=cmpbbr,clim=0.1*[-1 1],title='High-Low NE',clabel='\Deltar');
colorbar off;
title '';
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig9_D3.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

%% Fig NE_reg E

f = figure;
f_plotFC(mean(subAvg.FigE9.global_lowNE_HbT_reg,3),1,cmp=cmpvir,clim=[0 1],title='Low NE Ca++ Connectivity',clabel='r');
colorbar off;title '';
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig9_E1.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);
writetable(table(mean(subAvg.FigE9.global_Ca_vs_HbT_reg,3)), fullfile(fig_savePath, 'ExtDataFig9_E1.csv'));

f = figure;
f_plotFC(mean(subAvg.FigE9.global_highNE_HbT_reg,3),1,cmp=cmpvir,clim=[0 1],title='High NE Ca++ Connectivity',clabel='r');
colorbar off;title '';
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig9_E2.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);
writetable(table(mean(subAvg.FigE9.global_Ca_vs_HbT_reg,3)), fullfile(fig_savePath, 'ExtDataFig9_E2.csv'));

MouseID = {log(strcmp({log.GRAB},'GRAB_NE')).Mouse};
g1 = cat(3,NE_reg.global_lowNE_FC{strcmp({log.GRAB},'GRAB_NE')});
g1 = reshape(g1,144,[])';
g2 = cat(3,NE_reg.global_highNE_FC{strcmp({log.GRAB},'GRAB_NE')});
g2 = reshape(g2,144,[])';
[h,p] = f_lme(MouseID,g1,g2,0.05);

idx = logical(eye(12));
h(idx(:)) = 0;

f = figure;
f_plotFC(mean(subAvg.FigE9.global_highNE_HbT_reg,3)-mean(subAvg.FigE9.global_lowNE_HbT_reg,3),0,cmp=cmpbbr,clim=0.1*[-1 1],title='High-Low NE',clabel='\Deltar');
f_overlayStats_FC(reshape(h,12,12));
title('');colorbar off;
saveas(f, fullfile(fig_savePath, 'ExtDataFig9_E3.svg'));
writetable(table(mean(subAvg.FigE9.global_highNE_HbT_reg,3)-mean(subAvg.FigE9.global_lowNE_HbT_reg,3)), fullfile(fig_savePath, 'ExtDataFig9_E3.csv'));
writetable(table(reshape(p,12,12)), fullfile(fig_savePath, 'ExtDataFig9_E3_P.csv'));

f = figure;
f_plotFC(mean(subAvg.FigE9.global_highNE_HbT_reg,3)-mean(subAvg.FigE9.global_lowNE_HbT_reg,3),0,cmp=cmpbbr,clim=0.1*[-1 1],title='High-Low NE',clabel='\Deltar');
colorbar off;
title '';
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig9_E3.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

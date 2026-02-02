close all

% calculate subject averages

NE_order = order(NE_Idx);

subAvg.FigE4.XC_gfp_HD_HbT_allen = NaN(201,12,M);

for i = 1:M
    try 
        subAvg.FigE4.XC_gfp_HD_HbT_allen(:,:,i) = mean(cat(3,GRAB_FC.gfp_HD_vs_HbT_low{order(i).Runs}),3);
    end
end

plotBM = refBM;
plotBM(:,1:300) = NaN;

fig_savePath = fullfile(savePath,'ExtDataFig4');
[~, ~, ~] = mkdir(fig_savePath);

%% Extended data fig 4 E

t = 10:-0.1:-10;

f = figure;
meanSig = mean(subAvg.FigE4.XC_gfp_HD_HbT_allen(:,2,NE_Idx),3);
error = std(subAvg.FigE4.XC_gfp_HD_HbT_allen(:,2,NE_Idx),0,3)/sqrt(M_NE);
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
meanSig = mean(subAvg.FigE4.XC_gfp_HD_HbT_allen(:,5,NE_Idx),3);
error = std(subAvg.FigE4.XC_gfp_HD_HbT_allen(:,5,NE_Idx),0,3)/sqrt(M_NE);
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
meanSig = mean(subAvg.FigE4.XC_gfp_HD_HbT_allen(:,12,NE_Idx),3);
error = std(subAvg.FigE4.XC_gfp_HD_HbT_allen(:,12,NE_Idx),0,3)/sqrt(M_NE);
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
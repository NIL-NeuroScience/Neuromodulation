close all

load(fullfile(f_path,'Figures/results/HD/HD_metadata.mat'),'metadata');

mice = fieldnames(metadata);

subAvg.FigE10.XC_gfp_HD_HbT_low = NaN(201,numel(mice));
subAvg.FigE10.R_gfp_HD_HbT_low = cell(3,1);

for i = 1:3
    subAvg.FigE10.XC_gfp_HD_HbT_low(:,i) = mean([metadata.(mice{i}).XC_gfp_HD_HbT_low],2);
    subAvg.FigE10.R_gfp_HD_HbT_low{i} = mean(cat(3,metadata.(mice{i}).R_gfp_HD_HbT_low),3);
end

fig_savePath = fullfile(savePath,'ExtDataFig10');
[~, ~, ~] = mkdir(fig_savePath);

%% plot Fig HD A

run = 129;

pupil = Behavior.signals{run}(:,4);
t = 0.1:0.1:600;

f = figure;hold on;
plot(t,pupil,color=c_pupil);
plot([100,100],[0,1],'-k',lineWidth=2);
plot(179,pupil(1790),'o');
plot(243,pupil(2430),'o');
xlim([50, 350]);
box off;

pupil(1790)
pupil(2430)

idx = 500:3500;
saveas(f, fullfile(fig_savePath, 'ExtDataFig10_A.svg'));

T = table(t(idx)',pupil(idx), ...
    VariableNames={'Time','Pupil Diameter'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig10_A.csv'));

%% plot Fig HD C

runData = load(fullfile(f_path,'Figures/results/HD/Ca_masks.mat'));
runData = runData.data;

mask = runData.mask;
mask(isnan(mask)) = 0;
mask = bwboundaries(mask);

f = figure;hold on;
imagesc(runData.int_Ca);
axis image off;
clim([0,30000]);
colormap cmpinf;
set(gca,YDir='reverse');
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig10_C1.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

for i = 1:2
    plot(mask{i}(:,2),mask{i}(:,1),'-',color=0.7*[1,1,1],lineWidth=2);
end

set(gca,YDir='reverse');
saveas(f,fullfile(fig_savePath,'ExtDataFig10_C1.svg'));


f = figure;hold on;
imagesc(runData.std_Ca);
axis image off;
clim([0,5]);
colormap cmpinf;
set(gca,YDir='reverse');
exportgraphics(f, fullfile(fig_savePath,'ExtDataFig10_C2.jpg'),'Resolution',300,'BackgroundColor',[1 1 1]);

for i = 1:2
    plot(mask{i}(:,2),mask{i}(:,1),'-',color=0.7*[1,1,1],lineWidth=2);
end

set(gca,YDir='reverse');
saveas(f,fullfile(fig_savePath,'ExtDataFig10_C2.svg'));

%% plot Fig HD D

run = 138;
Ca = spectra.Ca{run};

t = 0.1:0.1:600;

f = figure;hold on;
plot(t,Ca(:,2),Color=0.7*c_Ca);
plot(t,Ca(:,3)-10,Color=c_Ca);
plot([100,100],[0, 5],'-k',lineWidth=2);
box off;

saveas(f, fullfile(fig_savePath, 'ExtDataFig10_D.svg'));

T = table(t',Ca(:,2),Ca(:,3), ...
    VariableNames={'Time','Ca_MOs','Ca_SSpll'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig10_D.csv'));

%% plot Fig HD E

mask = runData.mask;
mask(isnan(mask)) = 0;
std_Ca = runData.std_Ca(logical(mask(:)));
int_Ca = runData.int_Ca(logical(mask(:)));

lm = fitlm(int_Ca,std_Ca);
lm = table2array(lm.Coefficients);

skip = 10;
f = figure;hold on;
scatter(int_Ca(1:skip:end),std_Ca(1:skip:end),50,'filled',MarkerFaceAlpha=0.1,MarkerFaceColor=[0,0.7,0.7]);
plot([0,30000],lm(2,1)*[0,30000]+lm(1,1),'-k',lineWidth=2);
xlim([0,30000]);
ylim([0,5]);
box off;

saveas(f, fullfile(fig_savePath, 'ExtDataFig10_E.svg'));

T = table(int_Ca,std_Ca, ...
    VariableNames={'Ca_intensity','Ca_std_dev'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig10_E.csv'));

f_corr(int_Ca,std_Ca,1)

%% plot Fig HD F

data2P = load(fullfile(f_path,'Figures/results/HD/test2P.mat'));

t = (1:numel(data2P.dil))/15;

f = figure;hold on;
plot(t,data2P.dil,color=c_HbT);
plot(t,data2P.NE_sig-0.40,color=0.7*[1,1,1]);
plot([50,50],[0,0.1],'-k');
saveas(f, fullfile(fig_savePath, 'ExtDataFig10_F.svg'));

T = table(t',data2P.dil,data2P.NE_sig, ...
    VariableNames={'Time','Vessel Dilation','NE_mutant'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig10_F.csv'));

%% plot Fig HD G

dataMutant = load(fullfile(f_path,'Figures/results/HD/mutant_data_B6_248_250204.mat'));
dataMutant = dataMutant.mutant_data;

t = 0.1:0.1:600;
region = 12;

f = figure;hold on;
plot(t,dataMutant.HbT(:,region),color=c_HbT);
plot(t,dataMutant.gfp(:,region)-10,color=c_Ca);
plot(t,dataMutant.gfp_HD(:,region)-20,color=0.7*[1,1,1]);
plot([120,120],[0,5],'-k');
xlim([100,400]);

idx = 1000:4000;

saveas(f, fullfile(fig_savePath, 'ExtDataFig10_G.svg'));

T = table(t(idx)',dataMutant.HbT(idx,region),dataMutant.gfp(idx,region),dataMutant.gfp_HD(idx,region), ...
    VariableNames={'Time','HbT','gfp','gfp_HD'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig10_G.csv'));

%% plot Fig HD H

t = -5:1/15:5;

f = figure;
plot(t,xcorr(data2P.dil,data2P.NE_sig,75,'normalized'),color=0.7*[1,1,1]);
ylim(0.5*[-1,1]);
xlim([-5,5]);
box off;
saveas(f, fullfile(fig_savePath, 'ExtDataFig10_H.svg'));

T = table(t',xcorr(data2P.dil,data2P.NE_sig,75,'normalized'), ...
    VariableNames={'Time','r'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig10_H.csv'));

%% plot Fig HD I
t = 10:-0.1:-10;

f = figure;
meanSig1 = mean(subAvg.FigE2.XC_gfp_HD_HbT(:,NE_Idx),2);
error1 = std(subAvg.FigE2.XC_gfp_HD_HbT(:,NE_Idx),0,2)/sqrt(M_NE);
f_plotLineError(t,meanSig1,error1,color=c_GRAB);
meanSig2 = mean(subAvg.FigE2.XC_gfp_HD_HbT(:,ACh_Idx),2);
error2 = std(subAvg.FigE2.XC_gfp_HD_HbT(:,ACh_Idx),0,2)/sqrt(M_ACh);
f_plotLineError(t,meanSig2,error2,color=c_Orange);
meanSig3 = mean(subAvg.FigE10.XC_gfp_HD_HbT_low,2);
error3 = std(subAvg.FigE10.XC_gfp_HD_HbT_low,0,2)/sqrt(3);
f_plotLineError(t,meanSig3,error3,color=0.5*[1,1,1],lineWidth=3);
xlim([-5 5]);
ylim(0.6*[-1,1]);
xlabel('Time (s)');
ylabel('r');
set(gca,'FontSize',14);
title('x vs. HbT');
legend('','NE','','ACh','','Ca^2^+');
saveas(f, fullfile(fig_savePath, 'ExtDataFig10_I.svg'));

T = table(t',meanSig1,meanSig2,meanSig3,error1,error2,error3, ...
    VariableNames={'Time','mean_NE','mean_ACh','mean_NEmut','SEM_NE','SEM_ACh','SEM_NEmut'});
writetable(T, fullfile(fig_savePath, 'ExtDataFig10_I.csv'));

%%

f = figure;
f_plotMap(mean(subAvg.FigE2.R_gfp_HD_HbT(:,:,NE_Idx),3,'omitnan').*plotBM,cmp=cmpbbr,clim=1*[-1 1],title='NE vs. Ca^2^+',clabel='r');
for i = 1:12
    f_plotAllenRegion(i,2,linewidth=2,color=[0 0 0]);
end

%%
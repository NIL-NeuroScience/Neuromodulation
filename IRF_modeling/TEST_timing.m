%% IRF model

tic
[perf] = f_LR_varWeights(HbT,rfp_HD,gfp_HD,[-50,100],brain_mask,ds=32,method='fft');
toc

%% IRF model

tic
[perf,IRF] = f_2alphaDeconvolve(HbT,rfp_HD,gfp_HD,[-50,100],brain_mask,ds=32,method='direct');
toc

%%

tic
[IRFx1_SSp2.perf,IRFx1_SSp2.IRF] = TEST_1xIRF(HbT,rfp_HD,win,fs,brain_mask,4,[],4);
toc

%%
point = [197,108];
binSize = 30;
tmp_Ca = rfp_HD(point(1):point(1)+binSize,point(2):point(2)+binSize,:);
tmp_HbT = HbT(point(1):point(1)+binSize,point(2):point(2)+binSize,:);

tmpBM = ones(size(tmp_Ca,[1,2]));

%% old IRF model

tic
[IRFx1_SSp.perf,IRFx1_SSp.IRF] = f_1xIRF(tmp_HbT,tmp_Ca,win,fs,tmpBM,1,corrWin*fs,numThreads);
toc

%%

tic
[IRFx1_SSp2.perf,IRFx1_SSp2.IRF] = TEST_1xIRF(tmp_HbT,tmp_Ca,win,fs,tmpBM,1,corrWin*fs,numThreads);
toc

%% test

t0 = 0.1;
tau1 = 0.5;
tau2 = 0.53;
A = 1;
B = -1;

fun = @(params)f_hrf_cost_func(params(1),params(2),params(3),params(4),params(5),fs,win,Ca_mat,design_HbT);

options = optimset('MaxFunEvals',25000,'MaxIter',500,'Display','iter','Algorithm','active-set','FunValCheck','on','Display','off');
params = fmincon(fun,[t0, tau1, tau2, A, B],[],[],[],[],[0,0.01,0.01,-Inf,-Inf],[win(2),win(2),win(2),Inf,Inf],[],options);

%%

function J = f_hrf_cost_func(t0,tau1,tau2,A,B,sr,range,X_mat,y)
    hrf = f_alpha_IRF(t0,tau1,tau2,A,B,sr,range);
    
    conv_result = X_mat*hrf;
    J = sqrt(mean((y - conv_result).^2,'all'));
end

function IRF = f_alpha_IRF(t0,tau1,tau2,A,B,sr,range)
    hrf_l = range*sr;
    tr = (((hrf_l(1):hrf_l(2))/sr)-t0)';
    D = ((tr)./tau1).^3.*exp(-(tr)./tau1);
    D(tr<0) = 0;
    C = ((tr)./tau2).^3.*exp(-(tr)./tau2);
    C(tr<0) = 0;
    
    IRF = A.*D + B.*C;
end
%% View the best fitting parameters
accreRoot = '/gpfs22';
accreHome = '/home/middlepg';
accreScratch = '/scratch/middlepg';
if isdir(fullfile(accreScratch))
    matRoot = fullfile(accreRoot,accreHome,'m-files');
    modelRoot = fullfile(accreScratch,'perceptualchoice_stop_model');
    environ = 'accre';
else
    matRoot = '/Volumes/HD-1/Users/paulmiddlebrooks/matlab';
    modelRoot = '/Volumes/HD-1/Users/paulmiddlebrooks/perceptualchoice_stop_model';
    environ = 'local';
end

% Define variables
subj        = 1;
dt          = 1;
trialVar    = true;
optimScope  = 'all';
choiceMech  = 'race'; %{'race', 'li', 'ffi'};
stopMech    = 'li';
model       = 191;%[79 191 245 478];

fileStr.root   = fullfile(modelRoot,'data/2017-05-08/preproc01/subj%.2d/dt%d/%s/%s/');
fileStr.output = 'allFValAndBestX_%sTrials_model%.3d.txt';
% ===============================================================================================


if trialVar
    trialVarStr = 'trialvar';
else
    trialVarStr = 'notrialvar';
end

switch optimScope
    case 'go'
        modelArch = sprintf('c%s',choiceMech);
        stopMech = 'none';
    case {'stop','all'}
        modelArch = sprintf('c%s_i%s',choiceMech,stopMech);
end

% Loop over subjects and models
iSubj = subj;
iModel = model;


% Identify all final log files for this subject and model
allFinalLog = regexpdir(sprintf(fileStr.root,iSubj,dt,trialVarStr,modelArch),sprintf('finalLog_%sTrials_model%.3d_startVal.*_started.*.mat$',lower(optimScope),iModel))
allSAM      = regexpdir(sprintf(fileStr.root,iSubj,dt,trialVarStr,modelArch),sprintf('SAM_%sTrials_model%.3d_startVal.*_started.*.mat$',lower(optimScope),iModel));
nFile       = size(allSAM,1);


ds          = dataset({[1:nFile]','FileIndex'}, ...
    {cell(nFile,1),'FileNameLog'}, ...
    {cell(nFile,1),'FileNameSAM'}, ...
    {nan(nFile,1),'StartValueIndex'}, ...
    {nan(nFile,1),'FunctionValue'}, ...
    {cell(nFile,1),'ElapsedTime'}, ...
    {false(nFile,1),'ExitFlag'}, ...
    {cell(nFile,1),'BestX'}, ...
    {cell(nFile,1),'StartX'});

for iFile = 1:nFile
    
    % Get and enter final log filename
    [~,fileNameLog]         = fileparts(allFinalLog{iFile});
    ds.FileNameLog{iFile}   = fileNameLog;
    
    % Get and enter final SAM filename
    [~,fileNameSam]         = fileparts(allSAM{iFile});
    ds.FileNameSAM{iFile}   = fileNameSam;
    
    % Load file
    load(allFinalLog{iFile},'fVal','tElapse','exitFlag','X');
    
    % Enter optimization function value
    ds.FunctionValue(iFile) = fVal;
    
    % Enter elapsed time
    hours = floor(tElapse / 3600);
    tElapse = tElapse - hours * 3600;
    mins = floor(tElapse / 60);
    secs = tElapse - mins * 60;
    hms = sprintf('%02d:%02d:%05.2f', hours, mins, secs);
    ds.ElapsedTime{iFile} = hms;
    
    % Enter exit flag
    ds.ExitFlag(iFile) = exitFlag;
    
    % Enter best-fitting parameter values
    ds.BestX{iFile} = X;
    
    % Load SAM and identify and enter starting value
    settings  = regexp(fileNameLog, 'finalLog_(\w*)Trials_model(\w*)_startVal(\w*)_started(.*)', 'tokens');
    iStartVal = str2num(settings{1}{3});
    
    ds.StartValueIndex(iFile) = iStartVal;
    
    % Enter initial parameter values
    load(allSAM{iFile},'SAM');
    ds.StartX{iFile} = SAM.optim.x0(iStartVal,:);
    
end

% Sort dataset by function value
dSort = sortrows(ds,'FunctionValue')






%% Ensure errors and corrects match up to real data

subj                    = 1; %[1 2];
model                   = [79,191,245,478];
model                   = 245; %[478];
trialVar                = true;
simScope                = 'go';
choiceMech              = 'race';
stopMech                = 'none';
fileStr.root            = strcat(modelRoot,'/data/2017-05-08/preproc01/subj%.2d/dt5/%s/%s/');
doPlot                  = true;
doSave                  = true;
doStartParCluster       = false;

cd(fullfile(modelRoot,'src/code/2017-05-08/matlab'));

% Load the general SAM file
% =========================================================================
nameSAMGeneral    = 'SAM_%sTrials.mat';
modelArch = sprintf('c%s',choiceMech);
rootDir           = fileStr.root;
trialVarStr = 'trialvar';
load(fullfile(sprintf(rootDir,subj,trialVarStr,modelArch),sprintf(nameSAMGeneral,simScope)));


%%

% sum(simData.cnd == 1 & simData.acc == 1)
%
% sum(prd.nGoCCorr(1) + prd.nGoCCorr(2))
% sum(obs.nGoCCorr(1) + obs.nGoCCorr(2))
figure(10)

for iCnd = 1 : 3
    for iStm = 1 : 2
        for iAcc = 0 : 1
            % Correct
            % simRT = (simData.rt(simData.cnd == 1 & simData.acc == 1));
            % prdRT = ([prd.rtGoCCorr{3}'; prd.rtGoCCorr{4}']);
            % obsRT = ([obs.rtGoCCorr{3}; obs.rtGoCCorr{4}]);
            
            % Error
            % simRT = (simData.rt(simData.cnd == 2 & simData.acc == 0));
            % prdRT = ([prd.rtGoCError{3}'; prd.rtGoCError{4}']);
            % obsRT = ([obs.rtGoCCorr{3}; obs.rtGoCCorr{4}]);
            
            simRT = (simData.rt(simData.cnd == iCnd & simData.stm1 == iStm & simData.acc == iAcc));
            
            iCat = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('goTrial.*GO.*r%d.*c%d',iStm,iCnd), 'once')),prd.trialCat,'Uni',0)));
            switch iAcc
                case 0
                    prdRT = (prd.rtGoCError{iCat});
                    obsRT = (obs.rtGoCError{iCat});
                case 1
                    prdRT = (prd.rtGoCCorr{iCat});
                    obsRT = (obs.rtGoCCorr{iCat});
            end
            
            simRT = sort(simRT);
            prdRT = sort(prdRT);
            obsRT = sort(obsRT);
            iSimGoRT = nan(length(min(simRT) : max(simRT)), 1);
            iPrdGoRT = nan(length(min(prdRT) : max(prdRT)), 1);
            iObsGoRT = nan(length(min(obsRT) : max(obsRT)), 1);
            iRTIndex = 1;
            for iRT = min(simRT) : max(simRT)
                iSimGoRT(iRTIndex) = sum(simRT <= iRT) / length(simRT);
                iRTIndex = iRTIndex + 1;
            end
            iRTIndex = 1;
            for iRT = min(prdRT) : max(prdRT)
                iPrdGoRT(iRTIndex) = sum(prdRT <= iRT) / length(prdRT);
                iRTIndex = iRTIndex + 1;
            end
            iRTIndex = 1;
            for iRT = min(obsRT) : max(obsRT)
                iObsGoRT(iRTIndex) = sum(obsRT <= iRT) / length(obsRT);
                iRTIndex = iRTIndex + 1;
            end
            clf
            hold all;
            plot(min(simRT):max(simRT), iSimGoRT, 'color', 'k')
            pause
            plot(min(prdRT):max(prdRT), iPrdGoRT, '--', 'color', 'r')
            pause
            plot(min(obsRT):max(obsRT), iObsGoRT, '.-','color', 'b')
            
            pause
            
        end
    end
end


%%
cumProb           = SAM.optim.cost.stat.cumProb;
minBinSize        = SAM.optim.cost.stat.minBinSize;
dt = 10;

nTrialCat = size(prd, 1);
for iTrialCat = 1 : nTrialCat
    [prd.rtQGoCCorr{iTrialCat}, ...
        prd.cumProbGoCCorr{iTrialCat}, ...
        prd.cumProbDefectiveGoCCorr{iTrialCat}, ...
        prd.probMassGoCCorr{iTrialCat}, ...
        prd.probMassDefectiveGoCCorr{iTrialCat}] = ...
        sam_bin_data(prd.rtGoCCorr{iTrialCat},prd.pGoCCorr(iTrialCat),cumProb,minBinSize,dt);
    
    [prd.rtQGoCError{iTrialCat}, ...
        prd.cumProbGoCError{iTrialCat}, ...
        prd.cumProbDefectiveGoCError{iTrialCat}, ...
        prd.probMassGoCError{iTrialCat}, ...
        prd.probMassDefectiveGoCError{iTrialCat}] = ...
        sam_bin_data(prd.rtGoCError{iTrialCat},prd.pGoCError(iTrialCat),cumProb,minBinSize,dt);
end

%%
[nanmean(prd.rtGoCCorr{1}) nanmean(prd.rtGoCError{1}) nanmean(SAM.optim.obs.rtGoCCorr{1}) nanmean(SAM.optim.obs.rtGoCError{1});nanmean(prd.rtGoCCorr{2}) nanmean(prd.rtGoCError{2}) nanmean(SAM.optim.obs.rtGoCCorr{2}) nanmean(SAM.optim.obs.rtGoCError{2});]

[prd.nGoCError(1)/(prd.nGoCCorr(1) + prd.nGoCError(1)) SAM.optim.obs.nGoCError(1)/(SAM.optim.obs.nGoCCorr(1) + SAM.optim.obs.nGoCError(1));prd.nGoCError(2)/(prd.nGoCCorr(2) + prd.nGoCError(2)) SAM.optim.obs.nGoCError(2)/(SAM.optim.obs.nGoCCorr(2) + SAM.optim.obs.nGoCError(2)) ]

old
341.8528  404.7417  339.9475  484.1042
349.9977  398.4985  364.9365  353.4453

0.1084    0.0241
0.1292    0.0745


new
341.8528  404.7417  339.9475  353.4453
349.9977  398.4985  364.9365  484.1042

0.1084    0.0658
0.1292    0.0274

%%
iSsd = 16;
iRow = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('stopTrial.*ssd%d.*',iSsd), 'once')),SAM.optim.prd.trialCat,'Uni',0)));
SAM.optim.obs(iRow,1:16)


%%

ssdArray = unique(allData.ssd);
ssdArray(1) = [];
nSSD = nan(length(ssdArray),1);

for i = 1 : length(ssdArray)
    nSSD(i) = sum(allData.ssd == ssdArray(i));
end
disp([ssdArray, nSSD])

%%

go1Trial = allData.cnd == 1 & allData.ssd == 0 & allData.stm1 == 1 & allData.resp == 1;
goRT1 = allData.rt(go1Trial);

go2Trial = allData.cnd == 2 & allData.ssd == 0 & allData.stm1 == 1 & allData.resp == 1;
goRT2 = allData.rt(go2Trial);

go3Trial = allData.cnd == 3 & allData.ssd == 0 & allData.stm1 == 1 & allData.resp == 1;
goRT3 = allData.rt(go3Trial);

go4Trial = allData.cnd == 3 & allData.ssd == 0 & allData.stm1 == 2 & allData.resp == 2;
goRT4 = allData.rt(go4Trial);

go5Trial = allData.cnd == 2 & allData.ssd == 0 & allData.stm1 == 2 & allData.resp == 2;
goRT5 = allData.rt(go5Trial);

go6Trial = allData.cnd == 1 & allData.ssd == 0 & allData.stm1 == 2 & allData.resp == 2;
goRT6 = allData.rt(go6Trial);

stop1Trial = allData.cnd == 1 & allData.ssd ~= 0 & allData.stm1 == 1 & allData.resp == 1;
stopRT1 = allData.rt(stop1Trial);

stop2Trial = allData.cnd == 2 & allData.ssd ~= 0 & allData.stm1 == 1 & allData.resp == 1;
stopRT2 = allData.rt(stop2Trial);

stop3Trial = allData.cnd == 3 & allData.ssd ~= 0 & allData.stm1 == 1 & allData.resp == 1;
stopRT3 = allData.rt(stop3Trial);

stop4Trial = allData.cnd == 3 & allData.ssd ~= 0 & allData.stm1 == 2 & allData.resp == 2;
stopRT4 = allData.rt(stop4Trial);

stop5Trial = allData.cnd == 2 & allData.ssd ~= 0 & allData.stm1 == 2 & allData.resp == 2;
stopRT5 = allData.rt(stop5Trial);

stop6Trial = allData.cnd == 1 & allData.ssd ~= 0 & allData.stm1 == 2 & allData.resp == 2;
stopRT6 = allData.rt(stop6Trial);


figure(2)
hold all
plot([nanmean(goRT1) nanmean(goRT2) nanmean(goRT3) nanmean(goRT4) nanmean(goRT5) nanmean(goRT6)])
plot([nanmean(stopRT1) nanmean(stopRT2) nanmean(stopRT3) nanmean(stopRT4) nanmean(stopRT5) nanmean(stopRT6)])

[nanmean(goRT1) nanmean(goRT2) nanmean(goRT3) nanmean(goRT4) nanmean(goRT5) nanmean(goRT6)]
[nanmean(stopRT1) nanmean(stopRT2) nanmean(stopRT3) nanmean(stopRT4) nanmean(stopRT5) nanmean(stopRT6)]

%%
RT3 = nan(length(unique(allData.sess)),1);
RT6 = nan(length(unique(allData.sess)),1);
for i = 1 : length(RT3);
    stop3Trial = allData.sess == i & allData.cnd == 3 & allData.ssd ~= 0 & allData.stm1 == 1 & allData.resp == 1;
RT3(i) = nanmean(allData.rt(stop3Trial));
    stop3Trial = allData.sess == i & allData.cnd == 1 & allData.ssd ~= 0 & allData.stm1 == 2 & allData.resp == 2;
RT6(i) = nanmean(allData.rt(stop3Trial));
end



%%
for iSsd = 1 : 5
    iSsd
    [obs.nStopICorr(find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('stopTrial.*ssd%d_.*',iSsd), 'once')),obs.trialCat,'Uni',0)))),...
        obs.nStopIErrorCCorr(find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('stopTrial.*ssd%d_.*',iSsd), 'once')),obs.trialCat,'Uni',0))))]
    
end
%%

load('~/perceptualchoice_stop_model/data/2016-11-04/preproc01/subj02/data_subj02.mat')
data = data(data.ssd > 0, :);
    for cnd = 1:3
        
        rt1(cnd) = nanmean(data.rt(data.resp == 1 & data.stm1 == 1 & data.cnd == cnd))
        rt2(cnd) = nanmean(data.rt(data.resp == 2 & data.stm1 == 2 & data.cnd == cnd))
    end

    %%
 load('~/perceptualchoice_stop_model/data/2016-11-04/preproc01/subj02/data_subj02.mat')
data = data(data.ssd == 0, :);
       for cnd = 1:3
        
        rt1(cnd) = nanmean(data.rt(data.resp == 1 & data.stm1 == 1 & data.cnd == cnd))
        rt2(cnd) = nanmean(data.rt(data.resp == 2 & data.stm1 == 2 & data.cnd == cnd))
    end
%%

load('~/perceptualchoice_stop_model/data/2016-11-04/preproc01/subj01/data_subj01.mat')
data = data(data.ssd > 0, :);
    for cnd = 1:3
        
        rt1(cnd) = nanmean(data.rt(data.resp == 1 & data.stm1 == 1 & data.cnd == cnd))
        rt2(cnd) = nanmean(data.rt(data.resp == 2 & data.stm1 == 2 & data.cnd == cnd))
    end
%%

opt = ccm_options;
opt.ssd = 'collapse';

colorCoh = unique(trialData.targ1CheckerProp);
for i = 1 : length(unique(trialData.targ1CheckerProp))
    iCoh = colorCoh(i);
    opt.rightCheckerPct = iCoh * 100;

opt.outcome = {'stopIncorrectTarget'};
stop1 = ccm_trial_selection(trialData, opt);
pre1(i) = nanmean(trialData.responseOnset(stop1) - trialData.responseCueOn(stop1));

opt.outcome = {'stopIncorrectTarget', 'targetHoldAbort'};
stop2 = ccm_trial_selection(trialData, opt);
pre2(i) = nanmean(trialData.responseOnset(stop2) - trialData.responseCueOn(stop2));
end

%%

opt = ccm_options;
opt.ssd = 'none';

colorCoh = unique(trialData.targ1CheckerProp);
for i = 1 : length(unique(trialData.targ1CheckerProp))
    iCoh = colorCoh(i);
    opt.rightCheckerPct = iCoh * 100;

opt.outcome = {'goCorrectTarget'};
stop1 = ccm_trial_selection(trialData, opt);
pre1C(i) = nanmean(trialData.responseOnset(stop1) - trialData.responseCueOn(stop1));

opt.outcome = {'goCorrectTarget', 'targetHoldAbort'};
stop2 = ccm_trial_selection(trialData, opt);
pre2C(i) = nanmean(trialData.responseOnset(stop2) - trialData.responseCueOn(stop2));


opt.outcome = {'goCorrectDistractor'};
stop1 = ccm_trial_selection(trialData, opt);
pre1E(i) = nanmean(trialData.responseOnset(stop1) - trialData.responseCueOn(stop1));

opt.outcome = {'goCorrectDistractor', 'distractorHoldAbort'};
stop2 = ccm_trial_selection(trialData, opt);
pre2E(i) = nanmean(trialData.responseOnset(stop2) - trialData.responseCueOn(stop2));
end

%%
subject                 = [1 2];
model                   = [2,43,79,352,478];
architecture            = {'crace_ili','cli_ili','cffi_ili'};   
architecture            = {'crace'};   
dt                      = 5;
trialVar                = true;
simScope                = 'go';
fileStr.root            = '/Volumes/HD-1/Users/paulmiddlebrooks/perceptualchoice_stop_model/data/2016-11-04/preproc01/subj%.2d/dt%d/%s/%s/';
fileStr.src             = '/Volumes/HD-1/Users/paulmiddlebrooks/perceptualchoice_stop_model/src/code/2016-11-04/matlab/';
fileStr.result          = '/Volumes/HD-1/Users/paulmiddlebrooks/perceptualchoice_stop_model/results/2016-11-04/subj%.2d/dt%d/%s/%s/';
responseSide            = {'both'};
accuracy                = {'both'};
defective               = true;
savePlot                = true;

cd(fileStr.src)
for iRsp = 1 : length(responseSide)
    respSide = responseSide{iRsp};
    for iAcc = 1 : length(accuracy)
        accur = accuracy{iAcc};
plot_fits(subject,model,architecture,dt,trialVar,simScope,fileStr,respSide,accur,defective,savePlot); 
    end
end

%%  Check Xena's stop Rts
 load('~/perceptualchoice_stop_model/data/2017-05-08/preproc01/subj02/data_subj02.mat')

 for cnd = 1:3
        
        rt1(cnd) = nanmean(data.rt(data.resp == 1 & data.stm1 == 1 & data.cnd == cnd & data.ssd > 0))
        rt2(cnd) = nanmean(data.rt(data.resp == 2 & data.stm1 == 2 & data.cnd == cnd & data.ssd > 0))
       end

%%  Check Xena's stop Rts
 load('~/perceptualchoice_stop_model/data/raw/xena_behavior1.mat')

 opt = ccm_options;
opt.ssd = 'collapse';
trialData.targ1CheckerProp(trialData.targ1CheckerProp == .52) = .53;
colorCoh = unique(trialData.targ1CheckerProp);
for i = 1 : length(unique(trialData.targ1CheckerProp))
    iCoh = colorCoh(i);
    opt.rightCheckerPct = iCoh * 100;

% opt.outcome = {'stopIncorrectTarget'};
% stop1 = ccm_trial_selection(trialData, opt);
% pre1(i) = nanmean(trialData.responseOnset(stop1) - trialData.responseCueOn(stop1));

opt.outcome = {'stopIncorrectTarget', 'targetHoldAbort'};
stop2 = ccm_trial_selection(trialData, opt);
size(stop2)
pre2(i) = nanmean(trialData.responseOnset(stop2) - trialData.responseCueOn(stop2));
end

%%
data = allData;
for cnd = 1:3
        nTrial1(cnd) = sum(data.resp == 1 & data.stm1 == 1 & data.cnd == cnd & data.ssd > 0)
        nTrial2(cnd) = sum(data.resp == 2 & data.stm1 == 2 & data.cnd == cnd & data.ssd > 0)
        rt1(cnd) = nanmean(data.rt(data.resp == 1 & data.stm1 == 1 & data.cnd == cnd & data.ssd > 0))
        rt2(cnd) = nanmean(data.rt(data.resp == 2 & data.stm1 == 2 & data.cnd == cnd & data.ssd > 0))
       end

    
       %%
       for cnd = 1:3
       goCCorr1(cnd) = sum(data.cnd == cnd & data.resp ==1 & data.stm1 == 1 & data.ssd == 0)
       goCCorr2(cnd) = sum(data.cnd == cnd & data.resp ==2 & data.stm1 == 2 & data.ssd == 0)
       goCErr1(cnd) = sum(data.cnd == cnd & data.resp ==2 & data.stm1 == 1 & data.ssd == 0)
       goCErrr2(cnd) = sum(data.cnd == cnd & data.resp ==1 & data.stm1 == 2 & data.ssd == 0)
       stopICCorr1(cnd) = sum(data.cnd == cnd & data.resp ==1 & data.stm1 == 1 & data.ssd > 0)
       stopICCorr2(cnd) = sum(data.cnd == cnd & data.resp ==2 & data.stm1 == 2 & data.ssd > 0)
       stopICErr1(cnd) = sum(data.cnd == cnd & data.resp ==2 & data.stm1 == 1 & data.ssd > 0)
       stopICErr2(cnd) = sum(data.cnd == cnd & data.resp ==1 & data.stm1 == 2 & data.ssd > 0)
       stopCorrect(cnd) = sum(data.cnd == cnd & data.ssd > 0 & data.acc == 1 & isnan(data.resp))
       end

       %%  Check Xena's stop Rts
 load('~/perceptualchoice_stop_model/data/2017-05-08/preproc01/subj01/data_subj01.mat')

 
 
 %%
 
 
 subj                    = [3];
model                   = [79 478];
trialVar                = true;
dt                      = 1;
simScope                = 'all';
%choiceMech              = {'race','li','ffi'};
choiceMech              = {'ffi'};
stopMech                = 'li';
fileStr.root            = strcat(modelRoot,'/data/2017-05-18/preproc01/subj%.2d/dt%d/%s/%s/');
fileStr.fitRoot         = strcat(modelRoot,'/data/2017-05-18/preproc01/subj%.2d/dt%d/%s/%s/');
fileStr.result          = strcat(modelRoot,'/results/2017-05-18/subj%.2d/dt%d/%s/%s/');

    for jChoice = 1 : length(choiceMech)
        p = find_best_params(subj,trialVar,simScope,choiceMech{jChoice},stopMech,dt,model,fileStr)
    end

%%
    
subj                    = [2];
model                   = [478];
trialVar                = true;
dt                      = 1;
simScope                = 'all';
% choiceMech              = {'race','li','ffi'};
choiceMech              = {'ffi'};
stopMech                = 'li';
fileStr.root            = strcat(modelRoot,'/data/2017-05-17/preproc01/subj%.2d/dt%d/%s/%s/');
fileStr.fitRoot         = strcat(modelRoot,'/data/2017-05-17/preproc01/subj%.2d/dt%d/%s/%s/');
fileStr.result          = strcat(modelRoot,'/results/2017-05-17/subj%.2d/dt%d/%s/%s/');

    for jChoice = 1 : length(choiceMech)
        p = find_best_params(subj,trialVar,simScope,choiceMech{jChoice},stopMech,dt,model,fileStr)
    end
    
    
    
    
%%

% Make the activation function a matrix to ease user/coder code reading
kGoStopActAll = cell2mat(prd.dyn{kTrialCatStop}.stopICorr.goStim.targetGO.sY);

% Find maximum of activation functions after SSD
[maxGoStopAct, maxInd] = max(kGoStopActAll(:, iSsdArray(kSSD):end), [], 2);

% Use half the maximum to determine cancel time
halfMaxGoStopAct = maxGoStopAct / 2;

% Determine first index of activation function that falls below half-max.
halfMaxInd = find(kGoStopActAll(:,iSsdArray(kSSD)+maxInd:end) < halfMaxGoStopAct);

kCancelTime = iSsdArray(kSSD) + maxInd + halfMaxInd;


% arrayfun(@(x) find(x(:,iSsdArray(kSSD)+maxInd:end) < halfMaxGoStopAct), kGoStopActAll, 'uni', false)


% For each simulated trial, determine first index of activation function
% that falls below half-max
kCancelTime = nan(length(kGoStopActAll), 1);
for a = 1 : length(kGoStopActAll)
    kCancelTime(a) = find(kGoStopActAll(a,:) < halfMaxGoStopAct, 'first');
end
    
    
%%

subject                 = [1 3];
model                   = [79,352,478];
architecture            = {'crace_ili','cli_ili','cffi_ili'};
dt                      = 1;
trialVar                = true;
simScope                = 'all';
fileStr.root            = '/Volumes/HD-1/Users/paulmiddlebrooks/perceptualchoice_stop_model/data/2017-05-08/preproc01/subj%.2d/dt%d/%s/%s/';
fileStr.src             = '/Volumes/HD-1/Users/paulmiddlebrooks/perceptualchoice_stop_model/src/code/2017-05-08/matlab/';
fileStr.result          = '/Volumes/HD-1/Users/paulmiddlebrooks/perceptualchoice_stop_model/results/2017-05-08/subj%.2d/dt%d/%s/%s/';
fileStr.root            = '~/perceptualchoice_stop_model/data/2017-05-08/preproc01/subj%.2d/dt%d/%s/%s/';
fileStr.src             = '~/perceptualchoice_stop_model/src/code/2017-05-08/matlab/';
fileStr.result          = '~/perceptualchoice_stop_model/results/2017-05-08/subj%.2d/dt%d/%s/%s/';
responseSide            = {'both'};
accuracy                = {'correct'};
savePlot                = true;

if trialVar
    trialVarStr = 'trialvar';
else
    trialVarStr = 'notrialvar';
end

cd(fileStr.src)


% model                   = [79];
% architecture            = {'crace_ili'};


subject = 3;
model                   = 478;
architecture            = {'cffi_ili'};

switch subject
    case 1
        nSim = 14;
    case 3
        nSim = 65;
end



for kArch = 1 : length(architecture)
    for jMod = model
    modelName = ['model',num2str(jMod)];
    
     saveDir             = fullfile(sprintf(fileStr.root,subject,dt,trialVarStr,architecture{kArch}));
        fileName = sprintf('cancel_times_%s_%sSimulations', modelName, num2str(nSim));
        load(fullfile(saveDir, fileName))
        
    end
end
                        
%%


    
    
    
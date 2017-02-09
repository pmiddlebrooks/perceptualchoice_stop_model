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

fileStr.root   = fullfile(modelRoot,'data/2016-05-20/preproc01/subj%.2d/dt%d/%s/%s/');
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
fileStr.root            = strcat(modelRoot,'/data/2016-05-20/preproc01/subj%.2d/dt5/%s/%s/');
doPlot                  = true;
doSave                  = true;
doStartParCluster       = false;

cd(fullfile(modelRoot,'src/code/2016-05-20/matlab'));

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








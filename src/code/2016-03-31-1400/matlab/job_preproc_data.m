function job_preproc_data(rawDataDir,subject,sessions,preprocDataDir)
% JOB_PREPROC_DATA Preprocesses behavioral data from Paul's perceptual
% choice stop-signal task.
%
% DESCRIPTION
% This scripts processes the original task performance files. It involves the following steps:
% - Recoding of data, so that it is compatible with the stochasic accumulator modeling pipeline
% - Trials were removed, based on the following criteria
% 	* Trials in practice blocks
% 	*	Trials on which an omission error was made
% 	*	Trials with response times < 150 ms
%
% SYNTAX
% JOB_PREPROC_DATA(rawDataDir,preprocDataDir)
% rawDataDir        - directory where raw performance text files reside
% preprocDataDir    - directory where preprocessed data will be written
% .....................................................P....................
% Bram Zandbelt, bramzandbelt@gmail.com
% $Created : Wed 12 Mar 2014 11:11:29 CDT by bram
% $Modified: Wed 12 Mar 2014 12:30:20 CDT by bram

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. PROCESS INPUTS AND SPECIFY VARIABLES
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist(preprocDataDir,'dir') ~= 7
    mkdir(preprocDataDir)
end

% Assyming fixed stimulus onset time and duration (go-signal, stop-signal)
stmOns            = [0 NaN];
stmDur            = [1000 100];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. IMPORT RAW DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch subject
    case 'broca'
        switch sessions
            case 'behavior2'
        subjectInd = 1;
            case 'neural3'
        subjectInd = 1;
            otherwise
                error(fprintf('No sessions for %s called % yet', subject, sessions))
        end
    case 'xena'
        subjectInd = 2;
end

load(fullfile(rawDataDir,[subject,'_',sessions,'.mat']));

% Remove SSDs with low trial numbers?
% =========================================================================
nSSDMin         = 100; % minimum # of stop trials for a SSD to be included
truncateSSD     = true;
if truncateSSD
    ssdList = unique(trialData.ssd(~isnan(trialData.ssd), :));
    
    nSSD = nan(length(ssdList), 1);
    for i = 1 : length(ssdList)
        nSSD(i) = sum(trialData.ssd == ssdList(i));
    end
    removeSSD = ssdList(nSSD < nSSDMin);
    Lia = ismember(trialData.ssd, removeSSD);
    trialData = trialData(~Lia, :);
end




nObs = size(trialData,1);

% Pre-allocate arrays
% =========================================================================

allData           = dataset({nan(nObs,1),'subj'}, ...
    {nan(nObs,1),'sess'}, ...
    {nan(nObs,1),'block'}, ...
    {nan(nObs,1),'cnd'}, ...
    {nan(nObs,1),'stm1'}, ...
    {nan(nObs,1),'stm1Ons'}, ...
    {nan(nObs,1),'stm1Dur'}, ...
    {nan(nObs,1),'stm2'}, ...
    {nan(nObs,1),'stm2Ons'}, ...
    {nan(nObs,1),'stm2Dur'}, ...
    {nan(nObs,1),'iSSD'}, ...
    {nan(nObs,1),'ssd'}, ...
    {nan(nObs,1),'rsp1'}, ...
    {nan(nObs,1),'rsp2'}, ...
    {nan(nObs,1),'resp'}, ...
    {nan(nObs,1),'rt'}, ...
    {nan(nObs,1),'acc'});

% This code needs to be adjusted. Idealy, you make one big dataset
% variable for all subjects and all sessions, run preprocessing on it, and
% then split into subject-specific datasets. See below for how I did it for
% my code. This will be slightly different for you, depending on the format
% of your data files.
%
% Here is the code I use:
%
% % Identify task performance files
% % =========================================================================
% file = regexpdir(rawDataDir,'.*.mat$');
% %
% % % Import data
% % % =========================================================================
% for iFile = 1:size(file,1)
%
% %     load /Users/paulmiddlebrooks/Dropbox/PerceptualChoiceStopping/data/raw/broca_behavior2.mat
% %     load /scratch/middlepg/perceptualchoice_stop_model/data/raw/broca_behavior2.mat
%
%     % Extract path, name, extension of file
%     [~, n] = fileparts(file{iFile});
%
%     % Extract subject ID and session ID from file name
%     [id] = sscanf(n,'%s_behavior*');
%
%     % Import into dataset array
%     warning off
%     thisData = dataset('File',file{iFile},'HeaderLines',1,'VarNames', ...
%           {'cnd','block','stm1','stm2','iSSD','acc','resp','rt','ssd'});
%     warning on
% %     thisData.subj = id(1)*ones(size(thisData,1),1);
% %     thisData.sess = thisData.sessionTag;
%
%     % Concatenate dataset arrays
%     allData = [allData;thisData];
% end
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. FILL IN ALLDATA (NOW FOR BROCA ONLY)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Identify trials
% =========================================================================

iGoCCorr          = cell2mat(cellfun(@(in1) ~isempty(regexp(in1,'goCorrectTarget','once')),trialData.trialOutcome,'Uni',0));
iGoCError         = cell2mat(cellfun(@(in1) ~isempty(regexp(in1,'goCorrectDistractor','once')),trialData.trialOutcome,'Uni',0));
iStopICorr        = cell2mat(cellfun(@(in1) ~isempty(regexp(in1,'stopCorrect','once')),trialData.trialOutcome,'Uni',0));
iStopIErrorCCorr  = cell2mat(cellfun(@(in1) ~isempty(regexp(in1,'stopIncorrectTarget','once')),trialData.trialOutcome,'Uni',0));
iStopIErrorCError = cell2mat(cellfun(@(in1) ~isempty(regexp(in1,'stopIncorrectDistractor','once')),trialData.trialOutcome,'Uni',0));

iCorr             = iGoCCorr | iStopICorr;

% Aborted trials
iAbort            = cell2mat(cellfun(@(in1) ~isempty(regexp(in1,'Abort$','once')),trialData.trialOutcome,'Uni',0));

% Non-fixated trials
iNoFix            = cell2mat(cellfun(@(in1) ~isempty(regexp(in1,'^noFixation','once')),trialData.trialOutcome,'Uni',0));

% Omission error trials
iOmit             = cell2mat(cellfun(@(in1) ~isempty(regexp(in1,'goIncorrect','once')),trialData.trialOutcome,'Uni',0));

% Trials with different target angles
iTarget0          = trialData.targAngle == 0;
iTarget180        = trialData.targAngle == 180;

% Subject index
% =========================================================================
% May need to be recoded; now set to 1 for all rows
allData.subj(:)  = subjectInd;

% Session index
% =========================================================================
allData.sess(:)  = trialData.sessionTag;

% Task block index
% =========================================================================
% May need to be recoded; now set to 1 for all rows
allData.block(:) = 1;

% Condition index
% =========================================================================
switch subject
    case 'broca'
% Easy
allData.cnd(trialData.targ1CheckerProp == 0.43 | ...
    trialData.targ1CheckerProp == 0.57)       = 1;
% Intermediate
allData.cnd(trialData.targ1CheckerProp == 0.45 | ...
    trialData.targ1CheckerProp == 0.55)       = 2;
% Hard
allData.cnd(trialData.targ1CheckerProp == 0.47 | ...
    trialData.targ1CheckerProp == 0.53)       = 3;
    case 'xena'
% Easy
allData.cnd(trialData.targ1CheckerProp == 0.35 | ...
    trialData.targ1CheckerProp == 0.65)       = 1;
% Intermediate
allData.cnd(trialData.targ1CheckerProp == 0.42 | ...
    trialData.targ1CheckerProp == 0.58)       = 2;
% Hard
allData.cnd(trialData.targ1CheckerProp == 0.47 | ...
    trialData.targ1CheckerProp == 0.52 | ...
    trialData.targ1CheckerProp == 0.53)       = 3;
    otherwise
        error('Need to enter color coherence values for this subject, job_preproc_data line 193')
end
% Primary stimulus (go-signal) index
% =========================================================================
allData.stm1(trialData.targ1CheckerProp < 0.5 & (iGoCCorr | iStopIErrorCCorr)) = 1;
allData.stm1(trialData.targ1CheckerProp > 0.5 & (iGoCError | iStopIErrorCError)) = 1;
allData.stm1(trialData.targ1CheckerProp > 0.5 & (iGoCCorr | iStopIErrorCCorr)) = 2;
allData.stm1(trialData.targ1CheckerProp < 0.5 & (iGoCError | iStopIErrorCError)) = 2;

% Primary stimulus (go-signal) onset
% =========================================================================
% Assume identical onsets across trials (i.e. stimulus-aligned)
allData.stm1Ons(:)   = stmOns(1);

% Primary stimulus (go-signal) duration
% =========================================================================
% Assume identical durations across trials
allData.stm1Dur(:)   = stmDur(1);

% Secondary stimulus (stop-signal) index
% =========================================================================
allData.stm2(:)      = 0;
allData.stm2(trialData.stopSignalOn > 0) = 1;

% Secondary stimulus (stop-signal) onset
% =========================================================================
% Assume identical onsets across trials (i.e. stimulus-aligned)
allData.stm2Ons(:)   = stmOns(2);

% Secondary stimulus (stop-signal) duration
% =========================================================================
% Assume identical durations across trials
allData.stm2Dur(:)    = stmDur(2);

% Stop-signal delay index
% =========================================================================
ssdLevels             = unique(nonnans(trialData.ssd));
nSsd                  = numel(ssdLevels);

allData.iSSD(:)       = 0;
for iSsd = 1:nSsd
    allData.iSSD(trialData.ssd == ssdLevels(iSsd)) = iSsd;
end

% Stop-signal delay
% =========================================================================
allData.ssd(:)        = trialData.ssd;
allData.ssd(isnan(trialData.ssd)) = 0;

% Target response index primary stimulus
% =========================================================================
% Assume following stimulus-response mapping:
%
% Stimulus   Response
allData.rsp1(:)       = allData.stm1;

% Target response index secondary stimulus
% =========================================================================
% This one is not relevant for this task
allData.rsp2(:)      = allData.stm2;

% Response made: 1 is left, 2 is right
% =========================================================================

allData.resp(iTarget0   & iGoCCorr)             = 2;
allData.resp(iTarget0   & iGoCError)            = 2;

allData.resp(iTarget180 & iGoCCorr)             = 1;
allData.resp(iTarget180 & iGoCError)            = 1;

allData.resp(iTarget0   & iStopIErrorCCorr)     = 2;
allData.resp(iTarget0   & iStopIErrorCError)    = 2;

allData.resp(iTarget180 & iStopIErrorCCorr)     = 1;
allData.resp(iTarget180 & iStopIErrorCError)    = 1;

% Response time
% =========================================================================
allData.rt(:)        = trialData.responseOnset - trialData.responseCueOn;

% Accuracy
% =========================================================================
allData.acc(:)       = 0;
allData.acc(iCorr)   = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. REMOVE OUTLIER TRIALS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Remove aborted and non-fixated trials
% =========================================================================
allData(iAbort | iNoFix | iOmit,:) = [];

% Remove express saccades
% =========================================================================
allData(allData.stm2 == 0 & allData.rt < 100,:) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. SAVE DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subj                = unique(allData.subj);
nSubj               = length(subj);

for iSubj = 1:nSubj
    
    subjDir = fullfile(preprocDataDir,sprintf('subj%.2d',subj(iSubj)));
    if exist(subjDir) ~=7
        mkdir(subjDir)
    end
    
    data = allData(find(allData.subj == subj(iSubj)),:);
    fName = fullfile(subjDir,sprintf('data_subj%.2d.mat',subj(iSubj)));
    save(fName, 'data');
    
end
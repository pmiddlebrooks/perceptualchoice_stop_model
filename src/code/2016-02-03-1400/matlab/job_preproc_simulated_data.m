function data = job_preproc_simulated_data(prd, subjectInd, preprocDataDir)
% JOB_PREPROC_SIMULATED_DATA Reformats predicted data from SAM structure
% into preprocessed fromat (for parameter recovery testing)
%
% DESCRIPTION
% This scripts processes simulated predicted (prd) data.
% It involves the following steps:
% - Recoding of data, so that it is compatible with the stochasic accumulator modeling pipeline
% - Trials were removed, based on the following criteria
% 	* Trials in practice blocks
% 	*	Trials on which an omission error was made
% 	*	Trials with response times < 150 ms
%
% SYNTAX
% simData = JOB_PREPROC_DATA(prd)
% prd        - predicted data from simulated data based on given parameters
% simData    - simulated data reformated into preprocessed format (see
% job_preprocess_data.m)
% ..........................................................................


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. PROCESS INPUTS AND SPECIFY VARIABLES
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Assuming fixed stimulus onset time and duration (go-signal, stop-signal)
stmOns            = [0 NaN];
stmDur            = [1000 100];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. FILL IN ALLDATA (NOW FOR BROCA ONLY)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nObs = sum(prd.nTotal(:));


% Pre-allocate arrays
% =========================================================================

% Create an empty dataset allData to build iterating over factors
VarNames = {'subj','sess','block','cnd','stm1','stm1Ons','stm1Dur','stm2','stm2Ons','stm2Dur','iSSD','ssd','rsp1','rsp2','resp','rt','acc'};
allData = dataset([],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],[],'VarNames',VarNames);


% Loop through the factors to build out the simulated preprocessed data
% structure
% =========================================================================
nTrialCat = size(prd, 1);

for iTrialCat = 1 : nTrialCat
    
    iTotal = prd.nTotal(iTrialCat);
    
    % Preallocate a dataset to add to allData for each factor
    trialCatData           = dataset({nan(iTotal,1),'subj'}, ...
        {nan(iTotal,1),'sess'}, ...
        {nan(iTotal,1),'block'}, ...
        {nan(iTotal,1),'cnd'}, ...
        {nan(iTotal,1),'stm1'}, ...
        {nan(iTotal,1),'stm1Ons'}, ...
        {nan(iTotal,1),'stm1Dur'}, ...
        {nan(iTotal,1),'stm2'}, ...
        {nan(iTotal,1),'stm2Ons'}, ...
        {nan(iTotal,1),'stm2Dur'}, ...
        {nan(iTotal,1),'iSSD'}, ...
        {nan(iTotal,1),'ssd'}, ...
        {nan(iTotal,1),'rsp1'}, ...
        {nan(iTotal,1),'rsp2'}, ...
        {nan(iTotal,1),'resp'}, ...
        {nan(iTotal,1),'rt'}, ...
        {nan(iTotal,1),'acc'});
    
    
    % Subject index
    % =========================================================================
    % May need to be recoded; now set to 1 for all rows
    trialCatData.subj(:)  = subjectInd;
    
    % Session index
    % =========================================================================
    % Since this is simulated data, count all as a single session
    trialCatData.sess(:)  = 1;
    
    % Task block index
    % =========================================================================
    % May need to be recoded; now set to 1 for all rows
    trialCatData.block(:) = 1;
    
    
    % Is it a go or stop trial?
    if regexp(prd.trialCat{iTrialCat}, 'goTrial')
        goStop = 'go';
    elseif regexp(prd.trialCat{iTrialCat}, 'stopTrial')
        goStop = 'stop';
    end
    
    % Condition index
    % =========================================================================
    cnd = str2double(regexp(prd.trialCat{iTrialCat},'(?<=_c)\d+','match','once')); % take the number after the ",c" in the tag
    trialCatData.cnd(:) = cnd;
    
    % Primary stimulus (go-signal) index
    % =========================================================================
    stm1 = str2double(regexp(prd.trialCat{iTrialCat},'(?<=r)\d+','match','once')); % take the number after the ",c" in the tag
    trialCatData.stm1(:) = stm1;
    
    % Primary stimulus (go-signal) onset
    % =========================================================================
    % Assume identical onsets across trials (i.e. stimulus-aligned)
    trialCatData.stm1Ons(:)   = stmOns(1);
    
    % Primary stimulus (go-signal) duration
    % =========================================================================
    % Assume identical durations across trials
    trialCatData.stm1Dur(:)   = stmDur(1);
    
    % Secondary stimulus (stop-signal) index
    % =========================================================================
    switch goStop
        case 'go'
            trialCatData.stm2(:)      = 0;
        case 'stop'
            trialCatData.stm2(:) = 1;
    end
    
    % Secondary stimulus (stop-signal) onset
    % =========================================================================
    % Assume identical onsets across trials (i.e. stimulus-aligned)
    trialCatData.stm2Ons(:)   = stmOns(2);
    
    % Secondary stimulus (stop-signal) duration
    % =========================================================================
    % Assume identical durations across trials
    trialCatData.stm2Dur(:)    = stmDur(2);
    
    % Stop-signal delay value and index
    % =========================================================================
    switch goStop
        case 'go'
            trialCatData.iSSD(:) = 0;
            trialCatData.ssd(:) = 0;
        case 'stop'
            iSSD = cell2num(regexp(prd.trialCat{iTrialCat},'(?<=stopTrial_ssd)\d+','match'));
            trialCatData.iSSD(:) = iSSD;
            trialCatData.ssd(:) = prd.ssd(iTrialCat);
    end
    
    % Target response index primary stimulus
    % =========================================================================
    % Assume following stimulus-response mapping:
    %
    % Stimulus   Response
    trialCatData.rsp1       = trialCatData.stm1;
    
    % Target response index secondary stimulus
    % =========================================================================
    % This one is not relevant for this task
    trialCatData.rsp2      = trialCatData.stm2;
    
    
    % Response made, Response time, and Accuracy
    % =========================================================================
    if stm1 == 1
        cResp = 1;
        eResp = 2;
    else
        cResp = 2;
        eResp = 1;
    end

    trialCatData.resp(1 : prd.nGoCCorr(iTrialCat)) = cResp;
    trialCatData.resp(prd.nGoCCorr(iTrialCat)+1 : prd.nGoCCorr(iTrialCat) + prd.nGoCError(iTrialCat)) = eResp;
    
    trialCatData.rt(1 : prd.nGoCCorr(iTrialCat)) = prd.rtGoCCorr{iTrialCat};
    trialCatData.rt(prd.nGoCCorr(iTrialCat)+1 : prd.nGoCCorr(iTrialCat) + prd.nGoCError(iTrialCat)) = prd.rtGoCError{iTrialCat};
    
    trialCatData.acc(1 : prd.nGoCCorr(iTrialCat)) = 1;
    trialCatData.acc(prd.nGoCCorr(iTrialCat)+1 : prd.nGoCCorr(iTrialCat) + prd.nGoCError(iTrialCat)) = 0;
    
    
    % Concatenate the trialCatData with allData
    
    allData = [allData; trialCatData];
    
end % iTrialCat



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. REMOVE OUTLIER TRIALS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Remove express saccades
% =========================================================================
allData(allData.stm2 == 0 & allData.rt < 100,:) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5. SAVE DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    subjDir = fullfile(preprocDataDir,sprintf('subj%.2d',subjectInd));
    if exist(subjDir) ~=7
        mkdir(subjDir)
    end
    
    data = allData;
    fName = fullfile(subjDir,sprintf('data_subj%.2d.mat',subjectInd));
    save(fName, 'data');
    

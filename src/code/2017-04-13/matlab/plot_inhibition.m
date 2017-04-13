function plot_inhibition(subject,model,architecture,dt,trialVar,optimScope,fileStr,responseSide, accuracy, savePlot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. PROCESS INPUTS & SPECIFY VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
nSim = 1000;
ms2Std = 50;
% 1.1. Process inputs
% =========================================================================

if trialVar
    trialVarStr = 'trialvar';
else
    trialVarStr = 'notrialvar';
end

% 1.2. Specify dynamic variables
% =========================================================================

% Numbers
nSubject            = numel(subject);
nModel              = numel(model);
nArchitecture       = numel(architecture);

% 1.3. Specify static variables
% =========================================================================
rootDir             = fileStr.root;
nameFVal            = 'allFValAndBestX_%sTrials_model%.3d.txt';
modelStr            = {'zc','t0','v/ve','zc & t0','zc & v/ve','t0 & v/ve','zc, t0 & v/ve'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Plot the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop over subjects, architectures, and models
warning('off','stats:dataset:ModifiedVarnames')

% Set up the figure and panels
%         [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_figure;
standard_figure(1,1,'landscape');

%                 set_figure({1024,574,'pixels'},{'USLetter','landscape'},{'Helvetica',18});
p = panel();
p.pack({.1 .25 .25 .25 .15}, num2cell(repmat(1/6,1,6)));

annotation('textbox', [0 0.9 1 0.1], ...
    'String', sprintf('Subject %d, architecture %s',subject,architecture), ...
    'Interpreter','none', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');


% Display progress
fprintf(1,'Working on subject %d, architecture %s, and model %d. \n',subject,architecture,model);

% Load model-specific fits with cost function values

ds = dataset('File',fullfile(sprintf(rootDir,subject,dt,trialVarStr,architecture),sprintf(nameFVal,optimScope,model)));

% Load SAM
load(fullfile(sprintf(rootDir,subject,dt,trialVarStr,architecture),ds.FileNameSAM{1}),'SAM');

% Extract optimized X
iBestX = cell2mat(cellfun(@(in1) ~isempty(regexp(in1,'^BestX.*', 'once')),ds.Properties.VarNames,'Uni',0));
X = double(ds(1,iBestX));

% Alter the simulations for efficiency
SAM.sim.n = nSim;

% Specificy observations
obs = SAM.optim.obs;


% SSD index
nSSD = unique(obs.ssd(~isnan(obs.ssd)));
switch subject
    case 1
        iSsd = ceil(length(nSSD) * .65);
        iSsd = 4;
    case 2
        iSsd = 7;
end

% Get model predictions and costs
% [cost,altCost,prd] = sam_cost(X,SAM);
prd = sam_sim_expt('explore',X,SAM);

% Get rid of conditions for which ssd = 0 (an artifact of preporcessing)
obs(prd.ssd == 0 ,:) = [];
prd(prd.ssd == 0 ,:) = [];


% Plot observations and predictions
% =========================================================================
%                     sam_plot(SAM,prd);

% Plot it
plotit(SAM,prd,p,model,modelStr, responseSide, accuracy);

if savePlot
    saveDir             = fullfile(sprintf(fileStr.result,subject,dt,trialVarStr,architecture));
    if exist(saveDir,'dir') ~= 7
        mkdir(saveDir)
    end
    fileName = sprintf('Respond_%s_Accuracy_%s', responseSide, accuracy);
    print(gcf, fullfile(saveDir, fileName),'-dpdf', '-r300')
    print(gcf, fullfile(saveDir, fileName),'-dpng')
end
clear SAM prd















    function plotit(SAM,prd,p,model,modelStr, responseSide, accuracy)
        
        % Specify colors and line properties
        mapShades           	= [.1 .25 .4 .6 .75 .9];
        cndClr                  = ccm_colormap(mapShades);
        %         cndClr(1+length(mapShades)/2 : end, :) = flipud(cndClr(1+length(mapShades)/2 : end, :));
        
        colorCohArray = [1 2 3 6 5 4];
        
        conditionArray = 1:3;
        switch responseSide
            case 'both'
                responseArray = [1 2];
            case 'left'
                responseArray = 1;
            case 'right'
                responseArray = 2;
            otherwise
                error('responseSide need to be left, right, or both');
        end
        ssdArray = unique(prd.ssd(~isnan(prd.ssd)));
        ssdArray(ssdArray == 0) = [];
        
        
        
        goCCorrMrkObs           = 'o';
        goCCorrLnObs            = 'none';
        goCCorrMrkPrd           = 'none';
        goCCorrLnPrd            = '-';
        goCCorrLnWidth          = 2;
        
        goCErrorMrkObs          = '^';
        goCErrorLnObs           = 'none';
        goCErrorMrkPrd          = 'none';
        goCErrorLnPrd           = '-.';
        goCErrorLnWidth         = 1;
        
        stopIErrorCCorrMrkObs   = 's';
        stopIErrorCCorrLnObs    = 'none';
        stopIErrorCCorrMrkPrd   = 'none';
        stopIErrorCCorrLnPrd    = '--';
        stopIErrorCCorrLnWidth  = 2;
        
        stopICorrMrkObs         = 'o';
        stopICorrLnObs          = 'none';
        
        stopICorrMrkPrd         = 'none';
        stopICorrLnPrd          = '-';
        stopICorrLnWidth        = 2;
        
        
        % Get Model SSRT for whole session
        stopCat = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('stopTrial.*ssd.*'), 'once')),prd.trialCat,'Uni',0)));
        ssrtModel = mean(cell2mat(cellfun(@(x) reshape(x, [], 1), prd.rtStopICorr(stopCat), 'uni', false)));
        
        
        cancelTime = nan(size(prd, 1), 1);
        
        
        
        
        
        
        
        
        
        for iResponse = 1:length(responseArray)
            for iCondition = 1:length(conditionArray);
                
                iRsp = responseArray(iResponse);
                iCnd = conditionArray(iCondition);
                iCohInd = length(conditionArray) * (iRsp -1) + iCondition; % Which color coherence index is this?
                
                iColor = cndClr(iCohInd, :);
                
                
                % Identify the relevant rows in the dataset array
                iTrialCatGo = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('goTrial.*c%d.*GO.*r%d',iCnd,iRsp), 'once')),prd.trialCat,'Uni',0)));
                if isempty(iTrialCatGo)
                    iTrialCatGo = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('goTrial_c%d.*',iCnd), 'once')),prd.trialCat,'Uni',0)));
                end
                iGoRT = [obs.rtGoCCorr{iTrialCatGo}; obs.rtGoCError{iTrialCatGo}];
                
                iTrialCatStop = cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('stopTrial.*ssd.*c%d.*GO.*r%d',iCnd,iRsp), 'once')),prd.trialCat,'Uni',0));
                iSsdArray = obs.ssd(iTrialCatStop);
                iStopProbRespond = 1 - obs.pStopICorr(iTrialCatStop);
                iNStop = obs.nStopICorr(iTrialCatStop) + obs.nStopIErrorCCorr(iTrialCatStop) + obs.nStopIErrorCError(iTrialCatStop);
                
                
                
                % Identify the relevent indices of parameters in X
                xNames = SAM.model.variants.tree(model).XSpec.name.(optimScope).name;
                [tF, zcIndGOCorr] = ismember(sprintf('zc_{GO}'),xNames);
                [tF, t0Ind] = ismember(sprintf('t0_{GO}'),xNames);
                [tF, zcIndSTOP] = ismember(sprintf('zc_{STOP}'),xNames);
                [tF, z0IndSTOP] = ismember(sprintf('z0_{STOP}'),xNames);
                
                switch model
                    case 79
                        [tF, z0IndGOCorr] = ismember(sprintf('z0_{GO,r%i}', iRsp),xNames);
                    case 352
                        [tF, z0IndGOCorr] = ismember(sprintf('z0_{GO}'),xNames);
                        z0IndGOError = z0IndGOCorr;
                    case 478
                        [tF, z0IndGOCorr] = ismember(sprintf('z0_{GO,r%i}', iRsp),xNames);
                    otherwise
                        error(fprintf('Don''t have code yet for model %i\n', model));
                end
                
                % Get Baseline activity (i.e. standard devation to calculate cancel
                % time
                
                % Get latency-matched GO trials to compare with canceled STOP
                % trials
                [tF, zcIndGOCorr] = ismember(sprintf('zc_{GO}'),xNames);
                
                
                
                
                
                
                
                % SSRT w.r.t. SSD
                % ================================================
                p(2,colorCohArray(iCohInd)).select();
                p(2,colorCohArray(iCohInd)).hold('on');
                %                 p(2,colorCohArray(iCohInd)).title({sprintf('%d',round(cost))});
                
                
                % Observed SSRT estimate
                [fitParameters, lowestSSE] = Weibull(iSsdArray, iStopProbRespond, iNStop);
                ssrt = get_ssrt(iSsdArray, iStopProbRespond, iNStop, iGoRT, fitParameters);
                
                plot(iSsdArray, ssrt.integration,'Color',iColor,'Marker',stopICorrMrkObs,'LineStyle',stopICorrLnObs,'LineWidth',stopICorrLnWidth)
                
                % Predicted SSRTs
                meanSsrtPrd = nan(length(iSsdArray), 1);
                for kSSD = 1 : length(iSsdArray)
                    
                    kTrialCatStop = (prd.ssd == iSsdArray(kSSD)) & iTrialCatStop;
                    if sum(kTrialCatStop)
                        % Model SSRTs
                        meanSsrtPrd(kSSD) = mean(prd.rtStopICorr{kTrialCatStop});
                        %                         plot(ssdArray(kSSD), mean(prd.rtStopICorr{kTrialCatStop}),'Color','k','Marker',stopICorrMrkPrd,'LineStyle',stopICorrLnPrd,'LineWidth',stopICorrLnWidth)
                        %                 rtStopICorrPrdL{iCnd}         = cell2mat(cellfun(@(x) reshape(x, [], 1), prd.rtStopICorr(iTrialCatStop), 'uni', false));
                        
                        
                    end
                    
                    
                end
                plot(iSsdArray, meanSsrtPrd,'Color','k','Marker',stopICorrMrkPrd,'LineStyle',stopICorrLnPrd,'LineWidth',stopICorrLnWidth)
                
                
                
                
                
                
                
                
                
                
                % Model Cancel Times
                % ================================================
                
                
                % goRT indices
                goIndPrd = cellfun(@(in1) find(in1 > X(zcIndGOCorr), 1), prd.dyn{iTrialCatGo}.goCCorr.goStim.targetGO.sY, 'uni', false);
                goRTPrd = cellfun(@(in1, in2) in1(in2), prd.dyn{iTrialCatGo}.goCCorr.goStim.targetGO.sX, goIndPrd);
                % Predicted Cancel Times
                for kSSD = 1 : length(iSsdArray)
                    kTrialCatStop = (prd.ssd == iSsdArray(kSSD)) & iTrialCatStop;
                    
                    % Slow Go RTs
                    kGoSlowRTInd = goRTPrd >  prd.ssd(kTrialCatStop) + ssrtModel;
                    
                    % Only calculate cancel time if there are more than one
                    % trial to do so
                    if sum(kGoSlowRTInd) > 1
                        kGoSlowAct = mean(cell2mat(prd.dyn{iTrialCatGo}.goCCorr.goStim.targetGO.sY(kGoSlowRTInd)), 1)';
                        % Set values (NaNs) beyond RT to the threshold, to
                        % emulate ongoing spiking activity
                        %                     kGoSlowAct(isnan(kGoSlowAct)) = X(zcIndGOCorr);
                        kGoSlowAct(isnan(kGoSlowAct)) = max(kGoSlowAct);
                        kGoStopAct = mean(cell2mat(prd.dyn{kTrialCatStop}.stopICorr.goStim.targetGO.sY), 1)';
                        
                        kActDiff = kGoSlowAct - kGoStopAct;
                        kStd = std(kActDiff(1:prd.ssd(kTrialCatStop)));
                        
                        kStd2Ind = kActDiff > 2*kStd;
                        
                        
                        
                        % Look for a sequence of ms2Std ms for which the go sdf is 2
                        % std greater than the stop sdf.
                        % First whether the differential sdf was > 2*Std for the
                        % first ms2Std ms
                        if sum(kStd2Ind(1:ms2Std)) == ms2Std
                            cancelTime(kTrialCatStop) = 1;
                        else
                            % If it wasn't, determein whether there was a time
                            % after the checkerboard onset that the differential
                            % sdf was > 2*Std for at least ms2Std ms.
                            riseAbove2Std = find([0; diff(kStd2Ind)] == 1);
                            sinkBelow2Std = find([0; diff(kStd2Ind)] == -1);
                            if ~isempty(riseAbove2Std)
                                % Get rid of occasions for which the signals differ
                                % going into the epoch (and therefore they will
                                % cease to differ before they begin again to
                                % differ)
                                if ~isempty(sinkBelow2Std)
                                    sinkBelow2Std(sinkBelow2Std < riseAbove2Std(1)) = [];
                                end
                                
                                % If there's one more riseAbove2Std than sinkBelow2Std, the last riseAbove2Std
                                % will last until the end of the sdf: Add to
                                % singkBelowStd the end of the epoch
                                if length(riseAbove2Std) > length(sinkBelow2Std)
                                    sinkBelow2Std = [sinkBelow2Std; length(kActDiff)];
                                end
                                
                                % Now riseAbove2Std length should be equal. See if
                                % any of the riseAbove2Std streaks go longer than
                                % 50ms
                                ind = find(sinkBelow2Std - riseAbove2Std >= ms2Std, 1);
                                if ~isempty(ind)
                                    cancelTime(kTrialCatStop) = riseAbove2Std(ind);
                                end
                                % If they're not equal, the last riseAbove2Std
                                % will last until the end of the sdf: Add to
                                % singkBelowStd the end of the epoch
                            end
                        end
                    end
                    
                    % If there werwe ms2Std consecutive ms of sdfDiff > 2*Std,
                    % check whether the difference ever reached 6*Std
                    if ~isnan(cancelTime(kTrialCatStop))
                        
                        
                        % If it never reached 6Std, reset cancel time with NaN again
                        std6Ind = kActDiff(cancelTime(kTrialCatStop) : end) > 6*kStd;
                        if ~sum(std6Ind)
                            cancelTime(kTrialCatStop) = nan;
                        end
                    end
                    
                    
                    
                    
                    
                end
                % SSRT w.r.t. SSD
                % ================================================
                p(3,colorCohArray(iCohInd)).select();
                p(3,colorCohArray(iCohInd)).hold('on');
                %                 p(2,colorCohArray(iCohInd)).title({sprintf('%d',round(cost))});
                plot(iSsdArray, cancelTime(iTrialCatStop) - prd.ssd(iTrialCatStop) - ssrtModel,'Color','k','Marker',stopICorrMrkPrd,'LineStyle',stopICorrLnPrd,'LineWidth',stopICorrLnWidth)
                
                
                
                
                
                
                
                
                % SSD distributions
                % ================================================
                
                
                % ================================================
                p(4,colorCohArray(iCohInd)).select();
                p(4,colorCohArray(iCohInd)).hold('on');
                %                 p(2,colorCohArray(iCohInd)).title({sprintf('%d',round(cost))});
                plot(iSsdArray, obs.nTotal(iTrialCatStop)/max(obs.nTotal(iTrialCatStop)),'Color','k','Marker',stopICorrMrkObs,'LineStyle',stopICorrLnObs,'LineWidth',stopICorrLnWidth)
                
                
                
                 
                
                
            end  % iCondition
            
        end % iResponse
        
        
        
        
         
    end

end
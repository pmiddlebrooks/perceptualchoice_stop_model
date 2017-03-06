function plot_inhibition(subject,model,architecture,dt,trialVar,optimScope,fileStr,responseSide, accuracy, savePlot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. PROCESS INPUTS & SPECIFY VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
nSim = 50;
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
[cost,altCost,prd] = sam_cost(X,SAM);

% Get rid of conditions for which ssd = 0 (an artifact of preporcessing)
obs(prd.ssd == 0 ,:) = [];
prd(prd.ssd == 0 ,:) = [];


% Plot observations and predictions
% =========================================================================
%                     sam_plot(SAM,prd);

% Plot it
plotit(SAM,prd,p,model,modelStr,cost,altCost, responseSide, accuracy);

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















    function plotit(SAM,prd,p,model,modelStr,cost,altCost, responseSide, accuracy)
        
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
                

                

                % SSRT w.r.t. SSD
                % ================================================
                p(2,colorCohArray(iCohInd)).select();
                p(2,colorCohArray(iCohInd)).hold('on');
                p(2,colorCohArray(iCohInd)).title({sprintf('%d',round(cost))});
                
                
                % Observed SSRT estimate                
                        [fitParameters, lowestSSE] = Weibull(iSsdArray, iStopProbRespond, iNStop);
            ssrt = get_ssrt(iSsdArray, iStopProbRespond, iNStop, iGoRT, fitParameters);

                        plot(iSsdArray, ssrt.integration,'Color',iColor,'Marker',stopICorrMrkObs,'LineStyle',stopICorrLnObs,'LineWidth',stopICorrLnWidth)
                
                % Predicted SSRTs
                for kSSD = 1 : length(ssdArray)
                    
                    kTrialCatStop = (prd.ssd == ssdArray(kSSD)) & iTrialCatStop;
                    if sum(kTrialCatStop)
                        % Model SSRTs
                        
                        plot(ssdArray(kSSD), mean(prd.rtStopICorr{kTrialCatStop}),'Color','k','Marker',stopICorrMrkPrd,'LineStyle',stopICorrLnPrd,'LineWidth',stopICorrLnWidth)
                        %                 rtStopICorrPrdL{iCnd}         = cell2mat(cellfun(@(x) reshape(x, [], 1), prd.rtStopICorr(iTrialCatStop), 'uni', false));
                        
                        
                    end
                    
                    
                end
                
                
                
                
            end
        end
                
                
                
                % GoCCorr trials
                % -----------------------------------------------------------------
                if strcmp(accuracy, 'both') || strcmp(accuracy, 'correct')
                    rtGoCCorrObs     = obs.rtQGoCCorr{iTrialCatGo};
                    cumPGoCCorrObs   = obs.cumProbGoCCorr{iTrialCatGo};
                    
                    rtGoCCorrPrd     = prd.rtGoCCorr{iTrialCatGo};
                    cumPGoCCorrPrd   = cmtb_edf(prd.rtGoCCorr{iTrialCatGo}(:),prd.rtGoCCorr{iTrialCatGo}(:));
                    
                    
                    % Plot it
                    plot(rtGoCCorrObs,cumPGoCCorrObs,'Color',iColor,'Marker',goCCorrMrkObs,'LineStyle',goCCorrLnObs,'LineWidth',goCCorrLnWidth);
                    plot(rtGoCCorrPrd,cumPGoCCorrPrd,'Color',iColor,'Marker',goCCorrMrkPrd,'LineStyle',goCCorrLnPrd,'LineWidth',goCCorrLnWidth);
                    
                end
                
                % GoCError trials
                % -----------------------------------------------------------------
                if ~isempty(obs.rtQGoCError{iTrialCatGo}) && obs.nGoCError(iTrialCatGo) > 10 && (strcmp(accuracy, 'both') || strcmp(accuracy, 'error'))
                    
                    rtGoCErrorObs     = obs.rtQGoCError{iTrialCatGo};
                    cumPGoCErrorObs   = obs.cumProbGoCError{iTrialCatGo};
                    
                    rtGoCErrorPrd     = prd.rtGoCError{iTrialCatGo};
                    cumPGoCErrorPrd   = cmtb_edf(prd.rtGoCError{iTrialCatGo}(:),prd.rtGoCError{iTrialCatGo}(:));
                    
                    
                    % Plot it
                    iColor = cndClr(iCohInd, :);
                    plot(rtGoCErrorObs,cumPGoCErrorObs,'Color',iColor,'Marker',goCErrorMrkObs,'LineStyle',goCErrorLnObs,'LineWidth',goCErrorLnWidth);
                    plot(rtGoCErrorPrd,cumPGoCErrorPrd,'Color',iColor,'Marker',goCErrorMrkPrd,'LineStyle',goCErrorLnPrd,'LineWidth',goCErrorLnWidth);
                    
                end
                
                % Set axes
                switch subject
                    case 1
                        set(gca,'XLim',[200 700], ...
                            'XTick',100:100:700, ...
                            'YLim',[0 1], ...
                            'YTick',0:0.2:1)
                    case 2
                        set(gca,'XLim',[200 500], ...
                            'XTick',100:100:600, ...
                            'YLim',[0 1], ...
                            'YTick',0:0.2:1)
                    otherwise
                        error('Need to add axes limits for subject')
                end
                
                
                
                % RT distribution Stop trials
                %         ================================================
                p(3,colorCohArray(iCohInd)).select();
                p(3,colorCohArray(iCohInd)).hold('on');
                
                % Get rid of NaNs in obs.cumProbStopIErrorCCorr
                for j = 1 : size(obs, 1)
                    if isnan(obs.cumProbStopIErrorCCorr{j})
                        obs.cumProbStopIErrorCCorr{j} = [];
                    end
                end
                
                iRsp = responseArray(iResponse);
                iCnd = conditionArray(iCondition);
                
                % Identify the relevant rows in the dataset array
                iTrialCatStop = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('stopTrial.*ssd%d.*GO.*r%d.*c%d',iSsd,iRsp,iCnd), 'once')),prd.trialCat,'Uni',0)));
                
                if isempty(iTrialCatStop)
                    iTrialCatStop = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('stopTrial.*ssd%d.*c%d.*r%d.*',iSsd,iCnd,iRsp), 'once')),prd.trialCat,'Uni',0)));
                end
                
                % StopIErrorCCorr trials
                % -----------------------------------------------------------------
                rtStopIErrorCCorrObs   = obs.rtQStopIErrorCCorr{iTrialCatStop};
                cumPStopIErrorCCorrObs = obs.cumProbStopIErrorCCorr{iTrialCatStop};
                
                rtStopIErrorCCorrPrd     = prd.rtStopIErrorCCorr{iTrialCatStop};
                cumPStopIErrorCCorrPrd   = cmtb_edf(prd.rtStopIErrorCCorr{iTrialCatStop}(:),prd.rtStopIErrorCCorr{iTrialCatStop}(:));
                
                
                % Plot it
                iColor = cndClr(iCohInd, :);
                plot(rtStopIErrorCCorrObs,cumPStopIErrorCCorrObs,'Color',iColor,'Marker',stopIErrorCCorrMrkObs,'LineStyle',stopIErrorCCorrLnObs,'LineWidth',stopIErrorCCorrLnWidth);
                plot(rtStopIErrorCCorrPrd,cumPStopIErrorCCorrPrd,'Color',iColor,'Marker',stopIErrorCCorrMrkPrd,'LineStyle',stopIErrorCCorrLnPrd,'LineWidth',stopIErrorCCorrLnWidth);
                
                % StopICorr trials
                % -----------------------------------------------------------------
                rtStopICorrPrd     = prd.rtStopICorr{iTrialCatStop};
                cumPStopICorrPrd   = cmtb_edf(prd.rtStopICorr{iTrialCatStop}(:),prd.rtStopICorr{iTrialCatStop}(:));
                
                % Plot it
                iColor = 'r';%cndClr((length(conditionArray) * (iRsp -1) + iCondition), :);
                plot(rtStopICorrPrd,cumPStopICorrPrd,'Color',iColor,'Marker',stopICorrMrkPrd,'LineStyle',stopICorrLnPrd,'LineWidth',stopICorrLnWidth);
                
                
                % Set axes
                switch subject
                    case 1
                        set(gca,'XLim',[200 700], ...
                            'XTick',100:100:700, ...
                            'YLim',[0 1], ...
                            'YTick',0:0.2:1)
                    case 2
                        set(gca,'XLim',[200 500], ...
                            'XTick',100:100:600, ...
                            'YLim',[0 1], ...
                            'YTick',0:0.2:1)
                    otherwise
                        error('Need to add axes limits for subject')
                end
                
                
                
                % Inhibition function
                % ================================================
                p(4,colorCohArray(iCohInd)).select();
                p(4,colorCohArray(iCohInd)).hold('on');
                %     p(2,colorCohArray(iCohInd)).title({sprintf('Model %d',model(colorCohArray(iCohInd))), ...
                %                        modelStr{colorCohArray(iCohInd)}, ...
                %                        sprintf('\\chi^2 = %.1f',cost), ...
                %                        sprintf('BIC = %.1f',altCost)});
                
                
                iRsp = responseArray(iResponse);
                iCnd = conditionArray(iCondition);
                
                % Identify the relevant rows in the dataset array
                %                     iTrialCatStop = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('stopTrial.*ssd%d.*GO.*r%d.*c%d',iSsd,iRsp,iCnd), 'once')),prd.trialCat,'Uni',0)));
                %
                %                     if isempty(iTrialCatStop)
                %                         iTrialCatStop = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('stopTrial.*ssd%d.*c%d.*r%d.*',iSsd,iCnd,iRsp), 'once')),prd.trialCat,'Uni',0)));
                %                     end
                %                     iTrialCat = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('stopTrial.*GO:c%d',iCnd), 'once')),prd.trialCat,'Uni',0)));
                iTrialCatStop = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('stopTrial.*GO.*r%d.*c%d',iRsp,iCnd), 'once')),prd.trialCat,'Uni',0)));
                if isempty(iTrialCatStop)
                    iTrialCatStop = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('stopTrial.*c%d.*GO.*r%d',iCnd,iRsp), 'once')),prd.trialCat,'Uni',0)));
                end
                % Get the SSDs and observed and predicted response rates
                ssd = obs.ssd(iTrialCatStop);
                respRateObs = 1-obs.pStopICorr(iTrialCatStop);
                respRatePrd = 1-prd.pStopICorr(iTrialCatStop);
                
                % Plot them
                iColor = cndClr(iCohInd, :);
                plot(ssd,respRateObs,'Color',iColor,'Marker',stopIErrorCCorrMrkObs,'LineStyle',stopIErrorCCorrLnObs);
                plot(ssd,respRatePrd,'Color',iColor,'Marker',stopIErrorCCorrMrkPrd,'LineStyle',stopIErrorCCorrLnPrd);
                
                
                % Set axes
                %             switch subject
                %                 case 1
                set(gca,'XLim',[0 750], ...
                    'XTick',0:150:750, ...
                    'YLim',[0 1], ...
                    'YTick',0:0.2:1)
                %                 case 2
                %             set(gca,'XLim',[200 500], ...
                %                 'XTick',100:100:600, ...
                %                 'YLim',[0 1], ...
                %                 'YTick',0:0.2:1)
                %                 otherwise
                %                     error('Need to add axes limits for subject')
                %             end
%             end
%         end
        
    end

end
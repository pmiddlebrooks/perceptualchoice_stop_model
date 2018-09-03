
%%
% Which subject and set of sessions?
subjectName = 'joule';
neuralSet = true;
subSample = false;
if subSample
    addSubSample = 'subSample';
else
    addSubSample = [];
end
if neuralSet
    addData = 'neurophys';
else
     addData = 'behavior';
end

% Set subject-specific variables
switch subjectName
    case 'xena'
        subject = 2;
        modelDate = '2017-05-17';
        model = 478;
        architecture            = 'cffi_ili';
        sessionSet = 'behavior1';
    case 'broca'
        subject = 1;
        switch neuralSet
            case false
                modelDate = '2017-04-13';
                model                   = 79;
                architecture            = 'crace_ili';
                sessionSet = 'behavior2';
            case true
                modelDate = '2017-05-18';
                model                   = 478;
                architecture            = 'cli_ili';
                sessionSet = {...
                    'bp242n02';...
                    'bp244n02';...
                    'bp245n02';...
                    'bp246n02';...
                    'bp247n02'};
        end
    case 'joule'
        subject = 3;
        modelDate = '2017-05-08';
        model                   = 478;
        architecture            = 'cli_ili';
        sessionSet = {...
            'jp110n02';...
            'jp114n04';...
            'jp121n02';...
            'jp124n04';...
            'jp125n04'};
end



% Set Paths
accreRoot = '/gpfs22';
accreHome = '/home/middlepg';
accreScratch = '/scratch/middlepg';
if isdir(fullfile(accreScratch))
    matRoot = fullfile(accreRoot,accreHome,'m-files');
    modelRoot = fullfile(accreScratch,'perceptualchoice_stop_model');
    environ = 'accre';
else
    matRoot = '~/schalllab';
    modelRoot = '~/perceptualchoice_stop_model';
    environ = 'local';
end

addpath(genpath(fullfile(matRoot,'sam')));
addpath(genpath(fullfile(matRoot,'matlab_code_bbz')));
addpath(genpath(fullfile(matRoot,'matlab_file_exchange_tools')));
addpath(genpath(fullfile(matRoot,'cmtb')));
addpath(genpath(fullfile(modelRoot,'src/code',modelDate)));

cd(fullfile(modelRoot,'src/code',modelDate,'matlab'))



% Set simulation variables
trialVar                = true;
dt                      = 1;
optimScope                = 'all';
fileStr.root            = ['~/perceptualchoice_stop_model/data/',modelDate,'/preproc01/subj%.2d/dt%d/%s/%s/'];
fileStr.src             = ['~/perceptualchoice_stop_model/src/code/',modelDate,'/matlab/'];
fileStr.result          = ['~/perceptualchoice_stop_model/results/',modelDate,'/subj%.2d/dt%d/%s/%s/'];
savePlot                = true;




% Get observed inhibition data
opt = ccm_options;
opt.plotFlag = false;
opt.printPlot = false;
data = ccm_inhibition_population(subjectName, sessionSet, opt);


nSession = size(data.ssrtIntWeight, 1);







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. PROCESS INPUTS & SPECIFY VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

% 1.1. Process inputs
% =========================================================================
% nSim = 5000;
nSim = 100;

if trialVar
    trialVarStr = 'trialvar';
else
    trialVarStr = 'notrialvar';
end

% 1.2. Specify dynamic variables
% =========================================================================



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
                standard_figure(1,1,'landscape');

                %                 set_figure({1024,574,'pixels'},{'USLetter','landscape'},{'Helvetica',18});
                p = panel();
                nPlotCol = max(3, nModel);
                nPlotCol = 5;
                p.pack({.1 .3 .3 .3}, num2cell(repmat(1/nPlotCol,1,nPlotCol)));

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

% Specificy observations
obs = SAM.optim.obs;

% For SSRT-only analyses (i.e. this script), get rid of
    % conditions without and Stop RTs (Go trials or nStop = 0
    % trials);
    nStopRT = obs.nStopICorr;
    obs = obs(~isnan(nStopRT) & nStopRT > 0, :);

    
    
% Alter the simulations for efficiency
%             nSim = max(obs.nStopICorr);
firstStopCat = find(strncmp(obs.trialCat, 's', 1), 1);
if subSample
    nSim = max(obs.nTotal) * 2;
else
    nSim = round(mean(obs.nTotal(firstStopCat:end)));
    nSim = round(mean(data.nStopMean));
    nSim = round(mean(data.nStopMeanAll));
end
SAM.sim.n = nSim;

% return



% Loop through number of sessions to simulate
nCondition = length(data.pSignalArray)/2;
    conditionArray = 1:nCondition;
rtStopICorrPrdL = cell(nSession, length(conditionArray));
rtStopICorrPrdR = cell(nSession, length(conditionArray));
for iSim = 1 : nSession
    % Get model predictions and costs
    [cost,altCost,prd] = sam_cost(X,SAM);
    
    % For SSRT-only analyses (i.e. this script), get rid of
    % conditions without and Stop RTs (Go trials or nStop = 0
    % trials);
    prd = prd(~isnan(nStopRT) & nStopRT > 0, :);
    
    % Subsample the number of simulated trials to match the observed
    % numbers
    if subSample
    prd.rtStopICorr = cellfun(@(x) x(randperm(length(x))), prd.rtStopICorr, 'uni', false);
    prd.rtStopICorr = cellfun(@(x,y) x(1:y), prd.rtStopICorr, num2cell(obs.nStopICorr), 'uni', false);
    end
    
    
    
    %             Get SSRTs for each condition
    for iCondition = 1:length(conditionArray)
        % Get the categories
        iTrialCatStopL = cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('stopTrial.*ssd.*c%d.*GO.*r%d',iCondition,1), 'once')),prd.trialCat,'Uni',0));
        iTrialCatStopR = cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('stopTrial.*ssd.*c%d.*GO.*r%d',iCondition,2), 'once')),prd.trialCat,'Uni',0));
                
        % Get Stop RTs
        rtStopICorrPrdL{iSim, iCondition}         = cell2mat(cellfun(@(x) reshape(x, [], 1), prd.rtStopICorr(iTrialCatStopL), 'uni', false));
        rtStopICorrPrdR{iSim, iCondition}         = cell2mat(cellfun(@(x) reshape(x, [], 1), prd.rtStopICorr(iTrialCatStopR), 'uni', false));        
%         rtStopICorrPrdL{iSim, iCondition}         = cell2mat(cellfun(@mean, prd.rtStopICorr(iTrialCatStopL), 'uni', false));
%         rtStopICorrPrdR{iSim, iCondition}         = cell2mat(cellfun(@mean, prd.rtStopICorr(iTrialCatStopR), 'uni', false));        
    end
    
    
    
    
    
end


% Plot observations and predictions
% =========================================================================
%                     sam_plot(SAM,prd);

% Plot it
%                         plotit(SAM,prd,p,model,iModel,modelStr,cost,altCost);

stopICorrMrkPrd         = 'none';
stopICorrLnPrd          = '--';
stopICorrLnWidth        = 2;

stopICorrMrkObs         = 'o';
stopICorrLnObs          = 'none';




%     % SSRTS
%     % ================================================
    p(4,iModel).select();
    p(4,iModel).hold('on');

    % Predicted
%     plot(.1 + (1:nCondition), mean(cellfun(@mean, rtStopICorrPrdL), 1),'Color','r','Marker',stopICorrMrkPrd,'LineStyle',stopICorrLnPrd,'LineWidth',stopICorrLnWidth);
%     plot(.1 + (nCondition*2:-1:nCondition+1), mean(cellfun(@mean, rtStopICorrPrdR), 1),'Color','r','Marker',stopICorrMrkPrd,'LineStyle',stopICorrLnPrd,'LineWidth',stopICorrLnWidth);
    plot(.1 + (1:nCondition), mean(cellfun(@mean, rtStopICorrPrdL), 1),'Color','r','Marker','.','MarkerSize', 10,'LineStyle',stopICorrLnObs,'LineWidth',stopICorrLnWidth);
    plot(.1 + (nCondition*2:-1:nCondition+1), mean(cellfun(@mean, rtStopICorrPrdR), 1),'Color','r','Marker','.','MarkerSize', 10,'LineStyle',stopICorrLnObs,'LineWidth',stopICorrLnWidth);
    
    ciPrdL = 1.96 * std(cellfun(@mean, rtStopICorrPrdL), 1) / sqrt(nSession);
    ciPrdR = 1.96 * std(cellfun(@mean, rtStopICorrPrdR), 1) / sqrt(nSession);
    ciLinesY = [[mean(cellfun(@mean, rtStopICorrPrdL), 1) - ciPrdL, fliplr(mean(cellfun(@mean, rtStopICorrPrdR), 1) - ciPrdR)];...
        [mean(cellfun(@mean, rtStopICorrPrdL), 1) + ciPrdL , fliplr(mean(cellfun(@mean, rtStopICorrPrdR), 1) + ciPrdR)]];
    ciLinesX = .1 + repmat(1:nCondition*2, 2, 1);
    plot(ciLinesX, ciLinesY, 'Color','r','LineStyle','-','LineWidth',stopICorrLnWidth);
    
%     plot(1:nCondition,mean(cellfun(@mean, rtStopICorrPrdL), 1) + ciPrdL,'Color','r','Marker',stopICorrMrkPrd,'LineStyle',stopICorrLnPrd,'LineWidth',stopICorrLnWidth);
%     plot(1:nCondition,mean(cellfun(@mean, rtStopICorrPrdL), 1) - ciPrdL,'Color','r','Marker',stopICorrMrkPrd,'LineStyle',stopICorrLnPrd,'LineWidth',stopICorrLnWidth);
%     plot(nCondition*2:-1:nCondition+1,mean(cellfun(@mean, rtStopICorrPrdR), 1) + ciPrdR,'Color','r','Marker',stopICorrMrkPrd,'LineStyle',stopICorrLnPrd,'LineWidth',stopICorrLnWidth);
%     plot(nCondition*2:-1:nCondition+1,mean(cellfun(@mean, rtStopICorrPrdR), 1) - ciPrdR,'Color','r','Marker',stopICorrMrkPrd,'LineStyle',stopICorrLnPrd,'LineWidth',stopICorrLnWidth);
    

    % Observed
    plot(1:nCondition*2,mean(data.ssrtIntWeight, 1),'Color','k','Marker','.','MarkerSize', 10,'LineStyle',stopICorrLnObs,'LineWidth',stopICorrLnWidth);

    ciObs = 1.96 * std(data.ssrtIntWeight, 1) / sqrt(nSession);
    ciLinesY = [mean(data.ssrtIntWeight, 1) - ciObs; mean(data.ssrtIntWeight, 1) + ciObs];
    ciLinesX = repmat(1:nCondition*2, 2, 1);
    plot(ciLinesX, ciLinesY, 'Color','k','LineStyle','-','LineWidth',stopICorrLnWidth);
%     plot(1:nCondition*2,mean(data.ssrtIntWeight, 1) + ciObs,'Color','k','Marker',stopICorrMrkObs,'LineStyle',stopICorrLnObs,'LineWidth',stopICorrLnWidth);
%     plot(1:nCondition*2,mean(data.ssrtIntWeight, 1) - ciObs,'Color','k','Marker',stopICorrMrkObs,'LineStyle',stopICorrLnObs,'LineWidth',stopICorrLnWidth);
   % Set axes
    switch subject
        case 1
            switch neuralSet
                case true
            set(gca,'XLim',[.5 4.5], ...
                'XTick',1:4, ...
                'YLim',[0 100], ...
                'YTick',0:20:200)
                case false
            set(gca,'XLim',[.5 6.5], ...
                'XTick',1:6, ...
                'YLim',[0 100], ...
                'YTick',0:20:200)
            end
       case 2
            set(gca,'XLim',[.5 6.5], ...
                'XTick',1:6, ...
                'YLim',[0 120], ...
                'YTick',0:20:200)
        case 3
            set(gca,'XLim',[.5 4.5], ...
                'XTick',1:4, ...
                'YLim',[0 120], ...
                'YTick',0:20:200)
        otherwise
            error('Need to add axes limits for subject')
    end






if savePlot
    saveDir             = fullfile(sprintf(fileStr.result,subject,dt,trialVarStr,architecture));
    if exist(saveDir,'dir') ~= 7
        mkdir(saveDir)
    end
    fileName = sprintf('%s_%s_SSRT_Model_Observed_ConfidenceIntervals_nSim_%d_%s', subjectName, addData, nSim, addSubSample);
    print(gcf, fullfile(saveDir, fileName),'-dpdf', '-r300')
%     print(gcf, fullfile(saveDir, fileName),'-dpng')
end








% function plotit(SAM,prd,p,model,iModel,modelStr,cost,altCost)
%
% % Specify colors and line properties
% mapShades           	= [.1 .25 .4 .6 .75 .9];
% cndClr                  = ccm_colormap(mapShades);
% cndClr(1+length(mapShades)/2 : end, :) = flipud(cndClr(1+length(mapShades)/2 : end, :));
%
%
% conditionArray = 1:2;
%
% % goCCorrMrkObs           = 'o';
% % goCCorrLnObs            = 'none';
% % goCCorrMrkPrd           = 'none';
% % goCCorrLnPrd            = '-';
% % goCCorrLnWidth          = 2;
% %
% % goCErrorMrkObs          = '^';
% % goCErrorLnObs           = 'none';
% % goCErrorMrkPrd          = 'none';
% % goCErrorLnPrd           = '-.';
% % goCErrorLnWidth         = 1;
% %
% % stopIErrorCCorrMrkObs   = 's';
% % stopIErrorCCorrLnObs    = 'none';
% % stopIErrorCCorrMrkPrd   = 'none';
% % stopIErrorCCorrLnPrd    = '--';
% % stopIErrorCCorrLnWidth  = 2;
%
% stopICorrMrkPrd         = 'none';
% stopICorrLnPrd          = '--';
% stopICorrLnWidth        = 2;
%
% stopICorrMrkObs         = 'o';
% stopICorrLnObs          = 'none';
% stopICorrLnWidth        = 2;
%
%
%
% % if strcmp(optimScope, 'go') || strcmp(optimScope, 'all')
% %     % RT avgs Go trials
% %     % ================================================
% %     p(2,iModel).select();
% %     p(2,iModel).hold('on');
% %     p(2,iModel).title({sprintf('Model %d',model), ...
% %         sprintf('\\chi^2 = %.1f',cost), ...
% %         sprintf('BIC = %.1f',altCost)});
% %
% %
% %
% %     trialCatGoL = cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('goTrial.*GO.*r%d.*',1), 'once')),prd.trialCat,'Uni',0));
% %     trialCatGoR = cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('goTrial.*GO.*r%d.*',2), 'once')),prd.trialCat,'Uni',0));
% %
% %     rtGoCCorrObsL     = obs.rtGoCCorr(trialCatGoL);
% %     rtGoCCorrPrdL     = prd.rtGoCCorr(trialCatGoL);
% %
% %     plot(1:2,cellfun(@nanmean, rtGoCCorrObsL),'Color','k','Marker',goCCorrMrkObs,'LineStyle',goCCorrLnObs,'LineWidth',goCCorrLnWidth);
% %     plot(1:2,cellfun(@mean, rtGoCCorrPrdL),'Color','k','Marker',goCCorrMrkPrd,'LineStyle',goCCorrLnPrd,'LineWidth',goCCorrLnWidth);
% %
% %     rtGoCCorrObsR     = obs.rtGoCCorr(trialCatGoR);
% %     rtGoCCorrPrdR     = prd.rtGoCCorr(trialCatGoR);
% %
% %     plot(4:-1:3,cellfun(@nanmean, rtGoCCorrObsR),'Color','k','Marker',goCCorrMrkObs,'LineStyle',goCCorrLnObs,'LineWidth',goCCorrLnWidth);
% %     plot(4:-1:3,cellfun(@mean, rtGoCCorrPrdR),'Color','k','Marker',goCCorrMrkPrd,'LineStyle',goCCorrLnPrd,'LineWidth',goCCorrLnWidth);
% %
% %     % Set axes
% %     switch subject
% %         case 1
% %             set(gca,'XLim',[.5 4.5], ...
% %                 'XTick',1:4, ...
% %                 'YLim',[250 650], ...
% %                 'YTick',250:50:650)
% %         case 3
% %             set(gca,'XLim',[.5 4.5], ...
% %                 'XTick',1:4, ...
% %                 'YLim',[250 450], ...
% %                 'YTick',250:50:450)
% %         otherwise
% %             error('Need to add axes limits for subject')
% %     end
% %
% %
% %
% %     % Psychometric function
% %     % ================================================
% %     p(3,iModel).select();
% %     p(3,iModel).hold('on');
% %     %     p(2,iModel).title({sprintf('Model %d',model), ...
% %     %                        modelStr{iModel}, ...
% %     %                        sprintf('\\chi^2 = %.1f',cost), ...
% %     %                        sprintf('BIC = %.1f',altCost)});
% %
% %
% %
% %     pGoCCorrObsL     = 1 - obs.pGoCCorr(trialCatGoL);
% %     pGoCCorrPrdL     = 1 - prd.pGoCCorr(trialCatGoL);
% %
% %     plot(1:2,pGoCCorrObsL,'Color','k','Marker',goCCorrMrkObs,'LineStyle',goCCorrLnObs,'LineWidth',goCCorrLnWidth);
% %     %                         plot(1:3,pGoCCorrPrdL,'Color','k','Marker',goCCorrMrkPrd,'LineStyle',goCCorrLnPrd,'LineWidth',goCCorrLnWidth);
% %
% %     pGoCCorrObsR     = obs.pGoCCorr(trialCatGoR);
% %     pGoCCorrPrdR     = prd.pGoCCorr(trialCatGoR);
% %
% %     plot(4:-1:3,pGoCCorrObsR,'Color','k','Marker',goCCorrMrkObs,'LineStyle',goCCorrLnObs,'LineWidth',goCCorrLnWidth);
% %     plot(1:4,[pGoCCorrPrdL; flipud(pGoCCorrPrdR)],'Color','k','Marker',goCCorrMrkPrd,'LineStyle',goCCorrLnPrd,'LineWidth',goCCorrLnWidth);
% %
% %
% %     set(gca,'XLim',[.5 4.5], ...
% %         'XTick',1:4, ...
% %         'YLim',[0 1], ...
% %         'YTick',0:0.2:1)
% % end
%
%
%
%
%
% if strcmp(optimScope, 'stop') || strcmp(optimScope, 'all')
%
%
%     % RT avgs Stop trials
%     %         ================================================
%     p(2,iModel).select();
%     p(2,iModel).hold('on');
%
%     % Get rid of NaNs in obs.cumProbStopIErrorCCorr
%     for j = 1 : size(obs, 1)
%         if isnan(obs.cumProbStopIErrorCCorr{j})
%             obs.cumProbStopIErrorCCorr{j} = [];
%         end
%     end
%
%
%
%
%     for iCondition = 1:length(conditionArray)
%
%
%
%         iTrialCatStopL = cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('stopTrial.*ssd.*c%d.*GO.*r%d',iCondition,1), 'once')),prd.trialCat,'Uni',0));
%         iTrialCatStopR = cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('stopTrial.*ssd.*c%d.*GO.*r%d',iCondition,2), 'once')),prd.trialCat,'Uni',0));
%
%
%         % Get Stop RTs
%
% %         rtStopIErrorCCorrObsL{iCondition}   = cell2mat(obs.rtStopIErrorCCorr(iTrialCatStopL));
% %         rtStopIErrorCCorrPrdL{iCondition} 	= cell2mat(cellfun(@(x) reshape(x, [], 1), prd.rtStopIErrorCCorr(iTrialCatStopL), 'uni', false));
%         rtStopICorrPrdL{iCondition}         = cell2mat(cellfun(@(x) reshape(x, [], 1), prd.rtStopICorr(iTrialCatStopL), 'uni', false));
%
% %         rtStopIErrorCCorrObsR{iCondition}   = cell2mat(obs.rtStopIErrorCCorr(iTrialCatStopR));
% %         rtStopIErrorCCorrPrdR{iCondition} 	= cell2mat(cellfun(@(x) reshape(x, [], 1), prd.rtStopIErrorCCorr(iTrialCatStopR), 'uni', false));
%         rtStopICorrPrdR{iCondition}         = cell2mat(cellfun(@(x) reshape(x, [], 1), prd.rtStopICorr(iTrialCatStopR), 'uni', false));
%
%
% %         pStopIErrorCCorrObsL(iCondition)   = 1 - sum(obs.nStopIErrorCCorr(iTrialCatStopL)) / sum(obs.nStopIErrorCCorr(iTrialCatStopL) + obs.nStopIErrorCError(iTrialCatStopL));
% %         pStopIErrorCCorrPrdL(iCondition) 	= 1 - sum(prd.nStopIErrorCCorr(iTrialCatStopL)) / sum(prd.nStopIErrorCCorr(iTrialCatStopL) + prd.nStopIErrorCError(iTrialCatStopL));
% %
% %         pStopIErrorCCorrObsR(iCondition)   = sum(obs.nStopIErrorCCorr(iTrialCatStopR)) / sum(obs.nStopIErrorCCorr(iTrialCatStopR) + obs.nStopIErrorCError(iTrialCatStopR));
% %         pStopIErrorCCorrPrdR(iCondition) 	= sum(prd.nStopIErrorCCorr(iTrialCatStopR)) / sum(prd.nStopIErrorCCorr(iTrialCatStopR) + prd.nStopIErrorCError(iTrialCatStopR));
%
%
%
%
%     end
%
% %     plot(1:2,cellfun(@nanmean, rtStopIErrorCCorrObsL),'Color','r','Marker',stopIErrorCCorrMrkObs,'LineStyle',stopIErrorCCorrLnObs,'LineWidth',stopIErrorCCorrLnWidth);
% %     plot(1:2,cellfun(@mean, rtStopIErrorCCorrPrdL),'Color','r','Marker',stopIErrorCCorrMrkPrd,'LineStyle',stopIErrorCCorrLnPrd,'LineWidth',stopIErrorCCorrLnWidth);
% %
% %     plot(4:-1:3,cellfun(@nanmean, rtStopIErrorCCorrObsR),'Color','r','Marker',stopIErrorCCorrMrkObs,'LineStyle',stopIErrorCCorrLnObs,'LineWidth',stopIErrorCCorrLnWidth);
% %     plot(4:-1:3,cellfun(@mean, rtStopIErrorCCorrPrdR),'Color','r','Marker',stopIErrorCCorrMrkPrd,'LineStyle',stopIErrorCCorrLnPrd,'LineWidth',stopIErrorCCorrLnWidth);
% %
% %
% %
% %
% %     % Psychometric function
% %     % ================================================
% %     p(3,iModel).select();
% %     p(3,iModel).hold('on');
% %     %     p(2,iModel).title({sprintf('Model %d',model), ...
% %     %                        modelStr{iModel}, ...
% %     %                        sprintf('\\chi^2 = %.1f',cost), ...
% %     %                        sprintf('BIC = %.1f',altCost)});
% %
% %
% %
% %     pGoCCorrObsL     = 1 - obs.pGoCCorr(trialCatGoL);
% %     pGoCCorrPrdL     = 1 - prd.pGoCCorr(trialCatGoL);
% %
% %     %                         plot(1:3,pGoCCorrPrdL,'Color','k','Marker',goCCorrMrkPrd,'LineStyle',goCCorrLnPrd,'LineWidth',goCCorrLnWidth);
% %
% %     pGoCCorrObsR     = obs.pGoCCorr(trialCatGoR);
% %     pGoCCorrPrdR     = prd.pGoCCorr(trialCatGoR);
% %
% %     plot(1:2,pGoCCorrObsL,'Color','k','Marker',goCCorrMrkObs,'LineStyle',goCCorrLnObs,'LineWidth',goCCorrLnWidth);
% %     plot(4:-1:3,pGoCCorrObsR,'Color','k','Marker',goCCorrMrkObs,'LineStyle',goCCorrLnObs,'LineWidth',goCCorrLnWidth);
% %     plot(1:4,[pGoCCorrPrdL; flipud(pGoCCorrPrdR)],'Color','k','Marker',goCCorrMrkPrd,'LineStyle',goCCorrLnPrd,'LineWidth',goCCorrLnWidth);
% %
% %     plot(1:2,pStopIErrorCCorrObsL,'Color','r','Marker',stopIErrorCCorrMrkObs,'LineStyle',stopIErrorCCorrLnObs,'LineWidth',stopIErrorCCorrLnWidth);
% %     plot(4:-1:3,pStopIErrorCCorrObsR,'Color','r','Marker',stopIErrorCCorrMrkObs,'LineStyle',stopIErrorCCorrLnObs,'LineWidth',stopIErrorCCorrLnWidth);
% %     plot(1:4,[pStopIErrorCCorrPrdL, fliplr(pStopIErrorCCorrPrdR)],'Color','r','Marker',stopIErrorCCorrMrkPrd,'LineStyle',stopIErrorCCorrLnPrd,'LineWidth',stopIErrorCCorrLnWidth);
%
%
%
%
%
%
%     % SSRTS
%     % ================================================
%     p(4,iModel).select();
%     p(4,iModel).hold('on');
%
%     plot(1:2,cellfun(@mean, rtStopICorrPrdL),'Color','r','Marker',stopICorrMrkPrd,'LineStyle',stopICorrLnPrd,'LineWidth',stopICorrLnWidth);
%     plot(4:-1:3,cellfun(@mean, rtStopICorrPrdR),'Color','r','Marker',stopICorrMrkPrd,'LineStyle',stopICorrLnPrd,'LineWidth',stopICorrLnWidth);
%     plot(1:4,mean(data.ssrtIntWeight, 1),'Color','k','Marker',stopICorrMrkObs,'LineStyle',stopICorrLnObs,'LineWidth',stopICorrLnWidth);
%
%     % Set axes
%     switch subject
%         case 1
%             set(gca,'XLim',[.5 4.5], ...
%                 'XTick',1:4, ...
%                 'YLim',[0 100], ...
%                 'YTick',0:20:200)
%         case 3
%             set(gca,'XLim',[.5 4.5], ...
%                 'XTick',1:4, ...
%                 'YLim',[0 100], ...
%                 'YTick',0:20:200)
%         otherwise
%             error('Need to add axes limits for subject')
%     end
%
%
% end
% end



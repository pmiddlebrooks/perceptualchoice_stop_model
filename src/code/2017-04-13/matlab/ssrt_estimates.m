function ssrt_estimates(subject,model,architecture,dt,trialVar,optimScope,fileStr,responseSide, accuracy, defective, savePlot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. PROCESS INPUTS & SPECIFY VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
nSim = 3000;
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
for iSubject = 1:nSubject
    for iArchitecture = 1:nArchitecture
        
%         % Set up the figure and panels
%         %         [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_figure;
%         standard_figure(1,1,'landscape');
%         
%         %                 set_figure({1024,574,'pixels'},{'USLetter','landscape'},{'Helvetica',18});
%         p = panel();
%         p.pack({.1 .25 .25 .25 .15}, num2cell(repmat(1/nModel,1,nModel)));
%         
%         annotation('textbox', [0 0.9 1 0.1], ...
%             'String', sprintf('Subject %d, architecture %s',subject(iSubject),architecture{iArchitecture}), ...
%             'Interpreter','none', ...
%             'EdgeColor', 'none', ...
%             'HorizontalAlignment', 'center');
        
        for iModel = 1:nModel
            
            % Display progress
            fprintf(1,'Working on subject %d, architecture %s, and model %d. \n',subject(iSubject),architecture{iArchitecture},model(iModel));
            
            % Load model-specific fits with cost function values

            ds = dataset('File',fullfile(sprintf(rootDir,subject(iSubject),dt,trialVarStr,architecture{iArchitecture}),sprintf(nameFVal,optimScope,model(iModel))));

            % Load SAM
            load(fullfile(sprintf(rootDir,subject(iSubject),dt,trialVarStr,architecture{iArchitecture}),ds.FileNameSAM{1}),'SAM');
            
            % Extract optimized X
            iBestX = cell2mat(cellfun(@(in1) ~isempty(regexp(in1,'^BestX.*', 'once')),ds.Properties.VarNames,'Uni',0));
            X = double(ds(1,iBestX));
            
            % Alter the simulations for efficiency
            SAM.sim.n = nSim;
            
            % Specificy observations
            obs = SAM.optim.obs;
            
%             % SSD index
%             nSSD = unique(obs.ssd(~isnan(obs.ssd)));
%             switch subject(iSubject)
%                 case 1
%                     iSsd = ceil(length(nSSD) * .65);
%                     iSsd = 4;
%                 case 2
%                     iSsd = 7;
%             end
            
            % Get model predictions and costs
            [cost,altCost,prd] = sam_cost(X,SAM);

            
            
            % Identify trial categories 
            trialCatGoL = cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('goTrial.*GO.*r%d.*',1), 'once')),prd.trialCat,'Uni',0));
            trialCatGoR = cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('goTrial.*GO.*r%d.*',2), 'once')),prd.trialCat,'Uni',0));

            
            
            
            
            
            
            % Plot observations and predictions
            % =========================================================================
            %                     sam_plot(SAM,prd);
            
%             % Plot it
%             plotit(SAM,prd,p,model,iModel,modelStr,cost,altCost, responseSide, accuracy, defective);
%             
%             if savePlot
%                 saveDir             = fullfile(sprintf(fileStr.result,subject(iSubject),dt,trialVarStr,architecture{iArchitecture}));
%                 if exist(saveDir,'dir') ~= 7
%                     mkdir(saveDir)
%                 end
%                 fileName = sprintf('Respond_%s_Accuracy_%s', responseSide, accuracy);
%                 print(gcf, fullfile(saveDir, fileName),'-dpdf', '-r300')
%                 print(gcf, fullfile(saveDir, fileName),'-dpng')
%             end
%             clear SAM prd

        end
%         clf
    end
end








%     function plotit(SAM,prd,p,model,iModel,modelStr,cost,altCost, responseSide, accuracy, defective)
%         
%         % Specify colors and line properties
%         mapShades           	= [.1 .25 .4 .6 .75 .9];
%         cndClr                  = ccm_colormap(mapShades);
%         cndClr(1+length(mapShades)/2 : end, :) = flipud(cndClr(1+length(mapShades)/2 : end, :));
%         
%         
%         switch responseSide
%             case 'both'
%                 responseArray = [1 2];
%             case 'left'
%                 responseArray = 1;
%             case 'right'
%                 responseArray = 2;
%             otherwise
%                 error('responseSide need to be left, right, or both');
%         end
%         conditionArray = 1:3;
%         
%         goCCorrMrkObs           = 'o';
%         goCCorrLnObs            = 'none';
%         goCCorrMrkPrd           = 'none';
%         goCCorrLnPrd            = '-';
%         goCCorrLnWidth          = 2;
%         
%         goCErrorMrkObs          = '^';
%         goCErrorLnObs           = 'none';
%         goCErrorMrkPrd          = 'none';
%         goCErrorLnPrd           = '-.';
%         goCErrorLnWidth         = 1;
%         
%         stopIErrorCCorrMrkObs   = 's';
%         stopIErrorCCorrLnObs    = 'none';
%         stopIErrorCCorrMrkPrd   = 'none';
%         stopIErrorCCorrLnPrd    = '--';
%         stopIErrorCCorrLnWidth  = 2;
%         
%         stopICorrMrkPrd         = 'none';
%         stopICorrLnPrd          = ':';
%         stopICorrLnWidth        = 1;
%         
%         
%         
%         
%         if strcmp(optimScope, 'go') || strcmp(optimScope, 'all')
%             % RT distribution Go trials
%             % ================================================
%             p(2,iModel).select();
%             p(2,iModel).hold('on');
% %             p(2,iModel).title({sprintf('Model %d',model(iModel)), ...
% %                 sprintf('\\chi^2 = %.1f',cost), ...
% %                 sprintf('BIC = %.1f',altCost)});
%             p(2,iModel).title({sprintf('%d',round(cost))});
%             
%             for iResponse = 1:length(responseArray)
%                 for iCondition = 1:length(conditionArray);
%                     
%                     iRsp = responseArray(iResponse);
%                     iCnd = conditionArray(iCondition);
%                     
%                     % Identify the relevant rows in the dataset array
%                     %         iTrialCatGo = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('goTrial.*GO:c%d',iCnd), 'once')),prd.trialCat,'Uni',0)));
%                     iTrialCatGo = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('goTrial.*c%d.*GO.*r%d',iCnd,iRsp), 'once')),prd.trialCat,'Uni',0)));
%                     if isempty(iTrialCatGo)
%                         iTrialCatGo = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('goTrial_c%d.*',iCnd), 'once')),prd.trialCat,'Uni',0)));
%                     end
%                     
%                     % GoCCorr trials
%                     % -----------------------------------------------------------------
%                     if strcmp(accuracy, 'both') || strcmp(accuracy, 'correct')
%                         rtGoCCorrObs     = obs.rtQGoCCorr{iTrialCatGo};
%                         cumPGoCCorrObs   = obs.cumProbGoCCorr{iTrialCatGo};
%                         
%                         rtGoCCorrPrd     = prd.rtGoCCorr{iTrialCatGo};
%                         cumPGoCCorrPrd   = cmtb_edf(prd.rtGoCCorr{iTrialCatGo}(:),prd.rtGoCCorr{iTrialCatGo}(:));
%                         
%                         if defective
%                             cumPGoCCorrPrd   = cumPGoCCorrPrd * prd.pGoCCorr(iTrialCatGo);
%                             cumPGoCCorrObs   = obs.cumProbDefectiveGoCCorr{iTrialCatGo};
%                         end
%                         
%                         % Plot it
%                         iColor = cndClr((length(conditionArray) * (iRsp -1) + iCondition), :);
%                         plot(rtGoCCorrObs,cumPGoCCorrObs,'Color',iColor,'Marker',goCCorrMrkObs,'LineStyle',goCCorrLnObs,'LineWidth',goCCorrLnWidth);
%                         plot(rtGoCCorrPrd,cumPGoCCorrPrd,'Color',iColor,'Marker',goCCorrMrkPrd,'LineStyle',goCCorrLnPrd,'LineWidth',goCCorrLnWidth);
%                         
%                     end
%                     
%                     % GoCError trials
%                     % -----------------------------------------------------------------
%                     if ~isempty(obs.rtQGoCError{iTrialCatGo}) && obs.nGoCError(iTrialCatGo) > 10 && (strcmp(accuracy, 'both') || strcmp(accuracy, 'error'))
%                         
%                         rtGoCErrorObs     = obs.rtQGoCError{iTrialCatGo};
%                         cumPGoCErrorObs   = obs.cumProbGoCError{iTrialCatGo};
%                         
%                         rtGoCErrorPrd     = prd.rtGoCError{iTrialCatGo};
%                         cumPGoCErrorPrd   = cmtb_edf(prd.rtGoCError{iTrialCatGo}(:),prd.rtGoCError{iTrialCatGo}(:));
%                         
%                         if defective
%                             cumPGoCErrorPrd   = cumPGoCErrorPrd * prd.pGoCError(iTrialCatGo);
%                             cumPGoCErrorObs   = obs.cumProbDefectiveGoCError{iTrialCatGo};
%                         end
%                         
%                         % Plot it
%                         iColor = cndClr((length(conditionArray) * (iRsp -1) + iCondition), :);
%                         plot(rtGoCErrorObs,cumPGoCErrorObs,'Color',iColor,'Marker',goCErrorMrkObs,'LineStyle',goCErrorLnObs,'LineWidth',goCErrorLnWidth);
%                         plot(rtGoCErrorPrd,cumPGoCErrorPrd,'Color',iColor,'Marker',goCErrorMrkPrd,'LineStyle',goCErrorLnPrd,'LineWidth',goCErrorLnWidth);
%                         
%                     end
%                 end
%             end
%             
%             % Set axes
%             switch subject(iSubject)
%                 case 1
%                     set(gca,'XLim',[200 700], ...
%                         'XTick',100:100:700, ...
%                         'YLim',[0 1], ...
%                         'YTick',0:0.2:1)
%                 case 2
%                     set(gca,'XLim',[200 500], ...
%                         'XTick',100:100:600, ...
%                         'YLim',[0 1], ...
%                         'YTick',0:0.2:1)
%                 otherwise
%                     error('Need to add axes limits for subject')
%             end
%             if iModel > 1
%                 set(gca, 'YTick', [])
%             end
%             
%         end
%         
%         if strcmp(optimScope, 'stop') || strcmp(optimScope, 'all')
%             
%             
%             % RT distribution Stop trials
%             %         ================================================
%             p(3,iModel).select();
%             p(3,iModel).hold('on');
%             
%             % Get rid of NaNs in obs.cumProbStopIErrorCCorr
%             for j = 1 : size(obs, 1)
%                 if isnan(obs.cumProbStopIErrorCCorr{j})
%                     obs.cumProbStopIErrorCCorr{j} = [];
%                 end
%             end
%             for iResponse = 1:length(responseArray)
%                 for iCondition = 1:length(conditionArray);
%                     
%                     iRsp = responseArray(iResponse);
%                     iCnd = conditionArray(iCondition);
%                     
%                     % Identify the relevant rows in the dataset array
%                     iTrialCatStop = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('stopTrial.*ssd%d.*GO.*r%d.*c%d',iSsd,iRsp,iCnd), 'once')),prd.trialCat,'Uni',0)));
%                     
%                     if isempty(iTrialCatStop)
%                         iTrialCatStop = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('stopTrial.*ssd%d.*c%d.*r%d.*',iSsd,iCnd,iRsp), 'once')),prd.trialCat,'Uni',0)));
%                     end
%                     
%                     % StopIErrorCCorr trials
%                     % -----------------------------------------------------------------
%                     rtStopIErrorCCorrObs   = obs.rtQStopIErrorCCorr{iTrialCatStop};
%                     cumPStopIErrorCCorrObs = obs.cumProbStopIErrorCCorr{iTrialCatStop};
%                     
%                     rtStopIErrorCCorrPrd     = prd.rtStopIErrorCCorr{iTrialCatStop};
%                     cumPStopIErrorCCorrPrd   = cmtb_edf(prd.rtStopIErrorCCorr{iTrialCatStop}(:),prd.rtStopIErrorCCorr{iTrialCatStop}(:));
%                     
%                     if defective
%                         cumPStopIErrorCCorrPrd   = cumPStopIErrorCCorrPrd * prd.pStopIErrorCCorr(iTrialCatStop);
% %                         if sum(~isnan(cumPStopIErrorCCorrPrd)) == 0
% %                             cumPStopIErrorCCorrPrd = []
% %                         end
%                         cumPStopIErrorCCorrObs   = cumPStopIErrorCCorrObs * obs.pStopIErrorCCorr(iTrialCatStop);
% %                         cumPStopIErrorCCorrObs   = obs.cumProbDefectiveStopIErrorCCorr{iTrialCatGo}
%                     end
%                     
%                     % Plot it
%                     iColor = cndClr((length(conditionArray) * (iRsp -1) + iCondition), :);
%                     plot(rtStopIErrorCCorrObs,cumPStopIErrorCCorrObs,'Color',iColor,'Marker',stopIErrorCCorrMrkObs,'LineStyle',stopIErrorCCorrLnObs,'LineWidth',stopIErrorCCorrLnWidth);
%                     plot(rtStopIErrorCCorrPrd,cumPStopIErrorCCorrPrd,'Color',iColor,'Marker',stopIErrorCCorrMrkPrd,'LineStyle',stopIErrorCCorrLnPrd,'LineWidth',stopIErrorCCorrLnWidth);
%                     
%                     % StopICorr trials
%                     % -----------------------------------------------------------------
%                     rtStopICorrPrd     = prd.rtStopICorr{iTrialCatStop};
%                     cumPStopICorrPrd   = cmtb_edf(prd.rtStopICorr{iTrialCatStop}(:),prd.rtStopICorr{iTrialCatStop}(:));
%                     
%                     % Plot it
%                     iColor = 'r';%cndClr((length(conditionArray) * (iRsp -1) + iCondition), :);
%                     plot(rtStopICorrPrd,cumPStopICorrPrd,'Color',iColor,'Marker',stopICorrMrkPrd,'LineStyle',stopICorrLnPrd,'LineWidth',stopICorrLnWidth);
%                     
%                 end
%             end
%             
%             % Set axes
%             switch subject(iSubject)
%                 case 1
%                     set(gca,'XLim',[200 700], ...
%                         'XTick',100:100:700, ...
%                         'YLim',[0 1], ...
%                         'YTick',0:0.2:1)
%                 case 2
%                     set(gca,'XLim',[200 500], ...
%                         'XTick',100:100:600, ...
%                         'YLim',[0 1], ...
%                         'YTick',0:0.2:1)
%                 otherwise
%                     error('Need to add axes limits for subject')
%             end
%             
%             
%             
%             % Inhibition function
%             % ================================================
%             p(4,iModel).select();
%             p(4,iModel).hold('on');
%             %     p(2,iModel).title({sprintf('Model %d',model(iModel)), ...
%             %                        modelStr{iModel}, ...
%             %                        sprintf('\\chi^2 = %.1f',cost), ...
%             %                        sprintf('BIC = %.1f',altCost)});
%             
%             for iResponse = 1:length(responseArray)
%                 for iCondition = 1:length(conditionArray);
%                     
%                     iRsp = responseArray(iResponse);
%                     iCnd = conditionArray(iCondition);
%                     
%                     % Identify the relevant rows in the dataset array
% %                     iTrialCatStop = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('stopTrial.*ssd%d.*GO.*r%d.*c%d',iSsd,iRsp,iCnd), 'once')),prd.trialCat,'Uni',0)));
% %                     
% %                     if isempty(iTrialCatStop)
% %                         iTrialCatStop = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('stopTrial.*ssd%d.*c%d.*r%d.*',iSsd,iCnd,iRsp), 'once')),prd.trialCat,'Uni',0)));
% %                     end
% %                     iTrialCat = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('stopTrial.*GO:c%d',iCnd), 'once')),prd.trialCat,'Uni',0)));
%                     iTrialCatStop = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('stopTrial.*GO.*r%d.*c%d',iRsp,iCnd), 'once')),prd.trialCat,'Uni',0)));
%                      if isempty(iTrialCatStop)
%                                     iTrialCatStop = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('stopTrial.*c%d.*GO.*r%d',iCnd,iRsp), 'once')),prd.trialCat,'Uni',0)));
%                      end   
%                     % Get the SSDs and observed and predicted response rates
%                     ssd = obs.ssd(iTrialCatStop);
%                     respRateObs = 1-obs.pStopICorr(iTrialCatStop);
%                     respRatePrd = 1-prd.pStopICorr(iTrialCatStop);
%                     
%                     % Plot them
%                     iColor = cndClr((length(conditionArray) * (iRsp -1) + iCondition), :);
%                     plot(ssd,respRateObs,'Color',iColor,'Marker',stopIErrorCCorrMrkObs,'LineStyle',stopIErrorCCorrLnObs);
%                     plot(ssd,respRatePrd,'Color',iColor,'Marker',stopIErrorCCorrMrkPrd,'LineStyle',stopIErrorCCorrLnPrd);
%                     
%                 end
%             end
%             
%             % Set axes
%             %             switch subject(iSubject)
%             %                 case 1
%             set(gca,'XLim',[0 750], ...
%                 'XTick',0:150:750, ...
%                 'YLim',[0 1], ...
%                 'YTick',0:0.2:1)
%             %                 case 2
%             %             set(gca,'XLim',[200 500], ...
%             %                 'XTick',100:100:600, ...
%             %                 'YLim',[0 1], ...
%             %                 'YTick',0:0.2:1)
%             %                 otherwise
%             %                     error('Need to add axes limits for subject')
%             %             end
%             
%         end
%     end

end
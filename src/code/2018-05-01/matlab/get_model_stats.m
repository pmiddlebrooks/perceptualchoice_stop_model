function [goTable, goErrorTable, stopTable, stopErrorTable, ssrtTable] = get_model_stats(subject,model,architecture,dt,trialVar,optimScope,fileStr, savePlot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. PROCESS INPUTS & SPECIFY VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

% 1.1. Process inputs
% =========================================================================
nSim = 5000;
% nSim = 100;

% 1.0 hard-code which conditions and responses to plot for now
if subject == 2 || (subject == 1 && model == 79)
conditionArray       = 1:3; % 1 is easy, 2 is hard
else
conditionArray       = 1:2; % 1 is easy, 2 is hard
end

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



% Loop over subjects, architectures, and models
warning('off','stats:dataset:ModifiedVarnames')
for iSubject = 1:nSubject
    
    
        
         
        for iModel = 1:nModel
            
            % Display progress
            fprintf(1,'Working on subject %d, architecture %s, and model %d. \n',subject(iSubject),architecture,model(iModel));
            
            % Load model-specific fits with cost function values
            ds = dataset('File',fullfile(sprintf(rootDir,subject(iSubject),dt,trialVarStr,architecture),sprintf(nameFVal,optimScope,model(iModel))));
            
            % Load SAM
            load(fullfile(sprintf(rootDir,subject(iSubject),dt,trialVarStr,architecture),ds.FileNameSAM{1}),'SAM');
            
            
            % Extract optimized X
            iBestX = cell2mat(cellfun(@(in1) ~isempty(regexp(in1,'^BestX.*', 'once')),ds.Properties.VarNames,'Uni',0));
            X = double(ds(1,iBestX));
            
            % Alter the simulations for efficiency
            SAM.sim.n = nSim;
            
            % Specificy observations
            obs = SAM.optim.obs;
            
            
            % Get model predictions and costs
            [cost,altCost,prd] = sam_cost(X,SAM);
            
            % Plot observations and predictions
            % =========================================================================
            %                     sam_plot(SAM,prd);
            
            % Plot it
            getstats(SAM,prd,model,iModel,modelStr,cost,altCost);
            
%             if savePlot
%                 saveDir             = fullfile(sprintf(fileStr.result,subject(iSubject),dt,trialVarStr,architecture));
%                 if exist(saveDir,'dir') ~= 7
%                     mkdir(saveDir)
%                 end
%                 fileName = sprintf('Summary_Stats_Trials_Bins');
%                 print(gcf, fullfile(saveDir, fileName),'-dpdf', '-r300')
%                 print(gcf, fullfile(saveDir, fileName),'-dpng')
%                 dataFileName = sprintf('Summary_Stats_Trials_Bins_%d',model(iModel));
%                 save(fullfile(saveDir, dataFileName), 'prd', 'obs','cost','altCost')
%             end
    end
end








    function getstats(SAM,prd,model,iModel,modelStr,cost,altCost)
        
        
        
        if strcmp(optimScope, 'go') || strcmp(optimScope, 'all')
            % RT avgs Go trials
            % ================================================
            trialCatGoL = cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('goTrial.*GO.*r%d.*',1), 'once')),prd.trialCat,'Uni',0));
            trialCatGoR = cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('goTrial.*GO.*r%d.*',2), 'once')),prd.trialCat,'Uni',0));
            
%             rtGoCCorrObsL     = obs.rtGoCCorr(trialCatGoL);
            rtGoCCorrPrdL     = prd.rtGoCCorr(trialCatGoL)';
            rtGoCErrPrdL     = prd.rtGoCError(trialCatGoL)';
            
%             plot(1:2,cellfun(@nanmean, rtGoCCorrObsL),'Color','k','Marker',goCCorrMrkObs,'LineStyle',goCCorrLnObs,'LineWidth',goCCorrLnWidth);
%             plot(1:2,cellfun(@mean, rtGoCCorrPrdL),'Color','k','Marker',goCCorrMrkPrd,'LineStyle',goCCorrLnPrd,'LineWidth',goCCorrLnWidth);
            
%             rtGoCCorrObsR     = obs.rtGoCCorr(trialCatGoR);
            rtGoCCorrPrdR     = flipud(prd.rtGoCCorr(trialCatGoR))';
            rtGoCErrPrdR     = flipud(prd.rtGoCError(trialCatGoR))';
            
%             plot(4:-1:3,cellfun(@nanmean, rtGoCCorrObsR),'Color','k','Marker',goCCorrMrkObs,'LineStyle',goCCorrLnObs,'LineWidth',goCCorrLnWidth);
%             plot(4:-1:3,cellfun(@mean, rtGoCCorrPrdR),'Color','k','Marker',goCCorrMrkPrd,'LineStyle',goCCorrLnPrd,'LineWidth',goCCorrLnWidth);
            
            rtPrdCohGo = [rtGoCCorrPrdL rtGoCCorrPrdR];
            rtPrdMeanGo = cell2mat(cellfun(@mean, rtPrdCohGo, 'uni',0));
            rtPrdStdGo = cell2mat(cellfun(@std, rtPrdCohGo, 'uni',0));
            rtPrdSemGo = rtPrdStdGo ./ sqrt(cell2mat(cellfun(@length, rtPrdCohGo, 'uni',0)));

                        rtPrdCohGoErr = [rtGoCErrPrdL rtGoCErrPrdR];
            rtPrdMeanGoErr = cell2mat(cellfun(@mean, rtPrdCohGoErr, 'uni',0));
            rtPrdStdGoErr = cell2mat(cellfun(@std, rtPrdCohGoErr, 'uni',0));
            rtPrdSemGoErr = rtPrdStdGoErr ./ sqrt(cell2mat(cellfun(@length, rtPrdCohGoErr, 'uni',0)));

            
            goTable = array2table(round([rtPrdMeanGo; rtPrdStdGo; rtPrdSemGo]));
            goErrorTable = array2table(round([rtPrdMeanGoErr; rtPrdStdGoErr; rtPrdSemGoErr]));
        end
        
        
        
        
        
        if strcmp(optimScope, 'stop') || strcmp(optimScope, 'all')
            
            
            % RT avgs Stop trials
            %         ================================================
            
            % Get rid of NaNs in obs.cumProbStopIErrorCCorr
%             for j = 1 : size(obs, 1)
%                 if isnan(obs.cumProbStopIErrorCCorr{j})
%                     obs.cumProbStopIErrorCCorr{j} = [];
%                 end
%             end
            
            
            
            
            for iCondition = 1:length(conditionArray)
                
                
                
                iTrialCatStopL = cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('stopTrial.*ssd.*c%d.*GO.*r%d',iCondition,1), 'once')),prd.trialCat,'Uni',0));
                iTrialCatStopR = cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('stopTrial.*ssd.*c%d.*GO.*r%d',iCondition,2), 'once')),prd.trialCat,'Uni',0));
                
                
                % Get Stop RTs
                
%                 rtStopIErrorCCorrObsL{iCondition}   = cell2mat(obs.rtStopIErrorCCorr(iTrialCatStopL));
                rtStopIErrorCCorrPrdL{iCondition} 	= cell2mat(cellfun(@(x) reshape(x, [], 1), prd.rtStopIErrorCCorr(iTrialCatStopL), 'uni', false));
                rtStopIErrorCErrPrdL{iCondition} 	= cell2mat(cellfun(@(x) reshape(x, [], 1), prd.rtStopIErrorCError(iTrialCatStopL), 'uni', false));
                rtStopICorrPrdL{iCondition}         = cell2mat(cellfun(@(x) reshape(x, [], 1), prd.rtStopICorr(iTrialCatStopL), 'uni', false));
                
%                 rtStopIErrorCCorrObsR{iCondition}   = cell2mat(obs.rtStopIErrorCCorr(iTrialCatStopR));
                rtStopIErrorCCorrPrdR{iCondition} 	= cell2mat(cellfun(@(x) reshape(x, [], 1), prd.rtStopIErrorCCorr(iTrialCatStopR), 'uni', false));
                rtStopIErrorCErrPrdR{iCondition} 	= cell2mat(cellfun(@(x) reshape(x, [], 1), prd.rtStopIErrorCError(iTrialCatStopR), 'uni', false));
                rtStopICorrPrdR{iCondition}         = cell2mat(cellfun(@(x) reshape(x, [], 1), prd.rtStopICorr(iTrialCatStopR), 'uni', false));
                
                
%                 pStopIErrorCCorrObsL(iCondition)   = 1 - sum(obs.nStopIErrorCCorr(iTrialCatStopL)) / sum(obs.nStopIErrorCCorr(iTrialCatStopL) + obs.nStopIErrorCError(iTrialCatStopL));
%                 pStopIErrorCCorrPrdL(iCondition) 	= 1 - sum(prd.nStopIErrorCCorr(iTrialCatStopL)) / sum(prd.nStopIErrorCCorr(iTrialCatStopL) + prd.nStopIErrorCError(iTrialCatStopL));
%                 
%                 pStopIErrorCCorrObsR(iCondition)   = sum(obs.nStopIErrorCCorr(iTrialCatStopR)) / sum(obs.nStopIErrorCCorr(iTrialCatStopR) + obs.nStopIErrorCError(iTrialCatStopR));
%                 pStopIErrorCCorrPrdR(iCondition) 	= sum(prd.nStopIErrorCCorr(iTrialCatStopR)) / sum(prd.nStopIErrorCCorr(iTrialCatStopR) + prd.nStopIErrorCError(iTrialCatStopR));
%                 
                
                
                
            end
            
%             plot(1:2,cellfun(@nanmean, rtStopIErrorCCorrObsL),'Color','r','Marker',stopIErrorCCorrMrkObs,'LineStyle',stopIErrorCCorrLnObs,'LineWidth',stopIErrorCCorrLnWidth);
%             plot(1:2,cellfun(@mean, rtStopIErrorCCorrPrdL),'Color','r','Marker',stopIErrorCCorrMrkPrd,'LineStyle',stopIErrorCCorrLnPrd,'LineWidth',stopIErrorCCorrLnWidth);
%             
%             plot(4:-1:3,cellfun(@nanmean, rtStopIErrorCCorrObsR),'Color','r','Marker',stopIErrorCCorrMrkObs,'LineStyle',stopIErrorCCorrLnObs,'LineWidth',stopIErrorCCorrLnWidth);
%             plot(4:-1:3,cellfun(@mean, rtStopIErrorCCorrPrdR),'Color','r','Marker',stopIErrorCCorrMrkPrd,'LineStyle',stopIErrorCCorrLnPrd,'LineWidth',stopIErrorCCorrLnWidth);
            
            rtPrdCohStop = [rtStopIErrorCCorrPrdL fliplr(rtStopIErrorCCorrPrdR)];
            rtPrdMeanStop = cell2mat(cellfun(@mean, rtPrdCohStop, 'uni',0));
            rtPrdStdStop = cell2mat(cellfun(@std, rtPrdCohStop, 'uni',0));
            rtPrdSemStop = rtPrdStdStop ./ sqrt(cell2mat(cellfun(@length, rtPrdCohStop, 'uni',0)));
            
            rtPrdCohStopErr = [rtStopIErrorCErrPrdL fliplr(rtStopIErrorCErrPrdR)];
            rtPrdMeanStopErr = cell2mat(cellfun(@mean, rtPrdCohStopErr, 'uni',0));
            rtPrdStdStopErr = cell2mat(cellfun(@std, rtPrdCohStopErr, 'uni',0));
            rtPrdSemStopErr = rtPrdStdStopErr ./ sqrt(cell2mat(cellfun(@length, rtPrdCohStopErr, 'uni',0)));
            
            ssrtPrdCohStop = [rtStopICorrPrdL fliplr(rtStopICorrPrdR)];
            ssrtPrdMeanStop = cell2mat(cellfun(@mean, ssrtPrdCohStop, 'uni',0));
            ssrtPrdStdStop = cell2mat(cellfun(@std, ssrtPrdCohStop, 'uni',0));
            ssrtPrdSemStop = ssrtPrdStdStop ./ sqrt(cell2mat(cellfun(@length, ssrtPrdCohStop, 'uni',0)));
            
            stopTable = array2table(round([rtPrdMeanStop; rtPrdStdStop; rtPrdSemStop]));
            stopErrorTable = array2table(round([rtPrdMeanStopErr; rtPrdStdStopErr; rtPrdSemStopErr]));
            
            
            ssrtTable = array2table(round([ssrtPrdMeanStop; ssrtPrdStdStop; ssrtPrdSemStop]));
            
            
            
            % SSRTS
            % ================================================
            
%             plot(1:2,cellfun(@mean, rtStopICorrPrdL),'Color','r','Marker',stopICorrMrkPrd,'LineStyle',stopICorrLnPrd,'LineWidth',stopICorrLnWidth);
%             plot(4:-1:3,cellfun(@mean, rtStopICorrPrdR),'Color','r','Marker',stopICorrMrkPrd,'LineStyle',stopICorrLnPrd,'LineWidth',stopICorrLnWidth);
%             plot(1:4,mean(data.ssrtIntWeight, 1),'Color','k','Marker',stopICorrMrkObs,'LineStyle',stopICorrLnObs,'LineWidth',stopICorrLnWidth);
            
            
        end
    end

end
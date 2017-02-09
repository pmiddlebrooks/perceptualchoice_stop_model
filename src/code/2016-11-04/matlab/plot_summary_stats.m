function plot_summary_stats(subject,model,architecture,dt,trialVar,optimScope,fileStr, savePlot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. PROCESS INPUTS & SPECIFY VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

% 1.1. Process inputs
% =========================================================================
nSim = 2500;

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
        
        % Set up the figure and panels
        %         [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_figure;
        standard_figure(1,1,'landscape');
        
        %                 set_figure({1024,574,'pixels'},{'USLetter','landscape'},{'Helvetica',18});
        p = panel();
        p.pack({.1 .3 .3 .3}, num2cell(repmat(1/nModel,1,nModel)));
        
        annotation('textbox', [0 0.9 1 0.1], ...
            'String', sprintf('Subject %d, architecture %s',subject(iSubject),architecture{iArchitecture}), ...
            'Interpreter','none', ...
            'EdgeColor', 'none', ...
            'HorizontalAlignment', 'center');
        
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
            
            
            % Get model predictions and costs
            [cost,altCost,prd] = sam_cost(X,SAM);
            
            % Plot observations and predictions
            % =========================================================================
            %                     sam_plot(SAM,prd);
            
            % Plot it
            plotit(SAM,prd,p,model,iModel,modelStr,cost,altCost);
            
            if savePlot
                saveDir             = fullfile(sprintf(fileStr.result,subject(iSubject),dt,trialVarStr,architecture{iArchitecture}));
                if exist(saveDir,'dir') ~= 7
                    mkdir(saveDir)
                end
                fileName = sprintf('Summary_Stats');
                print(gcf, fullfile(saveDir, fileName),'-dpdf', '-r300')
                print(gcf, fullfile(saveDir, fileName),'-dpng')
            end
        end
    end
end








    function plotit(SAM,prd,p,model,iModel,modelStr,cost,altCost)
        
        % Specify colors and line properties
        mapShades           	= [.1 .25 .4 .6 .75 .9];
        cndClr                  = ccm_colormap(mapShades);
        cndClr(1+length(mapShades)/2 : end, :) = flipud(cndClr(1+length(mapShades)/2 : end, :));
        
        
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
        conditionArray = 1:3;
        
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
        
        stopICorrMrkPrd         = 'none';
        stopICorrLnPrd          = ':';
        stopICorrLnWidth        = 1;
        
        
        
        
        if strcmp(optimScope, 'go') || strcmp(optimScope, 'all')
            % RT avgs Go trials
            % ================================================
            p(2,iModel).select();
            p(2,iModel).hold('on');
            p(2,iModel).title({sprintf('Model %d',model(iModel)), ...
                sprintf('\\chi^2 = %.1f',cost), ...
                sprintf('BIC = %.1f',altCost)});
            
            
            
            trialCatGoL = cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('goTrial.*GO.*r%d.*',1), 'once')),prd.trialCat,'Uni',0));
            trialCatGoR = cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('goTrial.*GO.*r%d.*',2), 'once')),prd.trialCat,'Uni',0));
            
            rtGoCCorrObsL     = obs.rtGoCCorr(trialCatGoL);
            rtGoCCorrPrdL     = prd.rtGoCCorr(trialCatGoL);
            
            plot(1:3,cellfun(@nanmean, rtGoCCorrObsL),'Color','k','Marker',goCCorrMrkObs,'LineStyle',goCCorrLnObs,'LineWidth',goCCorrLnWidth);
            plot(1:3,cellfun(@mean, rtGoCCorrPrdL),'Color','k','Marker',goCCorrMrkPrd,'LineStyle',goCCorrLnPrd,'LineWidth',goCCorrLnWidth);
            
            rtGoCCorrObsR     = obs.rtGoCCorr(trialCatGoR);
            rtGoCCorrPrdR     = prd.rtGoCCorr(trialCatGoR);
            
            plot(6:-1:4,cellfun(@nanmean, rtGoCCorrObsR),'Color','k','Marker',goCCorrMrkObs,'LineStyle',goCCorrLnObs,'LineWidth',goCCorrLnWidth);
            plot(6:-1:4,cellfun(@mean, rtGoCCorrPrdR),'Color','k','Marker',goCCorrMrkPrd,'LineStyle',goCCorrLnPrd,'LineWidth',goCCorrLnWidth);
            
            % Set axes
            switch subject(iSubject)
                case 1
            set(gca,'XLim',[.5 6.5], ...
                'XTick',1:6, ...
                'YLim',[250 450], ...
                'YTick',250:50:450)
                case 2
           set(gca,'XLim',[.5 6.5], ...
                'XTick',1:6, ...
                'YLim',[200 400], ...
                'YTick',250:50:450)
                otherwise
                    error('Need to add axes limits for subject')
            end
            
            
            
            % Psychometric function
            % ================================================
            p(3,iModel).select();
            p(3,iModel).hold('on');
            %     p(2,iModel).title({sprintf('Model %d',model(iModel)), ...
            %                        modelStr{iModel}, ...
            %                        sprintf('\\chi^2 = %.1f',cost), ...
            %                        sprintf('BIC = %.1f',altCost)});
            
            
            
            pGoCCorrObsL     = 1 - obs.pGoCCorr(trialCatGoL);
            pGoCCorrPrdL     = 1 - prd.pGoCCorr(trialCatGoL);
            
            plot(1:3,pGoCCorrObsL,'Color','k','Marker',goCCorrMrkObs,'LineStyle',goCCorrLnObs,'LineWidth',goCCorrLnWidth);
            %                         plot(1:3,pGoCCorrPrdL,'Color','k','Marker',goCCorrMrkPrd,'LineStyle',goCCorrLnPrd,'LineWidth',goCCorrLnWidth);
            
            pGoCCorrObsR     = obs.pGoCCorr(trialCatGoR);
            pGoCCorrPrdR     = prd.pGoCCorr(trialCatGoR);
            
            plot(6:-1:4,pGoCCorrObsR,'Color','k','Marker',goCCorrMrkObs,'LineStyle',goCCorrLnObs,'LineWidth',goCCorrLnWidth);
            plot(1:6,[pGoCCorrPrdL; flipud(pGoCCorrPrdR)],'Color','k','Marker',goCCorrMrkPrd,'LineStyle',goCCorrLnPrd,'LineWidth',goCCorrLnWidth);
            
            
            set(gca,'XLim',[.5 6.5], ...
                'XTick',1:6, ...
                'YLim',[0 1], ...
                'YTick',0:0.2:1)
        end
        
        
        
        
        
        if strcmp(optimScope, 'stop') || strcmp(optimScope, 'all')
            
            
            % RT avgs Stop trials
            %         ================================================
            p(2,iModel).select();
            p(2,iModel).hold('on');
            
            % Get rid of NaNs in obs.cumProbStopIErrorCCorr
            for j = 1 : size(obs, 1)
                if isnan(obs.cumProbStopIErrorCCorr{j})
                    obs.cumProbStopIErrorCCorr{j} = [];
                end
            end
            
            
            
            
            for iCondition = 1:length(conditionArray);
                
%                 iTrialCatStopL = cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('stopTrial.*ssd.*GO.*r%d.*c%d',1,iCondition), 'once')),prd.trialCat,'Uni',0));
%                 iTrialCatStopR = cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('stopTrial.*ssd.*GO.*r%d.*c%d',2,iCondition), 'once')),prd.trialCat,'Uni',0));
                
                iTrialCatStopL = cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('stopTrial.*ssd.*c%d.*GO.*r%d',iCondition,1), 'once')),prd.trialCat,'Uni',0));
                iTrialCatStopR = cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('stopTrial.*ssd.*c%d.*GO.*r%d',iCondition,2), 'once')),prd.trialCat,'Uni',0));
                
                rtStopIErrorCCorrObsL{iCondition}   = cell2mat(obs.rtStopIErrorCCorr(iTrialCatStopL));
                rtStopIErrorCCorrPrdL{iCondition} 	= cell2mat(cellfun(@(x) reshape(x, [], 1), prd.rtStopIErrorCCorr(iTrialCatStopL), 'uni', false));
                rtStopICorrPrdL{iCondition}         = cell2mat(cellfun(@(x) reshape(x, [], 1), prd.rtStopICorr(iTrialCatStopL), 'uni', false));
                
                rtStopIErrorCCorrObsR{iCondition}   = cell2mat(obs.rtStopIErrorCCorr(iTrialCatStopR));
                rtStopIErrorCCorrPrdR{iCondition} 	= cell2mat(cellfun(@(x) reshape(x, [], 1), prd.rtStopIErrorCCorr(iTrialCatStopR), 'uni', false));
                rtStopICorrPrdR{iCondition}         = cell2mat(cellfun(@(x) reshape(x, [], 1), prd.rtStopICorr(iTrialCatStopR), 'uni', false));
                
                
                pStopIErrorCCorrObsL(iCondition)   = 1 - sum(obs.nStopIErrorCCorr(iTrialCatStopL)) / sum(obs.nStopIErrorCCorr(iTrialCatStopL) + obs.nStopIErrorCError(iTrialCatStopL));
                pStopIErrorCCorrPrdL(iCondition) 	= 1 - sum(prd.nStopIErrorCCorr(iTrialCatStopL)) / sum(prd.nStopIErrorCCorr(iTrialCatStopL) + prd.nStopIErrorCError(iTrialCatStopL));
                
                pStopIErrorCCorrObsR(iCondition)   = sum(obs.nStopIErrorCCorr(iTrialCatStopR)) / sum(obs.nStopIErrorCCorr(iTrialCatStopR) + obs.nStopIErrorCError(iTrialCatStopR));
                pStopIErrorCCorrPrdR(iCondition) 	= sum(prd.nStopIErrorCCorr(iTrialCatStopR)) / sum(prd.nStopIErrorCCorr(iTrialCatStopR) + prd.nStopIErrorCError(iTrialCatStopR));
                
                
                
            end
            
            plot(1:3,cellfun(@nanmean, rtStopIErrorCCorrObsL),'Color','r','Marker',stopIErrorCCorrMrkObs,'LineStyle',stopIErrorCCorrLnObs,'LineWidth',stopIErrorCCorrLnWidth);
            plot(1:3,cellfun(@mean, rtStopIErrorCCorrPrdL),'Color','r','Marker',stopIErrorCCorrMrkPrd,'LineStyle',stopIErrorCCorrLnPrd,'LineWidth',stopIErrorCCorrLnWidth);
            plot(1:3,cellfun(@mean, rtStopICorrPrdL),'Color','r','Marker',stopICorrMrkPrd,'LineStyle',stopICorrLnPrd,'LineWidth',stopICorrLnWidth);
            
            plot(6:-1:4,cellfun(@nanmean, rtStopIErrorCCorrObsR),'Color','r','Marker',stopIErrorCCorrMrkObs,'LineStyle',stopIErrorCCorrLnObs,'LineWidth',stopIErrorCCorrLnWidth);
            plot(6:-1:4,cellfun(@mean, rtStopIErrorCCorrPrdR),'Color','r','Marker',stopIErrorCCorrMrkPrd,'LineStyle',stopIErrorCCorrLnPrd,'LineWidth',stopIErrorCCorrLnWidth);
            plot(6:-1:4,cellfun(@mean, rtStopICorrPrdR),'Color','r','Marker',stopICorrMrkPrd,'LineStyle',stopICorrLnPrd,'LineWidth',stopICorrLnWidth);
            
            
            
            
            % Psychometric function
            % ================================================
            p(3,iModel).select();
            p(3,iModel).hold('on');
            %     p(2,iModel).title({sprintf('Model %d',model(iModel)), ...
            %                        modelStr{iModel}, ...
            %                        sprintf('\\chi^2 = %.1f',cost), ...
            %                        sprintf('BIC = %.1f',altCost)});
            
            
            
            pGoCCorrObsL     = 1 - obs.pGoCCorr(trialCatGoL);
            pGoCCorrPrdL     = 1 - prd.pGoCCorr(trialCatGoL);
            
            %                         plot(1:3,pGoCCorrPrdL,'Color','k','Marker',goCCorrMrkPrd,'LineStyle',goCCorrLnPrd,'LineWidth',goCCorrLnWidth);
            
            pGoCCorrObsR     = obs.pGoCCorr(trialCatGoR);
            pGoCCorrPrdR     = prd.pGoCCorr(trialCatGoR);
            
            plot(1:3,pGoCCorrObsL,'Color','k','Marker',goCCorrMrkObs,'LineStyle',goCCorrLnObs,'LineWidth',goCCorrLnWidth);
            plot(6:-1:4,pGoCCorrObsR,'Color','k','Marker',goCCorrMrkObs,'LineStyle',goCCorrLnObs,'LineWidth',goCCorrLnWidth);
            plot(1:6,[pGoCCorrPrdL; flipud(pGoCCorrPrdR)],'Color','k','Marker',goCCorrMrkPrd,'LineStyle',goCCorrLnPrd,'LineWidth',goCCorrLnWidth);
            
            plot(1:3,pStopIErrorCCorrObsL,'Color','r','Marker',stopIErrorCCorrMrkObs,'LineStyle',stopIErrorCCorrLnObs,'LineWidth',stopIErrorCCorrLnWidth);
            plot(6:-1:4,pStopIErrorCCorrObsR,'Color','r','Marker',stopIErrorCCorrMrkObs,'LineStyle',stopIErrorCCorrLnObs,'LineWidth',stopIErrorCCorrLnWidth);
            plot(1:6,[pStopIErrorCCorrPrdL, fliplr(pStopIErrorCCorrPrdR)],'Color','r','Marker',stopIErrorCCorrMrkPrd,'LineStyle',stopIErrorCCorrLnPrd,'LineWidth',stopIErrorCCorrLnWidth);
            
        end
    end

end
function plot_fits(subject,model,architecture,dt,trialVar,optimScope,fileStr,savePlot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. PROCESS INPUTS & SPECIFY VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% close all

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
        
        % Set up the figure and panels
        set_figure({1024,574,'pixels'},{'USLetter','landscape'},{'Helvetica',18});
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
            
            % Load model-specific SAM
            ds = dataset('File',fullfile(sprintf(rootDir,subject(iSubject),dt,trialVarStr,architecture{iArchitecture}),sprintf(nameFVal,optimScope,model(iModel))));
            
            % Load SAM
            load(fullfile(sprintf(rootDir,subject(iSubject),dt,trialVarStr,architecture{iArchitecture}),ds.FileNameSAM{1}),'SAM');
            
            
            % Extract optimized X
            iBestX = cell2mat(cellfun(@(in1) ~isempty(regexp(in1,'^BestX.*', 'once')),ds.Properties.VarNames,'Uni',0));
            X = double(ds(1,iBestX))
            
            % Get model predictions and costs
            [cost,altCost,prd] = sam_cost(X,SAM);
            
            % Plot observations and predictions
            % =========================================================================
%                     sam_plot(SAM,prd);
            
            % Plot it
            plotit(SAM,prd,p,model,iModel,modelStr,cost,altCost);
            
        end
    end
end


    function plotit(SAM,prd,p,model,iModel,modelStr,cost,altCost)
        
        % Specify colors and line properties
        cndClr                  = {[0 0.75 1],[1 0 0.5],[0.25 0 0.5]};
        
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
        
        % SSD index
        iSsd = 3;
        
        % Specificy observations
        obs = SAM.optim.obs;
        
        % Inhibition function
        p(2,iModel).select();
        p(2,iModel).hold('on');
        %     p(2,iModel).title({sprintf('Model %d',model(iModel)), ...
        %                        modelStr{iModel}, ...
        %                        sprintf('\\chi^2 = %.1f',cost), ...
        %                        sprintf('BIC = %.1f',altCost)});
        p(2,iModel).title({sprintf('Model %d',model(iModel)), ...
            sprintf('\\chi^2 = %.1f',cost), ...
            sprintf('BIC = %.1f',altCost)});
        
        for iRsp = 1:2
            for iCnd = 1:3
                
                % Identify the relevant rows in the dataset array
                iTrialCat = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('stopTrial.*GO:c%d',iCnd), 'once')),prd.trialCat,'Uni',0)));
                iTrialCat = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('stopTrial.*GO.*r%d.*c%d',iRsp,iCnd), 'once')),prd.trialCat,'Uni',0)));
                
                % Get the SSDs and observed and predicted response rates
                ssd = obs.ssd(iTrialCat);
                respRateObs = 1-obs.pStopICorr(iTrialCat);
                respRatePrd = 1-prd.pStopICorr(iTrialCat);
                
                % Plot them
                plot(ssd,respRateObs,'Color',cndClr{iCnd},'Marker',stopIErrorCCorrMrkObs,'LineStyle',stopIErrorCCorrLnObs);
                plot(ssd,respRatePrd,'Color',cndClr{iCnd},'Marker',stopIErrorCCorrMrkPrd,'LineStyle',stopIErrorCCorrLnPrd);
                
            end
        end
        
        % Set axes
        set(gca,'XLim',[0 750], ...
            'XTick',0:250:750, ...
            'YLim',[0 1], ...
            'YTick',0:0.2:1)
        
        
        % RT distribution Go trials
        p(3,iModel).select();
        p(3,iModel).hold('on');
        
        for iRsp = 1:2
            for iCnd = 1:3
                
                % Identify the relevant rows in the dataset array
                %         iTrialCatGo = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('goTrial.*GO:c%d',iCnd), 'once')),prd.trialCat,'Uni',0)));
                iTrialCatGo = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('goTrial.*GO.*r%d.*c%d',iRsp,iCnd), 'once')),prd.trialCat,'Uni',0)));
                if isempty(iTrialCatGo)
                    iTrialCatGo = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('goTrial_c%d.*',iCnd), 'once')),prd.trialCat,'Uni',0)));
                end
                
                % GoCCorr trials
                % -----------------------------------------------------------------
                rtGoCCorrObs     = obs.rtQGoCCorr{iTrialCatGo};
                cumPGoCCorrObs   = obs.cumProbGoCCorr{iTrialCatGo};
                
                rtGoCCorrPrd     = prd.rtGoCCorr{iTrialCatGo};
                cumPGoCCorrPrd   = cmtb_edf(prd.rtGoCCorr{iTrialCatGo}(:),prd.rtGoCCorr{iTrialCatGo}(:));
                
                % Plot it
                plot(rtGoCCorrObs,cumPGoCCorrObs,'Color',cndClr{iCnd},'Marker',goCCorrMrkObs,'LineStyle',goCCorrLnObs,'LineWidth',goCCorrLnWidth);
                plot(rtGoCCorrPrd,cumPGoCCorrPrd,'Color',cndClr{iCnd},'Marker',goCCorrMrkPrd,'LineStyle',goCCorrLnPrd,'LineWidth',goCCorrLnWidth);
                
                % GoCError trials
                % -----------------------------------------------------------------
                if ~isempty(obs.rtQGoCError{iTrialCatGo}) && obs.nGoCError(iTrialCatGo) > 10
                    
                    rtGoCErrorObs     = obs.rtQGoCError{iTrialCatGo};
                    cumPGoCErrorObs   = obs.cumProbGoCError{iTrialCatGo};
                    
                    rtGoCErrorPrd     = prd.rtGoCError{iTrialCatGo};
                    cumPGoCErrorPrd   = cmtb_edf(prd.rtGoCError{iTrialCatGo}(:),prd.rtGoCError{iTrialCatGo}(:));
                    
                    % Plot it
                    plot(rtGoCErrorObs,cumPGoCErrorObs,'Color',cndClr{iCnd},'Marker',goCErrorMrkObs,'LineStyle',goCErrorLnObs,'LineWidth',goCErrorLnWidth);
                    plot(rtGoCErrorPrd,cumPGoCErrorPrd,'Color',cndClr{iCnd},'Marker',goCErrorMrkPrd,'LineStyle',goCErrorLnPrd,'LineWidth',goCErrorLnWidth);
                    
                end
            end
        end
        
        % Set axes
        set(gca,'XLim',[0 1000], ...
            'XTick',0:500:1000, ...
            'YLim',[0 1], ...
            'YTick',0:0.2:1)
        
%         % RT distribution Stop trials
%         p(4,iModel).select();
%         p(4,iModel).hold('on');
%         
%         for iRsp = 1:2
%             for iCnd = 1:3
%                 
%                 % Identify the relevant rows in the dataset array
% %                 iTrialCatStop = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('stopTrial.*ssd%d.*GO:c%d',iSsd,iCnd), 'once')),prd.trialCat,'Uni',0)));
%                 iTrialCatStop = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('stopTrial.*ssd%d.*GO.*r%d.*c%d',iSsd,iRsp,iCnd), 'once')),prd.trialCat,'Uni',0)));
%                 
%                 if isempty(iTrialCatStop)
%                     iTrialCatStop = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('stopTrial.*ssd%d.*c%d.*',iSsd,iCnd), 'once')),prd.trialCat,'Uni',0)));
%                 end
%                 
%                 % StopIErrorCCorr trials
%                 % -----------------------------------------------------------------
%                 rtStopIErrorCCorrObs   = obs.rtQStopIErrorCCorr{iTrialCatStop};
%                 cumPStopIErrorCCorrObs = obs.cumProbStopIErrorCCorr{iTrialCatStop};
%                 
%                 rtStopIErrorCCorrPrd     = prd.rtStopIErrorCCorr{iTrialCatStop};
%                 cumPStopIErrorCCorrPrd   = cmtb_edf(prd.rtStopIErrorCCorr{iTrialCatStop}(:),prd.rtStopIErrorCCorr{iTrialCatStop}(:));
%                 
%                 % Plot it
%                 plot(rtStopIErrorCCorrObs,cumPStopIErrorCCorrObs,'Color',cndClr{iCnd},'Marker',stopIErrorCCorrMrkObs,'LineStyle',stopIErrorCCorrLnObs,'LineWidth',stopIErrorCCorrLnWidth);
%                 plot(rtStopIErrorCCorrPrd,cumPStopIErrorCCorrPrd,'Color',cndClr{iCnd},'Marker',stopIErrorCCorrMrkPrd,'LineStyle',stopIErrorCCorrLnPrd,'LineWidth',stopIErrorCCorrLnWidth);
%                 
%                 % StopICorr trials
%                 % -----------------------------------------------------------------
%                 rtStopICorrPrd     = prd.rtStopICorr{iTrialCatStop};
%                 cumPStopICorrPrd   = cmtb_edf(prd.rtStopICorr{iTrialCatStop}(:),prd.rtStopICorr{iTrialCatStop}(:));
%                 
%                 % Plot it
%                 plot(rtStopICorrPrd,cumPStopICorrPrd,'Color',cndClr{iCnd},'Marker',stopICorrMrkPrd,'LineStyle',stopICorrLnPrd,'LineWidth',stopICorrLnWidth);
%                 
%             end
%         end
        
        % Set axes
        set(gca,'XLim',[100 600], ...
            'XTick',100:100:600, ...
            'YLim',[0 1], ...
            'YTick',0:0.2:1)
        
    end

end
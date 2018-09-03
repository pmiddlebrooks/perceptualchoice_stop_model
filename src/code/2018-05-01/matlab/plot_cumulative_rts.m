function [obs, prd, obsRts, prdRts] = plot_cumulative_rts(subject,model,architecture,dt,trialVar,optimScope,fileStr,responseSide, coherence, accuracy, defective, addData, savePlot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. PROCESS INPUTS & SPECIFY VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all
nSim = 5000;
% nSim = 50;
% 1.1. Process inputs
% =========================================================================

switch subject
    case 1
        subjectName = 'broca';
    case 2
        subjectName = 'xena';
    case 3
        subjectName = 'joule';
end

if trialVar
    trialVarStr = 'trialvar';
else
    trialVarStr = 'notrialvar';
end

% 1.2. Specify dynamic variables
% =========================================================================

% Numbers
nModel              = numel(model);

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


% Display progress
fprintf(1,'Working on subject %d, architecture %s, and model %d. \n',subject,architecture,model);

% Load model-specific fits with cost function values

ds = dataset('File',fullfile(sprintf(rootDir,subject,dt,trialVarStr,architecture),sprintf(nameFVal,optimScope,model)));
%             ds = dataset('File',fullfile(sprintf(rootDir,subject,dt,trialVarStr,architecture),sprintf(nameFVal,optimScope,modelNum)));

% Load SAM
load(fullfile(sprintf(rootDir,subject,dt,trialVarStr,architecture),ds.FileNameSAM{1}),'SAM');

% Extract optimized X
iBestX = cell2mat(cellfun(@(in1) ~isempty(regexp(in1,'^BestX.*', 'once')),ds.Properties.VarNames,'Uni',0));
X = double(ds(1,iBestX));

% Alter the simulations for efficiency
SAM.sim.n = nSim;

% Specificy observations
obs = SAM.optim.obs;

ssdList = sort(unique(obs.ssd(~isnan(obs.ssd))));
ssdList(ssdList == 0) = [];
stopColorList = flipud(linspace(0.2, 0.9, length(ssdList))');

% Get model predictions and costs
[cost,altCost,prd] = sam_cost(X,SAM);

noSsd = prd.ssd == 0;
obs = obs(~noSsd, :);
prd = prd(~noSsd, :);


% Plot observations and predictions
% =========================================================================
%                     sam_plot(SAM,prd);


% Plot it
plotit(SAM,prd,model,modelStr,cost,altCost, responseSide, accuracy, defective);

% Export lists of RTs, conditions
obsRts = table(tableObsCnd, tableObsRsp, tableObsSsd, tableObsRt, 'VariableNames',{'Coherence','Response','Ssd','RT'});
prdRts = table(tablePrdCnd, tablePrdRsp, tablePrdSsd, tablePrdRt, 'VariableNames',{'Coherence','Response','Ssd','RT'});


if savePlot
    saveDir             = fullfile(sprintf(fileStr.result,subject,dt,trialVarStr,architecture));
    if exist(saveDir,'dir') ~= 7
        mkdir(saveDir)
    end
    fileName = sprintf('%s_Cumulative_rts_Respond_%s_Coherence_%s_Accuracy_%s', addData, responseSide, coherence, accuracy);
    print(gcf, fullfile(saveDir, fileName),'-dpdf', '-r300')
    
    % Export data
    writetable(obsRts,fullfile(saveDir,[subjectName,'_',addData,'_obsRts.csv']),'Delimiter',',')
    writetable(prdRts,fullfile(saveDir,[subjectName,'_',addData,'_prdRts.csv']),'Delimiter',',')
end
clear SAM








    function plotit(SAM,prd,model,modelStr,cost,altCost, responseSide, accuracy, defective)
        
        % Specify colors and line properties
        mapShades           	= [.1 .25 .4 .6 .75 .9];
        cndClr                  = ccm_colormap(mapShades);
        cndClr(1+length(mapShades)/2 : end, :) = flipud(cndClr(1+length(mapShades)/2 : end, :));
        
        
        switch responseSide
            case 'collapse'
                responseArray = 1;
            case 'both'
                responseArray = [1 2];
            case 'left'
                responseArray = 1;
            case 'right'
                responseArray = 2;
            otherwise
                error('responseSide need to be left, right, or both');
        end
        switch coherence
            case 'collapse'
                conditionArray = 1;
            case 'all'
                switch subject
                    case 1
                        if size(obs, 1) == 26
                            conditionArray = 1:2;
                        else
                            conditionArray = 1:3;
                        end
                    case 2
                        conditionArray = 1:3;
                    case 3
                        conditionArray = 1:2;
                end
        end
        
        % Set up the figure and panels
        %         [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_figure;
        standard_figure(1,1,'landscape');
        
        %                 set_figure({1024,574,'pixels'},{'USLetter','landscape'},{'Helvetica',18});
        p = panel();
        % nPlotCol = max(3, nModel);
        nPlotCol = 2;
        if length(conditionArray) == 2
            p.pack({.05, .45, .45, .05}, num2cell(repmat(1/nPlotCol,1,nPlotCol)));
        elseif length(conditionArray) == 3
            p.pack({.05, .3 .3 .3, .05}, num2cell(repmat(1/nPlotCol,1,nPlotCol)));
        elseif length(conditionArray) == 1
            p.pack({.05, .9, .05}, num2cell(repmat(1/nPlotCol,1,nPlotCol)));
        end
        annotation('textbox', [0 0.9 1 0.1], ...
            'String', sprintf('Subject %d, architecture %s',subject,architecture), ...
            'Interpreter','none', ...
            'EdgeColor', 'none', ...
            'HorizontalAlignment', 'center');
        
        
        
        
        goCCorrMrkObs           = 'none';
        goCCorrLnObs            = '-';
        goCCorrLnWidthObs          = 3;
        goCCorrMrkPrd           = 'none';
        goCCorrLnPrd            = '-';
        goCCorrLnWidthPrd          = 7;
        
        goCErrorMrkObs          = '^';
        goCErrorLnObs           = 'none';
        goCErrorMrkPrd          = 'none';
        goCErrorLnPrd           = '-.';
        goCErrorLnWidth         = 1;
        
        stopIErrorCCorrMrkObs   = 'none';
        stopIErrorCCorrLnObs    = '-';
        stopIErrorCCorrLnWidthObs  = 2;
        stopIErrorCCorrMrkPrd   = 'none';
        stopIErrorCCorrLnPrd    = '-';
        stopIErrorCCorrLnWidthPrd  = 5;
        
        stopICorrMrkPrd         = 'none';
        stopICorrLnPrd          = ':';
        stopICorrLnWidth        = 1;
        
        
        tableObsCnd = [];
        tableObsRsp = [];
        tableObsSsd = [];
        tableObsRt = [];
        
        tablePrdCnd = [];
        tablePrdRsp = [];
        tablePrdSsd = [];
        tablePrdRt = [];
        
        % RT distribution Go trials
        % ================================================
        
        if strcmp(optimScope, 'go') || strcmp(optimScope, 'all')
            
            %             p(2,iModel).title({sprintf('Model %d',model(iModel)), ...
            %                 sprintf('\\chi^2 = %.1f',cost), ...
            %                 sprintf('BIC = %.1f',altCost)});
            
            for iResponse = 1:length(responseArray)
                for iCondition = 1:length(conditionArray)
                    
                    iRsp = responseArray(iResponse);
                    iCnd = conditionArray(iCondition);
                    
                    p(iCondition+1, iResponse).select();
                    p(iCondition+1, iResponse).title({sprintf('Resp: %d, Coh: %d', iRsp, iCnd)});
                    p(iCondition+1, iResponse).hold('on');
                    
                    
                    if strcmp(coherence, 'collapse') && strcmp(responseSide, 'collapse')
                        iTrialCatGo = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('goTrial.*'), 'once')),prd.trialCat,'Uni',0)));
                    elseif strcmp(coherence, 'collapse') && strcmp(responseSide, 'both')
                        iTrialCatGo = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('goTrial.*GO.*r%d',iRsp), 'once')),prd.trialCat,'Uni',0)));
                    elseif strcmp(coherence, 'both') && strcmp(responseSide, 'collapse')
                        iTrialCatGo = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('goTrial.*c%d.*GO.*',iCnd), 'once')),prd.trialCat,'Uni',0)));
                    else
                        % Identify the relevant rows in the dataset array
                        iTrialCatGo = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('goTrial.*c%d.*GO.*r%d',iCnd,iRsp), 'once')),prd.trialCat,'Uni',0)));
                        if isempty(iTrialCatGo)
                            iTrialCatGo = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('goTrial_c%d.*',iCnd), 'once')),prd.trialCat,'Uni',0)));
                        end
                    end
                    
                    % GoCCorr trials
                    % -----------------------------------------------------------------
                    if strcmp(accuracy, 'both') || strcmp(accuracy, 'correct')
                        obs.rtGoCCorr = cellfun(@(x) reshape(x, length(x), 1), obs.rtGoCCorr, 'uni', false);
                        rtGoCCorrObs     = sort(cell2mat(obs.rtGoCCorr(iTrialCatGo)));
                        cumPGoCCorrObs   = cmtb_edf(rtGoCCorrObs,rtGoCCorrObs);
                        
                        prd.rtGoCCorr = cellfun(@(x) reshape(x, length(x), 1), prd.rtGoCCorr, 'uni', false);
                        rtGoCCorrPrd     = sort(cell2mat(prd.rtGoCCorr(iTrialCatGo)));
                        cumPGoCCorrPrd   = cmtb_edf(rtGoCCorrPrd,rtGoCCorrPrd);
                        
                        if defective
                            cumPGoCCorrPrd   = cumPGoCCorrPrd * prd.pGoCCorr(iTrialCatGo);
                            cumPGoCCorrObs   = obs.cumProbDefectiveGoCCorr{iTrialCatGo};
                        end
                        
                        % Plot it
                        iColor = cndClr((length(conditionArray) * (iRsp -1) + iCondition), :);
                        iColor = 'k';
                        plot(rtGoCCorrObs,cumPGoCCorrObs,'Color',iColor,'Marker',goCCorrMrkObs,'LineStyle',goCCorrLnObs,'LineWidth',goCCorrLnWidthObs);
                        plot(rtGoCCorrPrd,cumPGoCCorrPrd,'Color',iColor,'Marker',goCCorrMrkPrd,'LineStyle',goCCorrLnPrd,'LineWidth',goCCorrLnWidthPrd);
                        
                        % Build a table for export to CSV
                        tableObsRt = [tableObsRt; rtGoCCorrObs];
                        tableObsCnd = [tableObsCnd; repmat(iCnd, length(rtGoCCorrObs), 1)];
                        tableObsRsp = [tableObsRsp; repmat(iRsp, length(rtGoCCorrObs), 1)];
                        tableObsSsd = [tableObsSsd; nan(length(rtGoCCorrObs), 1)];
                        
                        tablePrdRt = [tablePrdRt; rtGoCCorrPrd];
                        tablePrdCnd = [tablePrdCnd; repmat(iCnd, length(rtGoCCorrPrd), 1)];
                        tablePrdRsp = [tablePrdRsp; repmat(iRsp, length(rtGoCCorrPrd), 1)];
                        tablePrdSsd = [tablePrdSsd; nan(length(rtGoCCorrPrd), 1)];
                        
                        
                        
                    end
                    
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
                end
            end
            
            
        end
        
        
        
        
        
        % RT distribution Stop trials
        % ================================================
        
        if strcmp(optimScope, 'stop') || strcmp(optimScope, 'all')
            
            
            % Get rid of NaNs in obs.cumProbStopIErrorCCorr
            for j = 1 : size(obs, 1)
                if isnan(obs.cumProbStopIErrorCCorr{j})
                    obs.cumProbStopIErrorCCorr{j} = [];
                end
            end
            for iResponse = 1:length(responseArray)
                for iCondition = 1:length(conditionArray)
                    p(iCondition+1, iResponse).select();
                    p(iCondition+1, iResponse).hold('on');
                    
                    
                    for iSsd = 1 : length(ssdList)
                        iRsp = responseArray(iResponse);
                        iCnd = conditionArray(iCondition);
                        
                        % Identify the relevant rows in the dataset array
                        if strcmp(coherence, 'collapse') && strcmp(responseSide, 'collapse')
                            iTrialCatStop = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('stopTrial.*ssd%d.*',iSsd), 'once')),prd.trialCat,'Uni',0)));
                        elseif strcmp(coherence, 'collapse') && strcmp(responseSide, 'both')
                            iTrialCatStop = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('stopTrial.*ssd%d.*GO.*r%d',iSsd,iRsp), 'once')),prd.trialCat,'Uni',0)));
                        elseif strcmp(coherence, 'both') && strcmp(responseSide, 'collapse')
                            iTrialCatStop = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('stopTrial.*ssd%d.*c%d.*GO.*',iSsd,iCnd), 'once')),prd.trialCat,'Uni',0)));
                        else
                            iTrialCatStop = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('stopTrial.*ssd%d.*GO.*r%d.*c%d',iSsd,iRsp,iCnd), 'once')),prd.trialCat,'Uni',0)));
                            
                            if isempty(iTrialCatStop)
                                iTrialCatStop = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('stopTrial.*ssd%d.*c%d.*r%d.*',iSsd,iCnd,iRsp), 'once')),prd.trialCat,'Uni',0)));
                            end
                        end
                        
                        % StopIErrorCCorr trials
                        % -----------------------------------------------------------------
                        obs.rtStopIErrorCCorr = cellfun(@(x) reshape(x, length(x), 1), obs.rtStopIErrorCCorr, 'uni', false);
                        rtStopIErrorCCorrObs   = sort(cell2mat(obs.rtStopIErrorCCorr(iTrialCatStop)));
                        if ~isempty(rtStopIErrorCCorrObs)
                            cumPStopIErrorCCorrObs = cmtb_edf(rtStopIErrorCCorrObs,rtStopIErrorCCorrObs);
                            
                            prd.rtStopIErrorCCorr = cellfun(@(x) reshape(x, length(x), 1), prd.rtStopIErrorCCorr, 'uni', false);
                            rtStopIErrorCCorrPrd     = sort(cell2mat(prd.rtStopIErrorCCorr(iTrialCatStop)));
                            if ~isempty(rtStopIErrorCCorrPrd)
                                cumPStopIErrorCCorrPrd   = cmtb_edf(rtStopIErrorCCorrPrd,rtStopIErrorCCorrPrd);
                                
                                if defective
                                    cumPStopIErrorCCorrPrd   = cumPStopIErrorCCorrPrd * prd.pStopIErrorCCorr(iTrialCatStop);
                                    %                         if sum(~isnan(cumPStopIErrorCCorrPrd)) == 0
                                    %                             cumPStopIErrorCCorrPrd = []
                                    %                         end
                                    cumPStopIErrorCCorrObs   = cumPStopIErrorCCorrObs * obs.pStopIErrorCCorr(iTrialCatStop);
                                    %                         cumPStopIErrorCCorrObs   = obs.cumProbDefectiveStopIErrorCCorr{iTrialCatGo}
                                end
                                
                                % Plot it
                                %                     iColor = cndClr((length(conditionArray) * (iRsp -1) + iCondition), :);
                                iColor = repmat(stopColorList(iSsd),1,3);
                                plot(rtStopIErrorCCorrObs,cumPStopIErrorCCorrObs,'Color',iColor,'Marker',stopIErrorCCorrMrkObs,'LineStyle',stopIErrorCCorrLnObs,'LineWidth',stopIErrorCCorrLnWidthObs);
                                plot(rtStopIErrorCCorrPrd,cumPStopIErrorCCorrPrd,'Color',iColor,'Marker',stopIErrorCCorrMrkPrd,'LineStyle',stopIErrorCCorrLnPrd,'LineWidth',stopIErrorCCorrLnWidthPrd);
                                
                                % Build a table for export to CSV
                                tableObsRt = [tableObsRt; rtStopIErrorCCorrObs];
                                tableObsCnd = [tableObsCnd; repmat(iCnd, length(rtStopIErrorCCorrObs), 1)];
                                tableObsRsp = [tableObsRsp; repmat(iRsp, length(rtStopIErrorCCorrObs), 1)];
                                tableObsSsd = [tableObsSsd; repmat(ssdList(iSsd), length(rtStopIErrorCCorrObs), 1)];
                                
                                tablePrdRt = [tablePrdRt; rtStopIErrorCCorrPrd];
                                tablePrdCnd = [tablePrdCnd; repmat(iCnd, length(rtStopIErrorCCorrPrd), 1)];
                                tablePrdRsp = [tablePrdRsp; repmat(iRsp, length(rtStopIErrorCCorrPrd), 1)];
                                tablePrdSsd = [tablePrdSsd; repmat(ssdList(iSsd), length(rtStopIErrorCCorrPrd), 1)];
                                
                            end
                        end
                        %                     % StopICorr trials
                        %                     % -----------------------------------------------------------------
                        %                     rtStopICorrPrd     = prd.rtStopICorr{iTrialCatStop};
                        %                     cumPStopICorrPrd   = cmtb_edf(prd.rtStopICorr{iTrialCatStop}(:),prd.rtStopICorr{iTrialCatStop}(:));
                        %
                        %                     % Plot it
                        %                     iColor = 'r';%cndClr((length(conditionArray) * (iRsp -1) + iCondition), :);
                        %                     plot(rtStopICorrPrd,cumPStopICorrPrd,'Color',iColor,'Marker',stopICorrMrkPrd,'LineStyle',stopICorrLnPrd,'LineWidth',stopICorrLnWidth);
                        
                    end
                    % Set axes
                    switch subject
                        case 1
                            set(gca,'XLim',[200 800], ...
                                'XTick',100:100:800, ...
                                'YLim',[0 1], ...
                                'YTick',0:0.2:1)
                        case 2
                            set(gca,'XLim',[200 500], ...
                                'XTick',100:100:500, ...
                                'YLim',[0 1], ...
                                'YTick',0:0.2:1)
                        case 3
                            set(gca,'XLim',[200 700], ...
                                'XTick',100:100:700, ...
                                'YLim',[0 1], ...
                                'YTick',0:0.2:1)
                        otherwise
                            error('Need to add axes limits for subject')
                    end
                    
                end
            end
            
            
            
            
            
        end
    end

end
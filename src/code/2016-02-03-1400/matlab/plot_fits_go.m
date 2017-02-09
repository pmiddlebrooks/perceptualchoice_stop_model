function plot_fits_go(subject,model,architecture,dt,trialVar,fileStr,responseSide, accuracy, savePlot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. PROCESS INPUTS & SPECIFY VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

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
optimScope          = 'go';
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
        figureHandle = 19 + iSubject * iArchitecture;
        
        % Set up the figure and panels
        %         [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_figure;
        standard_figure(1,1,'landscape', figureHandle);
        
        %                 set_figure({1024,574,'pixels'},{'USLetter','landscape'},{'Helvetica',18});
        p = panel();
        p.pack({.05 .475 .475}, num2cell(repmat(1/nModel,1,nModel)));
        
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
            
            % Print out otpimized X
            fprintf('Optimal parameters: \n')
            disp(X)
            
            % Get model predictions and costs
            [cost,altCost,prd] = sam_cost(X,SAM);
            
            % Plot observations and predictions
            % =========================================================================
            %                     sam_plot(SAM,prd);
            
            % Plot it
            plotit(SAM,prd,p,model,iModel,modelStr,cost,altCost, responseSide, accuracy);
            
            if savePlot
                saveDir             = fullfile(sprintf(fileStr.result,subject(iSubject),dt,trialVarStr,architecture{iArchitecture}));
                if exist(saveDir,'dir') ~= 7
                    mkdir(saveDir)
                end
                fileName = sprintf('GO_Respond_%s_Accuracy_%s', responseSide, accuracy);
                print(gcf, fullfile(saveDir, fileName),'-dpdf', '-r300')
                print(gcf, fullfile(saveDir, fileName),'-dpng')
            end
        end
    end
end


    function plotit(SAM,prd,p,model,iModel,modelStr,cost,altCost, responseSide, accuracy)
        
        % Specify colors and line properties
        cndClr                  = ccm_colormap([.1 .25 .4 .6 .75 .9]);
        cndClr(4:end,:) = flipud(cndClr(4:end,:));  % make the color map go [response 1 easy to hard : response2 easy to hard]
        %         cndClr                  = {[0 0.75 1],[1 0 0.5],[0.25 0 0.5]};
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
        
        
        % Specificy observations
        obs = SAM.optim.obs;
        
        
        
        % RT distribution Go trials: Correct Choices
        % ================================================
        if strcmp(accuracy, 'both') || strcmp(accuracy, 'correct')
            
            p(2,iModel).select();
            p(2,iModel).hold('on');
            p(2,iModel).title({sprintf('Model %d',model(iModel)), ...
                sprintf('\\chi^2 = %.1f',cost), ...
                sprintf('BIC = %.1f',altCost)});
            
            for iResponse = 1:length(responseArray)
                for iCondition = 1:length(conditionArray);
                    
                    iRsp = responseArray(iResponse);
                    iCnd = conditionArray(iCondition);
                    
                    % Identify the relevant rows in the dataset array
                    %         iTrialCatGo = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('goTrial.*GO:c%d',iCnd), 'once')),prd.trialCat,'Uni',0)));
                    iTrialCatGo = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('goTrial.*GO.*r%d.*c%d',iRsp,iCnd), 'once')),prd.trialCat,'Uni',0)));
                    if isempty(iTrialCatGo)
                        iTrialCatGo = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('goTrial_c%d.*',iCnd), 'once')),prd.trialCat,'Uni',0)));
                    end
                    
                    % GoCCorr trials
                    % -----------------------------------------------------------------
                    if strcmp(accuracy, 'both') || strcmp(accuracy, 'correct')
                        rtGoCCorrObs     = obs.rtQGoCCorr{iTrialCatGo};
                        cumPGoCCorrObs   = obs.cumProbGoCCorr{iTrialCatGo};
                        
                        rtGoCCorrPrd     = prd.rtGoCCorr{iTrialCatGo};
                        cumPGoCCorrPrd   = cmtb_edf(prd.rtGoCCorr{iTrialCatGo}(:),prd.rtGoCCorr{iTrialCatGo}(:));
                        
                        % Plot it
                        iColor = cndClr((length(conditionArray) * (iRsp -1) + iCondition), :);
                        plot(rtGoCCorrObs,cumPGoCCorrObs,'Color',iColor,'Marker',goCCorrMrkObs,'LineStyle',goCCorrLnObs,'LineWidth',goCCorrLnWidth);
                        plot(rtGoCCorrPrd,cumPGoCCorrPrd,'Color',iColor,'Marker',goCCorrMrkPrd,'LineStyle',goCCorrLnPrd,'LineWidth',goCCorrLnWidth);
                        
                        
                    end
                end
                
                % Set axes
                switch subject(iSubject)
                    case 1
                        set(gca,'XLim',[200 700], ...
                            'XTick',100:100:700, ...
                            'YLim',[0 1], ...
                            'YTick',0:0.2:1)
                    case 2
                        set(gca,'XLim',[150 550], ...
                            'XTick',150:100:550, ...
                            'YLim',[0 1], ...
                            'YTick',0:0.2:1)
                    otherwise
                        error('Need to add axes limits for subject')
                end
                
            end
        end
        
        
        % RT distribution Go trials: Error Choices
        % ================================================
        if  strcmp(accuracy, 'both') || strcmp(accuracy, 'error')
            
            p(3,iModel).select();
            p(3,iModel).hold('on');
            
            for iResponse = 1:length(responseArray)
                for iCondition = 1:length(conditionArray);
                    
                    iRsp = responseArray(iResponse);
                    iCnd = conditionArray(iCondition);
                    
                    % Identify the relevant rows in the dataset array
                    %         iTrialCatGo = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('goTrial.*GO:c%d',iCnd), 'once')),prd.trialCat,'Uni',0)));
                    iTrialCatGo = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('goTrial.*GO.*r%d.*c%d',iRsp,iCnd), 'once')),prd.trialCat,'Uni',0)));
                    if isempty(iTrialCatGo)
                        iTrialCatGo = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('goTrial_c%d.*',iCnd), 'once')),prd.trialCat,'Uni',0)));
                    end
                    
                    
                    % GoCError trials
                    % -----------------------------------------------------------------
                    if ~isempty(obs.rtQGoCError{iTrialCatGo}) && obs.nGoCError(iTrialCatGo) > 10
                        
                        rtGoCErrorObs     = obs.rtQGoCError{iTrialCatGo};
                        cumPGoCErrorObs   = obs.cumProbGoCError{iTrialCatGo};
                        
                        rtGoCErrorPrd     = prd.rtGoCError{iTrialCatGo};
                        cumPGoCErrorPrd   = cmtb_edf(prd.rtGoCError{iTrialCatGo}(:),prd.rtGoCError{iTrialCatGo}(:));
                        
                        % Plot it
                        iColor = cndClr((length(conditionArray) * (iRsp -1) + iCondition), :);
                        plot(rtGoCErrorObs,cumPGoCErrorObs,'Color',iColor,'Marker',goCErrorMrkObs,'LineStyle',goCErrorLnObs,'LineWidth',goCErrorLnWidth);
                        plot(rtGoCErrorPrd,cumPGoCErrorPrd,'Color',iColor,'Marker',goCErrorMrkPrd,'LineStyle',goCErrorLnPrd,'LineWidth',goCErrorLnWidth);
                        
                    end
                end
            end
            
            % Set axes
            switch subject(iSubject)
                case 1
                    set(gca,'XLim',[200 700], ...
                        'XTick',100:100:700, ...
                        'YLim',[0 1], ...
                        'YTick',0:0.2:1)
                case 2
                    set(gca,'XLim',[150 550], ...
                        'XTick',150:100:550, ...
                        'YLim',[0 1], ...
                        'YTick',0:0.2:1)
                otherwise
                    error('Need to add axes limits for subject')
            end
            
        end
        
    end
end


function plot_simulation(subject,model,architecture,dt,trialVar,fileStr, savePlot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. PROCESS INPUTS & SPECIFY VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1.0 hard-code which conditions and responses to plot for now
close all
conditionArray       = 2;
responseArray   = 1;
accuracy        = 'correct';
responseSide    = 'right';

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


% 1.4 Hard code some plot parameters
% =========================================================================
rtLim = 400;  % how far out do draw starting point and threshold
zcLim = 150;  % how far out do draw starting point and threshold
modelLnStyle = {'-','-','-'};
modelLnWidth = [3 3 3];
modelLnClr = [0 0 0; .2 .8 .4; ([139 69 19]/255)];
z0LnWidth = 2;
zcLnWidth = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Plot the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop over subjects, architectures, and models
warning('off','stats:dataset:ModifiedVarnames')
for iSubject = 1:nSubject
    figureHandle = 19 + iSubject;
    
    % Set up the figure and panels
    %         [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_figure;
    standard_figure(1,1,'landscape', figureHandle);
    
    %                 set_figure({1024,574,'pixels'},{'USLetter','landscape'},{'Helvetica',18});
    p = panel();
    p.pack({.05 .475 .475}, num2cell(repmat(1/nModel,1,nModel)));
    
    
    for iArchitecture = 1:nArchitecture
        p(2,iArchitecture).select();
        p(2,iArchitecture).hold('on');
        p(2,iArchitecture).title({sprintf('Architecture %s',architecture{iArchitecture})});
        
        for iModel = 1:nModel
            modelNum = model(iModel);
            
            
            % Display progress
            fprintf(1,'Working on subject %d, architecture %s, and model %d. \n',subject(iSubject),architecture{iArchitecture},modelNum);
            
            % Load model-specific fits with cost function values
            ds = dataset('File',fullfile(sprintf(rootDir,subject(iSubject),dt,trialVarStr,architecture{iArchitecture}),sprintf(nameFVal,optimScope,modelNum)));
            
            % Load SAM
            load(fullfile(sprintf(rootDir,subject(iSubject),dt,trialVarStr,architecture{iArchitecture}),ds.FileNameSAM{1}),'SAM');
            
            
            % Extract optimized X
            iBestX = cell2mat(cellfun(@(in1) ~isempty(regexp(in1,'^BestX.*', 'once')),ds.Properties.VarNames,'Uni',0));
            X = double(ds(1,iBestX));
            
            
            
            
            % Print out otpimized X
            fprintf('Optimal parameters: \n')
            disp(X)
            
            prd = sam_sim_expt('explore',X,SAM);
            
            
            
            
            
            
            % ***************************************************************************
            % 2.2. Plot it
            % ***************************************************************************
            %             plotit(SAM,prd,p,model,iModel,modelStr,cost,altCost, responseSide, accuracy);
            
            
            
            % RT distribution Go trials: Correct Choices
            % ================================================
            if strcmp(accuracy, 'both') || strcmp(accuracy, 'correct')
                
                
                for iResponse = 1:length(responseArray)
                    for iCondition = 1:length(conditionArray);
                        
                        iRsp = responseArray(iResponse);
                        iCnd = conditionArray(iCondition);
                        
                        % Identify the relevent indices of parameters in X
                        xNames = SAM.model.variants.tree(modelNum).XSpec.name.go.name;
                        [tF, vInd] = ismember(sprintf('v_{GO,c%i}', iCnd),xNames);
                        [tF, z0Ind] = ismember(sprintf('z0_{GO,r%i}', iRsp),xNames);
                        switch modelNum
                            case 60
                                [tF, zcInd] = ismember(sprintf('zc_{GO}'),xNames);
                                [tF, t0Ind] = ismember(sprintf('t0_{GO}'),xNames);
                            case 137
                                [tF, zcInd] = ismember(sprintf('zc_{GO,r%i}', iRsp),xNames);
                                [tF, t0Ind] = ismember(sprintf('t0_{GO}'),xNames);
                            case 167
                                [tF, zcInd] = ismember(sprintf('zc_{GO}'),xNames);
                                [tF, t0Ind] = ismember(sprintf('t0_{GO,r%i}', iRsp),xNames);
                            otherwise
                                error(fprintf('Don''t have code yet for model %i\n', modelNum));
                        end
                        
                        
                        
                        % Identify the relevant rows in the dataset array
                        %         iTrialCatGo = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('goTrial.*GO:c%d',iCnd), 'once')),prd.trialCat,'Uni',0)));
                        iTrialCatGo = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('goTrial.*GO.*r%d.*c%d',iRsp,iCnd), 'once')),prd.trialCat,'Uni',0)));
                        if isempty(iTrialCatGo)
                            iTrialCatGo = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('goTrial_c%d.*',iCnd), 'once')),prd.trialCat,'Uni',0)));
                        end
                        
                        % GoCCorr trials
                        % -----------------------------------------------------------------
                        if strcmp(accuracy, 'both') || strcmp(accuracy, 'correct')
                            
                            % plot Average unit dynamics
                            dynTime = prd.dyn{iTrialCatGo}.goCCorr.goStim.targetGO.qX;
                            dynAct = prd.dyn{iTrialCatGo}.goCCorr.goStim.targetGO.qY;
                            plot(dynTime, dynAct, 'Color',modelLnClr(iModel,:),'LineStyle',modelLnStyle{iModel},'LineWidth',modelLnWidth(iModel))
                            
                            % plot starting point z0 and threshold zc
                            plot([0 rtLim],[X(z0Ind) X(z0Ind)], 'Color',modelLnClr(iModel,:),'LineStyle',modelLnStyle{iModel},'LineWidth',z0LnWidth)
                            plot([0 rtLim],[X(zcInd) X(zcInd)], 'Color',modelLnClr(iModel,:),'LineStyle',modelLnStyle{iModel},'LineWidth',zcLnWidth)
                        end
                    end
                    
                    
                end
                % Set axes
                switch subject(iSubject)
                    case 1
                        set(gca,'XLim',[0 rtLim], ...
                            'XTick',0:100:rtLim, ...
                            'YLim',[0 120], ...
                            'YTick',0:10:100)
                    case 2
                        set(gca,'XLim',[0 rtLim], ...
                            'XTick',0:100:rtLim, ...
                            'YLim',[0 120], ...
                            'YTick',0:10:100)
                    otherwise
                        error('Need to add axes limits for subject')
                end
            end
            
            
            % RT distribution Go trials: Error Choices
            % ================================================
            if  strcmp(accuracy, 'both') || strcmp(accuracy, 'error')
                
                p(3,iArchitecture).select();
                p(3,iArchitecture).hold('on');
                
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
                        set(gca,'XLim',[0 rtLim], ...
                            'XTick',0:100:rtLim, ...
                            'YLim',[0 100], ...
                            'YTick',0:10:100)
                    case 2
                        set(gca,'XLim',[0 rtLim], ...
                            'XTick',0:100:rtLim, ...
                            'YLim',[0 100], ...
                            'YTick',0:10:100)
                    otherwise
                        error('Need to add axes limits for subject')
                end
                
            end
            
            
            
            
            
            % Save plot
            if savePlot
                saveDir             = fullfile(sprintf(fileStr.result,subject(iSubject),dt,trialVarStr,'simulations'));
                if exist(saveDir,'dir') ~= 7
                    mkdir(saveDir)
                end
                fileName = sprintf('GO_Simulation_Respond_%s_Accuracy_%s', responseSide, accuracy);
                print(gcf, fullfile(saveDir, fileName),'-dpdf', '-r300')
                print(gcf, fullfile(saveDir, fileName),'-dpng')
            end
        end
    end
end



end


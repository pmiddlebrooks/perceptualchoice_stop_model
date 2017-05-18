function plot_simulation(subject,model,architecture,dt,trialVar,fileStr, savePlot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. PROCESS INPUTS & SPECIFY VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1.0 hard-code which conditions and responses to plot for now
nSim = 100;
close all
conditionArray       = 2;
responseArray   = 1;
accuracy        = 'both';
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
optimScope          = 'all';
rootDir             = fileStr.root;
nameFVal            = 'allFValAndBestX_%sTrials_model%.3d.txt';
modelStr            = {'zc','t0','v/ve','zc & t0','zc & v/ve','t0 & v/ve','zc, t0 & v/ve'};


% 1.4 Hard code some plot parameters
% =========================================================================
rtLim = 500;  % how far out do draw starting point and threshold
zcLim = 150;  % how far out do draw starting point and threshold
unitGoCorr = 1;
unitGoError = 2;
unitStop = 3;
unitLnStyle = {'-','-','-'};
unitLnWidth = [4 4 4];
unitLnClr = [.1 .7 .3; .6 .6 .6; ([255 50 19]/255)];
% unitLnSingleClr = [.3 .9 .5; .8 .8 .8; ([255 120 90]/255)];
unitLnSingleClr = [.2 .9 .3; .5 .9 .5; ([255 70 30]/255)];
z0LnWidth = 2;
zcLnWidth = 2;






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Plot the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop over subjects, architectures, and models
warning('off','stats:dataset:ModifiedVarnames')
for iSubject = 1:nSubject
    figureHandle = 19 + iSubject;
    
    switch subject(iSubject)
        case 1
            iSsdIndCancel = 3;
            iSsdIndNoncancel = 4;
        case 2
            iSsdIndCancel = 7;
    end
    %     switch subject(iSubject)
    %         case 1
    %             iSsdInd = 5;
    %             iSsdInd = 3;
    %         case 2
    %             iSsdInd = 16;
    %     end
    
    % Set up the figure and panels
    %         [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_figure;
    standard_figure(1,1,'landscape', figureHandle);
    
    %                 set_figure({1024,574,'pixels'},{'USLetter','landscape'},{'Helvetica',18});
    p = panel();
    p.pack({.04 .32 .32 .32}, 2);
    %     p.pack({.05 .475 .475}, num2cell(repmat(1/nModel,1,nModel)));
    
    
    for iArchitecture = 1:nArchitecture
        
        
        
        for iModel = 1:nModel
            modelNum = model(iModel);
            
            
            % Set up the figure and panels
            %         [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_figure;
            standard_figure(1,1,'landscape', figureHandle);
            
            %                 set_figure({1024,574,'pixels'},{'USLetter','landscape'},{'Helvetica',18});
            p = panel();
            p.pack({.04 .32 .32 .32}, 2);
            
            
            
            
            % Display progress
            fprintf(1,'Working on subject %d, architecture %s, and model %d. \n',subject(iSubject),architecture{iArchitecture},modelNum);
            
            % Load model-specific fits with cost function values
            ds = dataset('File',fullfile(sprintf(rootDir,subject(iSubject),dt,trialVarStr,architecture{iArchitecture}),sprintf(nameFVal,optimScope,modelNum)));
            
            % Load SAM
            load(fullfile(sprintf(rootDir,subject(iSubject),dt,trialVarStr,architecture{iArchitecture}),ds.FileNameSAM{1}),'SAM');
            
            % Specificy observations
            obs = SAM.optim.obs;
            
            % Extract optimized X
            iBestX = cell2mat(cellfun(@(in1) ~isempty(regexp(in1,'^BestX.*', 'once')),ds.Properties.VarNames,'Uni',0));
            X = double(ds(1,iBestX));
            
            
            
            
            % Print out otpimized X
            fprintf('Optimal parameters: \n')
            disp(X)
            
            % Alter the simulations for efficiency
            SAM.sim.n = nSim;
            
            prd = sam_sim_expt('explore',X,SAM);
            
            
            
            
            
            
            % ***************************************************************************
            % 2.2. Plot it
            % ***************************************************************************
            
            
            
            % RT distribution Go trials: Correct Choices
            % ================================================
            if strcmp(accuracy, 'both') || strcmp(accuracy, 'correct')
                
                
                
                iRspCorr = responseArray;
                iRspError = setxor(iRspCorr, [1 2]);
                iCnd = conditionArray;
                
                
                % Identify the relevent indices of parameters in X
                xNames = SAM.model.variants.tree(modelNum).XSpec.name.(optimScope).name;
                [tF, zcIndGOCorr] = ismember(sprintf('zc_{GO}'),xNames);
                zcIndGOError = zcIndGOCorr;
                [tF, t0Ind] = ismember(sprintf('t0_{GO}'),xNames);
                [tF, zcIndSTOP] = ismember(sprintf('zc_{STOP}'),xNames);
                [tF, z0IndSTOP] = ismember(sprintf('z0_{STOP}'),xNames);
                
                switch modelNum
                    case 79
                        [tF, z0IndGOCorr] = ismember(sprintf('z0_{GO,r%i}', iRspCorr),xNames);
                        [tF, z0IndGOError] = ismember(sprintf('z0_{GO,r%i}', iRspError),xNames);
                    case 352
                        [tF, z0IndGOCorr] = ismember(sprintf('z0_{GO}'),xNames);
                        z0IndGOError = z0IndGOCorr;
                    case 478
                        [tF, z0IndGOCorr] = ismember(sprintf('z0_{GO,r%i}', iRspCorr),xNames);
                        [tF, z0IndGOError] = ismember(sprintf('z0_{GO,r%i}', iRspError),xNames);
                    otherwise
                        error(fprintf('Don''t have code yet for model %i\n', modelNum));
                end
                
                
                
                
                
                
                
         
                
                
                
                
                
                
                % CANCELED STOP TRIAL
                
                
                 % Identify the relevant rows in the dataset array
                switch optimScope
                    case {'go'}
                        iTrialCatGo = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('goTrial.*GO.*r%d.*c%d',iRspCorr,iCnd), 'once')),prd.trialCat,'Uni',0)));
                        if isempty(iTrialCatGo)
                            iTrialCatGo = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('goTrial_c%d.*',iCnd), 'once')),prd.trialCat,'Uni',0)));
                        end
                    case {'all'}
                        iTrialCatGo = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('stopTrial.*ssd%d.*GO.*r%d.*c%d',iSsdIndCancel,iRspCorr,iCnd), 'once')),prd.trialCat,'Uni',0)));
                        
                        if isempty(iTrialCatGo)
                            iTrialCatGo = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('stopTrial.*ssd%d.*c%d.*r%d.*',iSsdIndCancel,iCnd,iRspCorr), 'once')),prd.trialCat,'Uni',0)));
                        end
                end
                
                iSsd = prd.ssd(iTrialCatGo);

                
                
                % GoCCorr unit
                % -----------------------------------------------------------------
                if strcmp(accuracy, 'both') || strcmp(accuracy, 'correct')
                    
                    p(2,1).select();
                    p(2,1).hold('on');
                    p(2,1).title({sprintf('CANCELED: %s',architecture{iArchitecture})});
                    set_the_axes
                    
                    % plot  unit dynamics
                    dynTimeGoCorr = prd.dyn{iTrialCatGo}.stopICorr.goStim.targetGO.sX;
                    dynActGoCorr = prd.dyn{iTrialCatGo}.stopICorr.goStim.targetGO.sY;
                    cellfun(@(x,y) plot(x,y, 'Color',unitLnSingleClr(unitGoCorr,:),'LineStyle',unitLnStyle{unitGoCorr},'LineWidth',1), dynTimeGoCorr,dynActGoCorr, 'uni', false)
                    
                    % plot starting point z0 and threshold zc
%                     plot([0 rtLim],[X(z0IndGOCorr) X(z0IndGOCorr)], 'Color',unitLnClr(unitGoCorr,:),'LineStyle',unitLnStyle{unitGoCorr},'LineWidth',z0LnWidth)
                    plot([0 rtLim],[X(zcIndGOCorr) X(zcIndGOCorr)], 'Color',unitLnClr(unitGoCorr,:),'LineStyle',unitLnStyle{unitGoCorr},'LineWidth',zcLnWidth)
                end
                
                
                
                % GoCError unit
                % -----------------------------------------------------------------
                if strcmp(accuracy, 'both')
                    
                    p(3,1).select();
                    p(3,1).hold('on');
                    set_the_axes
                    
                    % plot  unit dynamics
                    dynTimeGoError = prd.dyn{iTrialCatGo}.stopICorr.goStim.nonTargetGO.sX;
                    dynActGoError = prd.dyn{iTrialCatGo}.stopICorr.goStim.nonTargetGO.sY;
                    cellfun(@(x,y) plot(x,y, 'Color',unitLnSingleClr(unitGoError,:),'LineStyle',unitLnStyle{unitGoError},'LineWidth',1), dynTimeGoError,dynActGoError, 'uni', false)
                    
                    % plot starting point z0 and threshold zc
%                     plot([0 rtLim],[X(z0IndGOError) X(z0IndGOError)], 'Color',unitLnClr(unitGoError,:),'LineStyle',unitLnStyle{unitGoError},'LineWidth',z0LnWidth)
                    plot([0 rtLim],[X(zcIndGOError) X(zcIndGOError)], 'Color',unitLnClr(unitGoError,:),'LineStyle',unitLnStyle{unitGoError},'LineWidth',zcLnWidth)
                end
                
                
                
                % Stop unit
                % -----------------------------------------------------------------
                
                p(4,1).select();
                p(4,1).hold('on');
                set_the_axes
                
                % plot  unit dynamics
                dynTimeStop = prd.dyn{iTrialCatGo}.stopICorr.goStim.targetSTOP.sX;
                dynActStop = prd.dyn{iTrialCatGo}.stopICorr.goStim.targetSTOP.sY;
                cellfun(@(x,y) plot(x(iSsd+1:end),y(iSsd+1:end), 'Color',unitLnSingleClr(unitStop,:),'LineStyle',unitLnStyle{unitStop},'LineWidth',1), dynTimeStop,dynActStop, 'uni', false)
                dynTimeStop = nanmean(cell2mat(dynTimeStop));
                dynActStop = nanmean(cell2mat(dynActStop));
                
                % plot starting point z0 and threshold zc, and SSD
%                 plot([iSsd iSsd],[0 1.1*max([X(zcIndGOCorr), X(zcIndSTOP)])], 'Color',unitLnClr(unitStop,:),'LineStyle',unitLnStyle{unitStop},'LineWidth',z0LnWidth)
%                 plot([iSsd+1 rtLim],[X(z0IndSTOP) X(z0IndSTOP)], 'Color',unitLnClr(unitStop,:),'LineStyle',unitLnStyle{unitStop},'LineWidth',z0LnWidth)
                plot([iSsd+1 rtLim],[X(zcIndSTOP) X(zcIndSTOP)], 'Color',unitLnClr(unitStop,:),'LineStyle',unitLnStyle{unitStop},'LineWidth',zcLnWidth)
                
                
                
                % Plot mean funcitons
                %                 plot(nanmean(cell2mat(dynTimeGoCorr)), nanmean(cell2mat(dynActGoCorr)), 'Color',unitLnClr(unitGoCorr,:),'LineStyle',unitLnStyle{unitGoCorr},'LineWidth',unitLnWidth(unitGoCorr))
                %                 plot(nanmean(cell2mat(dynTimeGoError)), nanmean(cell2mat(dynActGoError)), 'Color',unitLnClr(unitGoError,:),'LineStyle',unitLnStyle{unitGoError},'LineWidth',unitLnWidth(unitGoError))
                %                 plot(dynTimeStop(iSsd+1:end), dynActStop(iSsd+1:end), 'Color',unitLnClr(unitStop,:),'LineStyle',unitLnStyle{unitStop},'LineWidth',unitLnWidth(unitStop))
                
                
                
                
                
                
                
     
                
                
                
                
                
                
                
                % NONCANCELED STOP TRIAL
                
                
                % Identify the relevant rows in the dataset array
                switch optimScope
                    case {'go'}
                        iTrialCatGo = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('goTrial.*GO.*r%d.*c%d',iRspCorr,iCnd), 'once')),prd.trialCat,'Uni',0)));
                        if isempty(iTrialCatGo)
                            iTrialCatGo = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('goTrial_c%d.*',iCnd), 'once')),prd.trialCat,'Uni',0)));
                        end
                    case {'all'}
                        iTrialCatGo = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('stopTrial.*ssd%d.*GO.*r%d.*c%d',iSsdIndNoncancel,iRspCorr,iCnd), 'once')),prd.trialCat,'Uni',0)));
                        
                        if isempty(iTrialCatGo)
                            iTrialCatGo = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('stopTrial.*ssd%d.*c%d.*r%d.*',iSsdIndNoncancel,iCnd,iRspCorr), 'once')),prd.trialCat,'Uni',0)));
                        end
                end
                
                iSsd = prd.ssd(iTrialCatGo);

                
                
                % GoCCorr unit
                % -----------------------------------------------------------------
                if strcmp(accuracy, 'both') || strcmp(accuracy, 'correct')
                    
                    p(2,2).select();
                    p(2,2).hold('on');
                    p(2,2).title({sprintf('NONCANCELED: %s',architecture{iArchitecture})});
                    set_the_axes
                    
                    % plot Average unit dynamics
                    dynTimeGoCorr = prd.dyn{iTrialCatGo}.stopIErrorCCorr.goStim.respGO.sX;
                    dynActGoCorr = prd.dyn{iTrialCatGo}.stopIErrorCCorr.goStim.respGO.sY;
                    cellfun(@(x,y) plot(x,y, 'Color',unitLnSingleClr(unitGoCorr,:),'LineStyle',unitLnStyle{unitGoCorr},'LineWidth',1), dynTimeGoCorr,dynActGoCorr, 'uni', false)
                    
                    % plot starting point z0 and threshold zc
%                     plot([0 rtLim],[X(z0IndGOCorr) X(z0IndGOCorr)], 'Color',unitLnClr(unitGoCorr,:),'LineStyle',unitLnStyle{unitGoCorr},'LineWidth',z0LnWidth)
                    plot([0 rtLim],[X(zcIndGOCorr) X(zcIndGOCorr)], 'Color',unitLnClr(unitGoCorr,:),'LineStyle',unitLnStyle{unitGoCorr},'LineWidth',zcLnWidth)
                end
                
                
                
                % GoCError unit
                % -----------------------------------------------------------------
                if strcmp(accuracy, 'both')
                    
                    p(3,2).select();
                    p(3,2).hold('on');
                    set_the_axes
                    
                    % plot Average unit dynamics
                    dynTimeGoError = prd.dyn{iTrialCatGo}.stopIErrorCCorr.goStim.nonTargetGO.sX;
                    dynActGoError = prd.dyn{iTrialCatGo}.stopIErrorCCorr.goStim.nonTargetGO.sY;
                    cellfun(@(x,y) plot(x,y, 'Color',unitLnSingleClr(unitGoError,:),'LineStyle',unitLnStyle{unitGoError},'LineWidth',1), dynTimeGoError,dynActGoError, 'uni', false)
                    
                    % plot starting point z0 and threshold zc
%                     plot([0 rtLim],[X(z0IndGOError) X(z0IndGOError)], 'Color',unitLnClr(unitGoError,:),'LineStyle',unitLnStyle{unitGoError},'LineWidth',z0LnWidth)
                    plot([0 rtLim],[X(zcIndGOError) X(zcIndGOError)], 'Color',unitLnClr(unitGoError,:),'LineStyle',unitLnStyle{unitGoError},'LineWidth',zcLnWidth)
                end
                
                
                % Stop unit
                % -----------------------------------------------------------------
                
                p(4,2).select();
                p(4,2).hold('on');
                set_the_axes
                
                % plot Average unit dynamics
                dynTimeStop = prd.dyn{iTrialCatGo}.stopIErrorCCorr.goStim.targetSTOP.sX;
                dynActStop = prd.dyn{iTrialCatGo}.stopIErrorCCorr.goStim.targetSTOP.sY;
                cellfun(@(x,y) plot(x(iSsd+1:end),y(iSsd+1:end), 'Color',unitLnSingleClr(unitStop,:),'LineStyle',unitLnStyle{unitStop},'LineWidth',1), dynTimeStop,dynActStop, 'uni', false)
                dynTimeStop = nanmean(cell2mat(dynTimeStop));
                dynActStop = nanmean(cell2mat(dynActStop));
                
                % plot starting point z0 and threshold zc
%                 plot([iSsd iSsd],[0 1.1*max([X(zcIndGOCorr), X(zcIndSTOP)])], 'Color',unitLnClr(unitStop,:),'LineStyle',unitLnStyle{unitStop},'LineWidth',z0LnWidth)
%                 plot([iSsd+1 rtLim],[X(z0IndSTOP) X(z0IndSTOP)], 'Color',unitLnClr(unitStop,:),'LineStyle',unitLnStyle{unitStop},'LineWidth',z0LnWidth)
                plot([iSsd+1 rtLim],[X(zcIndSTOP) X(zcIndSTOP)], 'Color',unitLnClr(unitStop,:),'LineStyle',unitLnStyle{unitStop},'LineWidth',zcLnWidth)
                
                
                
                % Plot mean funcitons
                % Normalize the
                % Get index of last non-nan value from each trial
                %                 lastInd = cellfun(@(x) find(isnan(x), 1)-1, dynActGoCorr);
                %
                %                 normFactor = repmat(mean(lastInd) ./ lastInd, 1, length(dynActGoCorr{1}));
                %
                %                 normFn = normFactor .* cell2mat(dynTimeGoCorr);
                %                 normFn = repmat(lastInd .* normFactor, (size(dynActGoCorr, 1)), length(dynActGoCorr{1}));
                %                 for i = 1 : size(dynActGoCorr, 1)
                %                     normFn(i, 1:lastInd(i)) = dynTimeGoCorr{i}(1:lastInd(i)) * normFactor(i);
                %                 end
                %                 plot(nanmean(normFn), nanmean(cell2mat(dynActGoCorr)), 'Color',unitLnClr(unitGoCorr,:),'LineStyle',unitLnStyle{unitGoCorr},'LineWidth',unitLnWidth(unitGoCorr));
                
                
                %                 plot(nanmean(cell2mat(dynTimeGoCorr)), nanmean(cell2mat(dynActGoCorr)), 'Color',unitLnClr(unitGoCorr,:),'LineStyle',unitLnStyle{unitGoCorr},'LineWidth',unitLnWidth(unitGoCorr))
                %                 plot(nanmean(cell2mat(dynTimeGoError)), nanmean(cell2mat(dynActGoError)), 'Color',unitLnClr(unitGoError,:),'LineStyle',unitLnStyle{unitGoError},'LineWidth',unitLnWidth(unitGoError))
                %                 plot(dynTimeStop(iSsd+1:end), dynActStop(iSsd+1:end), 'Color',unitLnClr(unitStop,:),'LineStyle',unitLnStyle{unitStop},'LineWidth',unitLnWidth(unitStop))
                
                
                
                
                %                 % Set axes
                %                 switch subject(iSubject)
                %                     case 1
                %                         set(gca,'XLim',[0 rtLim], ...
                %                             'XTick',0:100:rtLim, ...
                %                             'YLim',[0 1.1*max([X(zcIndGOCorr), X(zcIndSTOP)])], ...
                %                             'YTick',0:10:100)
                %                     case 2
                %                         set(gca,'XLim',[0 rtLim], ...
                %                             'XTick',0:100:rtLim, ...
                %                             'YLim',[0 1.1*max([X(zcIndGOCorr), X(zcIndSTOP)])], ...
                %                             'YTick',0:10:100)
                %                     otherwise
                %                         error('Need to add axes limits for subject')
                %                 end
                %             end
                %
                
                %             % RT distribution Go trials: Error Choices
                %             % ================================================
                %             if  strcmp(accuracy, 'both') || strcmp(accuracy, 'error')
                %
                %                 p(3,iArchitecture).select();
                %                 p(3,iArchitecture).hold('on');
                %
                %                 for iResponse = 1:length(responseArray)
                %                     for iCondition = 1:length(conditionArray);
                %
                %                         iRsp = responseArray(iResponse);
                %                         iCnd = conditionArray(iCondition);
                %
                %                         % Identify the relevant rows in the dataset array
                %                         %         iTrialCatGo = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('goTrial.*GO:c%d',iCnd), 'once')),prd.trialCat,'Uni',0)));
                %                         iTrialCatGo = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('goTrial.*GO.*r%d.*c%d',iRsp,iCnd), 'once')),prd.trialCat,'Uni',0)));
                %                         if isempty(iTrialCatGo)
                %                             iTrialCatGo = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('goTrial_c%d.*',iCnd), 'once')),prd.trialCat,'Uni',0)));
                %                         end
                %
                %
                %                         % GoCError trials
                %                         % -----------------------------------------------------------------
                %
                %                             rtGoCErrorObs     = obs.rtQGoCError{iTrialCatGo};
                %                             cumPGoCErrorObs   = obs.cumProbGoCError{iTrialCatGo};
                %
                %                             rtGoCErrorPrd     = prd.rtGoCError{iTrialCatGo};
                %                             cumPGoCErrorPrd   = cmtb_edf(prd.rtGoCError{iTrialCatGo}(:),prd.rtGoCError{iTrialCatGo}(:));
                %
                %                             % Plot it
                %                             iColor = cndClr((length(conditionArray) * (iRsp -1) + iCondition), :);
                %                             plot(rtGoCErrorObs,cumPGoCErrorObs,'Color',iColor,'Marker',goCErrorMrkObs,'LineStyle',goCErrorLnObs,'LineWidth',goCErrorLnWidth);
                %                             plot(rtGoCErrorPrd,cumPGoCErrorPrd,'Color',iColor,'Marker',goCErrorMrkPrd,'LineStyle',goCErrorLnPrd,'LineWidth',goCErrorLnWidth);
                %
                %                     end
                %                 end
                %
                %                 % Set axes
                %                 switch subject(iSubject)
                %                     case 1
                %                         set(gca,'XLim',[0 rtLim], ...
                %                             'XTick',0:100:rtLim, ...
                %                             'YLim',[0 120], ...
                %                             'YTick',0:10:100)
                %                     case 2
                %                         set(gca,'XLim',[0 rtLim], ...
                %                             'XTick',0:100:rtLim, ...
                %                             'YLim',[0 120], ...
                %                             'YTick',0:10:100)
                %                     otherwise
                %                         error('Need to add axes limits for subject')
                %                 end
                %
                %             end
                
                
                
                
                set(gcf,'renderer','painters')
                % Save plot
                if savePlot
                    saveDir             = fullfile(sprintf(fileStr.result,subject(iSubject),dt,trialVarStr,'simulations'));
                    if exist(saveDir,'dir') ~= 7
                        mkdir(saveDir)
                    end
                    fileName = sprintf('%s_%s_Simulation_Respond_%s_Accuracy_%s_SSD_%d', architecture{1}, optimScope, responseSide, accuracy, iSsd);
                    print(gcf, fullfile(saveDir, fileName),'-dpdf', '-r300')
                    print(gcf, fullfile(saveDir, fileName),'-depsc', '-r300')
                end
            end
        end
    end
end



    function set_the_axes
        % Set axes
        switch subject(iSubject)
            case 1
                set(gca,'XLim',[0 rtLim], ...
                    'XTick',0:100:rtLim, ...
                    'YLim',[0 1.1*max([X(zcIndGOCorr), X(zcIndSTOP)])], ...
                    'YTick',0:10:100)
            case 2
                set(gca,'XLim',[0 rtLim], ...
                    'XTick',0:100:rtLim, ...
                    'YLim',[0 1.1*max([X(zcIndGOCorr), X(zcIndSTOP)])], ...
                    'YTick',0:10:100)
            otherwise
                error('Need to add axes limits for subject')
        end
    end

end






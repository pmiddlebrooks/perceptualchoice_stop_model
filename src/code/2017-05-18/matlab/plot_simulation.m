function plot_simulation(subject,model,architecture,dt,trialVar,fileStr, savePlot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. PROCESS INPUTS & SPECIFY VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1.0 hard-code which conditions and responses to plot for now
close all
conditionArray       = [1,2]; % 1 is easy, 2 is hard
accuracy        = 'both';
responseSide    = {'left','right'};
nSim = 100;

newOrRepeatSimulation = 'new';
rngSeedId = 0;

DO_STOPS = false;

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
rtLim = 700;  % how far out do draw starting point and threshold
zcLim = 150;  % how far out do draw starting point and threshold
unitGoCorr = 1;
unitGoError = 2;
unitStop = 3;
unitLnStyle = {'-','-','-'};
unitLnWidth = [4 4 4];
unitLnClr = [.1 .7 .3; .6 .6 .6; ([255 50 19]/255)];
unitLnSingleClr = [.3 .9 .5; .5 .5 .5; ([255 120 90]/255)];
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
            iSsdIndCancel = 4;
            iSsdIndNoncancel = 4;
        case 3
            iSsdIndCancel = 3;
            iSsdIndNoncancel = 3;
    end
    
    
    % Set up the figure and panels
    %         [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_figure;
    standard_figure(1,1,'landscape', figureHandle);
    
    %                 set_figure({1024,574,'pixels'},{'USLetter','landscape'},{'Helvetica',18});
    %     p.pack({.05 .475 .475}, num2cell(repmat(1/nModel,1,nModel)));
    
    
    for iArchitecture = 1:nArchitecture
        
        for iModel = 1:nModel
            modelNum = model(iModel);
            
            
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
            
            switch newOrRepeatSimulation
                case 'new'
                    SAM.sim.rng.id = randi(100);  % set a random rng for a new simulation
                case 'repeat'
                    SAM.sim.rng.id = rngSeedId;
            end
            
            prd = sam_sim_expt('explore',X,SAM);
            
            
            
            
            
            
            % Loop through resonse side and conditions
            
            for iRspInd = 1:length(responseSide)
                
                if strcmp(responseSide{iRspInd}, 'right')
                    responseArray   = 2;
                elseif strcmp(responseSide{iRspInd}, 'left')
                    responseArray   = 1;
                end
                
                for iCndInd = 1 : length(conditionArray)
                    

                    clf
                        p = panel();
    p.pack({.04 .32 .32 .32}, 2);

                    % ***************************************************************************
                    % 2.2. Plot it
                    % ***************************************************************************
                    
                    
                    
                    % RT distribution Go trials: Correct Choices
                    % ================================================
                    
                    
                    
                    iRspCorr = responseArray;
                    iRspError = setxor(iRspCorr, [1 2]);
                    iCnd = conditionArray(iCndInd);
                    
                    
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
                    
                    
                    
                    
                    % CORRECT CHOICE GO TRIAL
                    
                    
                    % Identify the relevant rows in the dataset array
                    iTrialCatGo = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('goTrial.*GO.*r%d.*c%d',iRspCorr,iCnd), 'once')),prd.trialCat,'Uni',0)));
                    if isempty(iTrialCatGo)
                        iTrialCatGo = find(cell2mat(cellfun(@(in1) ~isempty(regexp(in1,sprintf('goTrial_c%d.*',iCnd), 'once')),prd.trialCat,'Uni',0)));
                    end
                    goTrialCat = iTrialCatGo;
                    
                    
                    % GoCCorr unit
                    % -----------------------------------------------------------------
                    if strcmp(accuracy, 'both') || strcmp(accuracy, 'correct')
                        
                        p(4,1).select();
                        p(4,1).hold('on');
                        p(4,1).title({sprintf('GO CORRECT: %s Mean RT: %d',architecture{iArchitecture}, round(mean(prd.rtGoCCorr{iTrialCatGo})))});
                        set_the_axes('go')
                        
                        % plot  unit dynamics
                        dynTimeGoCorr = prd.dyn{iTrialCatGo}.goCCorr.goStim.targetGO.sX;
                        dynActGoCorr = prd.dyn{iTrialCatGo}.goCCorr.goStim.targetGO.sY;
                        cellfun(@(x,y) plot(x,y, 'Color',unitLnSingleClr(unitGoCorr,:),'LineStyle',unitLnStyle{unitGoCorr},'LineWidth',1), dynTimeGoCorr,dynActGoCorr, 'uni', false)
                        
                        % plot starting point z0 and threshold zc
                        %                     plot([0 rtLim],[X(z0IndGOCorr) X(z0IndGOCorr)], 'Color',unitLnClr(unitGoCorr,:),'LineStyle',unitLnStyle{unitGoCorr},'LineWidth',z0LnWidth)
                        plot([0 rtLim],[X(zcIndGOCorr) X(zcIndGOCorr)], 'Color',unitLnClr(unitGoCorr,:),'LineStyle',unitLnStyle{unitGoCorr},'LineWidth',zcLnWidth)
                        
                        fprintf('\nMean RT GO Correct: %d\n', round(mean(prd.rtGoCCorr{iTrialCatGo})))
                        
                    end
                    
                    % GoCError unit
                    % -----------------------------------------------------------------
                    if strcmp(accuracy, 'both')
                        
                        %                     % plot  unit dynamics
                        dynTimeGoError = prd.dyn{iTrialCatGo}.goCCorr.goStim.nonTargetGO.sX;
                        dynActGoError = prd.dyn{iTrialCatGo}.goCCorr.goStim.nonTargetGO.sY;
                        cellfun(@(x,y) plot(x,y, 'Color',unitLnSingleClr(unitGoError,:),'LineStyle',unitLnStyle{unitGoError},'LineWidth',1), dynTimeGoError,dynActGoError, 'uni', false)
                    end
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    % ERROR CHOICE GO TRIAL
                    
                    p(4,2).select();
                    p(4,2).hold('on');
                    p(4,2).title({sprintf('GO Error: %s Mean RT: %d',architecture{iArchitecture}, round(mean(prd.rtGoCError{iTrialCatGo})))});
                    set_the_axes('go')
 if isfield(prd.dyn{iTrialCatGo}, 'goCError')                   
                    % GoCCorr unit
                    % -----------------------------------------------------------------
                    if strcmp(accuracy, 'both') || strcmp(accuracy, 'correct')
                        
                        % plot  unit dynamics
                        dynTimeGoCorr = prd.dyn{iTrialCatGo}.goCError.goStim.targetGO.sX;
                        dynActGoCorr = prd.dyn{iTrialCatGo}.goCError.goStim.targetGO.sY;
                        cellfun(@(x,y) plot(x,y, 'Color',unitLnSingleClr(unitGoCorr,:),'LineStyle',unitLnStyle{unitGoCorr},'LineWidth',1), dynTimeGoCorr,dynActGoCorr, 'uni', false)
                        
                        % plot starting point z0 and threshold zc
                        %                     plot([0 rtLim],[X(z0IndGOCorr) X(z0IndGOCorr)], 'Color',unitLnClr(unitGoCorr,:),'LineStyle',unitLnStyle{unitGoCorr},'LineWidth',z0LnWidth)
                        plot([0 rtLim],[X(zcIndGOCorr) X(zcIndGOCorr)], 'Color',unitLnClr(unitGoCorr,:),'LineStyle',unitLnStyle{unitGoCorr},'LineWidth',zcLnWidth)
                        
                        
                    end
                    
                    % GoCError unit
                    % -----------------------------------------------------------------
                    if strcmp(accuracy, 'both')
                        
                        %                     % plot  unit dynamics
                        dynTimeGoError = prd.dyn{iTrialCatGo}.goCError.goStim.nonTargetGOError.sX;
                        dynActGoError = prd.dyn{iTrialCatGo}.goCError.goStim.nonTargetGOError.sY;
                        cellfun(@(x,y) plot(x,y, 'Color',unitLnSingleClr(unitGoError,:),'LineStyle',unitLnStyle{unitGoError},'LineWidth',1), dynTimeGoError,dynActGoError, 'uni', false)
                        
                        fprintf('\nMean RT GO Error: %d\n', round(mean(prd.rtGoCError{iTrialCatGo})))
                    end
                    
 end                   
                    
                    fprintf('Respond_%s_Accuracy_%s_Coherence_%d:\n', responseSide{iRspInd}, accuracy, iCnd)
                    prd.modelMat{iTrialCatGo}.Z0
                    prd.modelMat{iTrialCatGo}.V
                    
                    
                    
                    
                    
                    
                    
                        
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
                        
                        
                     if DO_STOPS
                       
                        % GoCCorr unit
                        % -----------------------------------------------------------------
                        if strcmp(accuracy, 'both') || strcmp(accuracy, 'correct')
                            
                            p(2,1).select();
                            p(2,1).hold('on');
                            p(2,1).title({sprintf('CANCELED: %s',architecture{iArchitecture})});
                            set_the_axes('go')
                            
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
                            
                            %                     p(3,1).select();
                            %                     p(3,1).hold('on');
                            %                     set_the_axes('go')
                            %
                            %                     % plot  unit dynamics
                            %                     dynTimeGoError = prd.dyn{iTrialCatGo}.stopICorr.goStim.nonTargetGO.sX;
                            %                     dynActGoError = prd.dyn{iTrialCatGo}.stopICorr.goStim.nonTargetGO.sY;
                            %                     cellfun(@(x,y) plot(x,y, 'Color',unitLnSingleClr(unitGoError,:),'LineStyle',unitLnStyle{unitGoError},'LineWidth',1), dynTimeGoError,dynActGoError, 'uni', false)
                            %
                            %                     % plot starting point z0 and threshold zc
                            %                     %                     plot([0 rtLim],[X(z0IndGOError) X(z0IndGOError)], 'Color',unitLnClr(unitGoError,:),'LineStyle',unitLnStyle{unitGoError},'LineWidth',z0LnWidth)
                            %                     plot([0 rtLim],[X(zcIndGOError) X(zcIndGOError)], 'Color',unitLnClr(unitGoError,:),'LineStyle',unitLnStyle{unitGoError},'LineWidth',zcLnWidth)
                        end
                        
                        
                        
                        % Stop unit
                        % -----------------------------------------------------------------
                        
                        p(3,1).select();
                        p(3,1).hold('on');
                        set_the_axes('stop')
                        
                        % plot  unit dynamics
                        dynTimeStop = prd.dyn{iTrialCatGo}.stopICorr.goStim.targetSTOP.sX;
                        dynActStop = prd.dyn{iTrialCatGo}.stopICorr.goStim.targetSTOP.sY;
                        cellfun(@(x,y) plot(x(iSsd+1:end),y(iSsd+1:end), 'Color',unitLnSingleClr(unitStop,:),'LineStyle',unitLnStyle{unitStop},'LineWidth',1), dynTimeStop,dynActStop, 'uni', false)
                        dynTimeStop = nanmean(cell2mat(dynTimeStop));
                        dynActStop = nanmean(cell2mat(dynActStop));
                        
                        % plot starting point z0 and threshold zc, and SSD
                        %                     plot([iSsd iSsd],[0 1.1*max([X(zcIndGOCorr), X(zcIndSTOP)])], 'Color',unitLnClr(unitStop,:),'LineStyle',unitLnStyle{unitStop},'LineWidth',z0LnWidth)
                        %                     plot([iSsd+1 rtLim],[X(z0IndSTOP) X(z0IndSTOP)], 'Color',unitLnClr(unitStop,:),'LineStyle',unitLnStyle{unitStop},'LineWidth',z0LnWidth)
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
                            set_the_axes('go')
                            
                            
                            % plot Average unit dynamics
                            dynTimeGoCorr = prd.dyn{iTrialCatGo}.stopIErrorCCorr.goStim.respGO.sX;
                            dynActGoCorr = prd.dyn{iTrialCatGo}.stopIErrorCCorr.goStim.respGO.sY;
                            cellfun(@(x,y) plot(x,y, 'Color',unitLnSingleClr(unitGoCorr,:),'LineStyle',unitLnStyle{unitGoCorr},'LineWidth',1), dynTimeGoCorr,dynActGoCorr, 'uni', false)
                            
                            fprintf('\nMean RT Noncanceled: %d\n', round(mean(prd.rtStopIErrorCCorr{iTrialCatGo})))
                            % plot starting point z0 and threshold zc
                            %                     plot([0 rtLim],[X(z0IndGOCorr) X(z0IndGOCorr)], 'Color',unitLnClr(unitGoCorr,:),'LineStyle',unitLnStyle{unitGoCorr},'LineWidth',z0LnWidth)
                            plot([0 rtLim],[X(zcIndGOCorr) X(zcIndGOCorr)], 'Color',unitLnClr(unitGoCorr,:),'LineStyle',unitLnStyle{unitGoCorr},'LineWidth',zcLnWidth)
                        end
                        
                        
                        
                        % GoCError unit
                        % -----------------------------------------------------------------
                        if strcmp(accuracy, 'both')
                            
                            %                     p(3,2).select();
                            %                     p(3,2).hold('on');
                            %                     set_the_axes('go')
                            %
                            %                     % plot Average unit dynamics
                            %                     dynTimeGoError = prd.dyn{iTrialCatGo}.stopIErrorCCorr.goStim.nonTargetGO.sX;
                            %                     dynActGoError = prd.dyn{iTrialCatGo}.stopIErrorCCorr.goStim.nonTargetGO.sY;
                            %                     cellfun(@(x,y) plot(x,y, 'Color',unitLnSingleClr(unitGoError,:),'LineStyle',unitLnStyle{unitGoError},'LineWidth',1), dynTimeGoError,dynActGoError, 'uni', false)
                            %
                            %                     % plot starting point z0 and threshold zc
                            %                     %                     plot([0 rtLim],[X(z0IndGOError) X(z0IndGOError)], 'Color',unitLnClr(unitGoError,:),'LineStyle',unitLnStyle{unitGoError},'LineWidth',z0LnWidth)
                            %                     plot([0 rtLim],[X(zcIndGOError) X(zcIndGOError)], 'Color',unitLnClr(unitGoError,:),'LineStyle',unitLnStyle{unitGoError},'LineWidth',zcLnWidth)
                        end
                        
                        
                        
                        % Stop unit
                        % -----------------------------------------------------------------
                        
                        p(3,2).select();
                        p(3,2).hold('on');
                        set_the_axes('stop')
                        
                        % plot Average unit dynamics
                        dynTimeStop = prd.dyn{iTrialCatGo}.stopIErrorCCorr.goStim.targetSTOP.sX;
                        dynActStop = prd.dyn{iTrialCatGo}.stopIErrorCCorr.goStim.targetSTOP.sY;
                        cellfun(@(x,y) plot(x(iSsd+1:end),y(iSsd+1:end), 'Color',unitLnSingleClr(unitStop,:),'LineStyle',unitLnStyle{unitStop},'LineWidth',1), dynTimeStop,dynActStop, 'uni', false)
                        dynTimeStop = nanmean(cell2mat(dynTimeStop));
                        dynActStop = nanmean(cell2mat(dynActStop));
                        
                        % plot starting point z0 and threshold zc
                        %                     plot([iSsd iSsd],[0 1.1*max([X(zcIndGOCorr), X(zcIndSTOP)])], 'Color',unitLnClr(unitStop,:),'LineStyle',unitLnStyle{unitStop},'LineWidth',z0LnWidth)
                        %                     plot([iSsd+1 rtLim],[X(z0IndSTOP) X(z0IndSTOP)], 'Color',unitLnClr(unitStop,:),'LineStyle',unitLnStyle{unitStop},'LineWidth',z0LnWidth)
                        plot([iSsd+1 rtLim],[X(zcIndSTOP) X(zcIndSTOP)], 'Color',unitLnClr(unitStop,:),'LineStyle',unitLnStyle{unitStop},'LineWidth',zcLnWidth)
                        
                        
                        
                    end
                    
                    
                    
                    set(gcf,'renderer','painters')
                    % Save plot
                    if savePlot
                        saveDir             = fullfile(sprintf(fileStr.result,subject(iSubject),dt,trialVarStr,'simulations'));
                        if exist(saveDir,'dir') ~= 7
                            mkdir(saveDir)
                        end
                        fileName = sprintf('%s_%s_Simulation_Respond_%s_Accuracy_%s_Coherence_%d_SSD_%d_RT_%d', architecture{1}, optimScope, responseSide{iRspInd}, accuracy, iCnd, iSsd, round(mean(prd.rtGoCCorr{goTrialCat})));
                        print(gcf, fullfile(saveDir, fileName),'-dpdf', '-r300')
                        %                     print(gcf, fullfile(saveDir, fileName),'-depsc', '-r300')
                    end
                    
                end
            end
        end
    end
end


    function set_the_axes(goOrStopUnit)
        % Set axes
        switch goOrStopUnit
            case 'go'
                switch subject(iSubject)
                    case 1
                        set(gca,'XLim',[0 rtLim], ...
                            'XTick',0:100:rtLim, ...
                            'YLim',[0 max([X(zcIndGOCorr), X(zcIndSTOP)])], ...  % 'YLim',[0 1.1*max([X(zcIndGOCorr), X(zcIndSTOP)])]
                            'YTick',0:40:160)
                    case 3
                        set(gca,'XLim',[0 rtLim], ...
                            'XTick',0:100:rtLim, ...
                            'YLim',[0 1.1*max([X(zcIndGOCorr), X(zcIndSTOP)])], ...
                            'YTick',0:10:100)
                    otherwise
                        error('Need to add axes limits for subject')
                end
            case 'stop'
                switch subject(iSubject)
                    case 1
                        set(gca,'XLim',[0 rtLim], ...
                            'XTick',0:100:rtLim, ...
                            'YLim',[0 100], ...  % 'YLim',[0 1.1*max(X(zcIndSTOP))], ...
                            'YTick',0:10:50)
                    case 3
                        set(gca,'XLim',[0 rtLim], ...
                            'XTick',0:100:rtLim, ...
                            'YLim',[0 1.1*max(X(zcIndSTOP))], ...
                            'YTick',0:10:30)
                    otherwise
                        error('Need to add axes limits for subject')
                end
        end
        
    end
end


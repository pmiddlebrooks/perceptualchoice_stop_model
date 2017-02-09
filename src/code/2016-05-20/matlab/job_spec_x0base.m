function job_spec_x0base(subj,trialVar,optimScope,choiceMech,stopMech,iModel,fileStr,x0Base,doPlot,doSave,doStartParCluster)
%
% INPUTS
% subj              vector of subject indices
% trialVar          true or false, influences t0, z0, and eta
% optimScope        'go','stop', or 'all'
% choiceMech        'race', 'ffi', or 'li'
% stopMech          'race', 'bi', or 'li'
% iModel            scalar of the model for which to specify userSpec*.txt
% fileStr           struct containing formated strings linking to
% *rootDir          - root directory
% * bestXGO         - path to best-fitting GO parameters (optimScope == 'stop')
% * bestXSTOP       -
% x0Base            cell array containing starting values. cell indices should correspond to subject indices
% doPlot            true or false, whether to plot observed and predicted distributions or not
% doSave            true or false, whether to save parameters
% doStartParCluster true or false, whether to start the parallel cluster
%                   - set to false, when
%                     * not plotting observed and predicted distributions
%                     * parallel cluster is already started
%
% SYNTAX
% job_spec_x0base(subj,trialVar,optimScope,choiceMech,stopMech,iModel,fileStr,x0Base,doPlot,doSave,doStartParCluster)
%

if trialVar
    trialVarStr = 'trialvar';
else
    trialVarStr = 'notrialvar';
end

switch optimScope
    case 'go'
        modelArch = sprintf('c%s',choiceMech);
        stopMech = 'none';
    case 'stop'
        modelArch = sprintf('c%s_i%s',choiceMech,stopMech);
        modelArchGO = sprintf('c%s',choiceMech);
    case 'all'
        modelArch = sprintf('c%s_i%s',choiceMech,stopMech);
end

rootDir           = fileStr.root;
switch optimScope
    case 'stop'
        bestXGO       = fileStr.bestXGO;
    case 'all'
        bestXSTOP     = fileStr.bestXSTOP;
end
nameSAMGeneral    = 'SAM_%sTrials.mat';
nameSAMModel      = 'SAM_%sTrials_model%.3d.mat';
nameUserSpecX     = 'userSpecX_%sTrials_model%.3d.txt';

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% #. START THE PARALLEL POOL
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if doStartParCluster
    
    % Load SAM file for first subject, to get access to compCluster field
    load(fullfile(sprintf(rootDir,subj(1),trialVarStr,modelArch),sprintf(nameSAMGeneral,optimScope)));
    
    % In unique directory to prevent collision of parallel jobs
    % e.g. see: http://www.mathworks.com/matlabcentral/answers/97141-why-am-i-unable-to-start-a-local-matlabpool-from-multiple-matlab-sessions-that-use-a-shared-preferen
    c = parcluster();
    if isfield(SAM,'compCluster')
        c.NumWorkers = SAM.compCluster.nProcessors;
    else
        c.NumWorkers = 1;
    end
    [~,homeDir] = system('echo $HOME');
    homeDir = strtrim(homeDir);
    release = version('-release')
    tempDir = fullfile(homeDir,'.matlab','local_cluster_jobs',release);
    if exist(tempDir) ~= 7
        mkdir(tempDir)
    end
    t = tempname(tempDir);
    mkdir(t);
    c.JobStorageLocation=t;
    tWait = 1+60*rand();
    pause(tWait);
    myPool = parpool(c);
    
    clear SAM
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FINETUNING BY HAND
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for iSubjInd = 1 : length(subj)
    iSubj = subj(iSubjInd);
    
    % Load the general SAM file
    % =========================================================================
    load(fullfile(sprintf(rootDir,iSubj,trialVarStr,modelArch),sprintf(nameSAMGeneral,optimScope)));
    
    % If fitting all trials, load best-fitting GO parameters and add STOP parameters
    % =========================================================================
    switch optimScope
        case 'go'
            
            % Determine number of parameter categories
            nXCat = SAM.model.XCat.n;
            
            % Pre-allocate X
            X = nan(1,sum(SAM.model.variants.tree(iModel).XSpec.n.nCatClass(1,:)));
            
            for iXCat = 1:nXCat
                
                if ~SAM.model.XCat.included(iXCat)
                    % Not included => use valExcluded
                    % =======================================================
                    
                    % Identify the column indices for this paramater in the all model
                    iColX = SAM.model.variants.tree(iModel).XSpec.i.go.iCatClass{1,iXCat};
                    
                    % Extract the default parameter value and put in X
                    X(iColX) = SAM.model.XCat.valExcluded(iXCat);
                    
                else
                    
                    % Identify the column indices for this paramater in the all model
                    iColX = SAM.model.variants.tree(iModel).XSpec.i.go.iCatClass{1,iXCat};
                    
                    X(iColX) = x0Base{iSubj}{iXCat};
                    
                end
            end
            
        case 'stop'
            
            % Load best-fitting GO parameters
            ds = dataset('File',fullfile(sprintf(rootDir,iSubj,trialVarStr,modelArchGO),sprintf(bestXGO,iModel)),'HeaderLines',0);
            
            % Determine number of parameter categories
            nXCat = SAM.model.XCat.n;
            
            % Pre-allocate X
            X = nan(1,SAM.model.variants.tree(iModel).XSpec.n.n);
            
            for iXCat = 1:nXCat
                
                if ~SAM.model.XCat.included(iXCat)
                    % Not included => use valExcluded
                    % =======================================================
                    
                    % Identify the column indices for this paramater in the stop model
                    iColX = SAM.model.variants.tree(iModel).XSpec.i.stop.iCatClass{1,iXCat};
                    
                    % Extract the default parameter value and put in X
                    X(iColX) = SAM.model.XCat.valExcluded(iXCat);
                    
                elseif SAM.model.XCat.included(iXCat) && ~SAM.model.XCat.classSpecific(iXCat)
                    % Included, not class-specific
                    % - GO parameters: use best-fitting GO values
                    % - STOP parameters: use best-fitting GO values
                    % =======================================================
                    
                    % Identify the column indices for this paramater in the go model
                    iColBestXGO = SAM.model.variants.tree(iModel).XSpec.i.go.iCatClass{iXCat};
                    
                    % Identify the column indices for this paramater in the stop model
                    iColXGO = SAM.model.variants.tree(iModel).XSpec.i.stop.iCatClass{1,iXCat};
                    iColXSTOP = SAM.model.variants.tree(iModel).XSpec.i.stop.iCatClass{2,iXCat};
                    
                    % Construct the column names based on the go model indices
                    colNames    = arrayfun(@(in1) sprintf('BestX_%d',in1),iColBestXGO,'Uni',0);
                    
                    % Extract the best-fitting GO parameters and put in X
                    X(iColXGO) = cell2mat(arrayfun(@(in1) double(ds(1,in1)),colNames,'Uni',0));
                    X(iColXSTOP) = cell2mat(arrayfun(@(in1) double(ds(1,in1)),colNames,'Uni',0));
                    
                elseif SAM.model.XCat.included(iXCat) && SAM.model.XCat.classSpecific(iXCat)
                    % Included, class-specific
                    % - GO parameters: use best-fitting GO values
                    % - STOP parameters: use user-specified values
                    % =======================================================
                    
                    % Identify the column indices for this paramater in the go model
                    iColBestXGO = SAM.model.variants.tree(iModel).XSpec.i.go.iCatClass{iXCat};
                    
                    % Identify the column indices for this paramater in the stop model
                    iColXGO = SAM.model.variants.tree(iModel).XSpec.i.stop.iCatClass{1,iXCat};
                    iColXSTOP = SAM.model.variants.tree(iModel).XSpec.i.stop.iCatClass{2,iXCat};
                    
                    % Construct the column names based on the go model indices
                    colNames    = arrayfun(@(in1) sprintf('BestX_%d',in1),iColBestXGO,'Uni',0);
                    
                    % Extract the best-fitting GO parameters and put in X
                    X(iColXGO) = cell2mat(arrayfun(@(in1) double(ds(1,in1)),colNames,'Uni',0));
                    
                    if ~isempty(iColXSTOP)
                        % Parameter ve is not included if there is just one STOP unit
                        X(iColXSTOP) = x0Base{iSubj}{iXCat};
                    end
                end
            end
            
        case 'all'
            
            % Load best-fitting GO and STOP parameters
            ds = dataset('File',fullfile(sprintf(bestXSTOP,iSubj,trialVarStr,modelArch,iModel)),'HeaderLines',0);
            
            % Extract optimized X
            iBestX = cell2mat(cellfun(@(in1) ~isempty(regexp(in1,'^BestX.*', 'once')),ds.Properties.VarNames,'Uni',0));
            X = double(ds(1,iBestX));
            
    end
    
    % Save X
    % =========================================================================
    if doSave
        save(fullfile(sprintf(rootDir,iSubj,trialVarStr,modelArch),sprintf(nameUserSpecX,optimScope,iModel)),'-ascii','-double','-tabs','X');
    end
    
    if doPlot
        
        % Specify model-specific SAM structure
        % =========================================================================
        SAM                   = sam_spec_job_specific(SAM,iModel);
        
        SAM.sim.n = 400;
        % Predict RTs
        % =========================================================================
        [cost,altCost,prd] = sam_cost(X,SAM);
        
        % Plot observations and predictions
        % =========================================================================
        
        switch optimScope
            case 'go'
              % Do nothing. Only need to downsample stop trial categories for plotting  
            case {'stop'}
                % Plot the trial categories with the top nPlot number of trials
                nPlot = 10;
                [s ind] = sort(SAM.optim.obs.nStopIErrorCCorr);
%                 [s ind] = sort(SAM.optim.obs.nTotal);
                indPlot = ind(end - (nPlot-1) : end);
                prd = prd(indPlot,:);
                SAM.optim.obs = SAM.optim.obs(indPlot,:);
        end  
        sam_plot(SAM,prd);
        
        
        fprintf(1,'BIC = %.2f, Chi^2 = %.2f \n',sum(cell2mat(prd.bic)),sum(cell2mat(prd.chiSquare)));
        
        %     switch lower(SAM.optim.cost.stat.stat)
        %         case 'bic'
        %             fprintf(1,'BIC = %.2f, Chi^2 = %.2f \n',cost,altCost)
        %         case 'chisquare'
        %             fprintf(1,'Chi^2 = %.2f, BIC = %.2f \n',cost,altCost)
        %     end
        

    end
    
end

% Shut down the parallel pool
% =========================================================================
if doStartParCluster
    delete(myPool);
end

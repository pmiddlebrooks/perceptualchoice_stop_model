function job_simulate_data(subject,model,architecture,dt,trialVar,optimScope,fileStr,nSim)

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
saveDir             = fileStr.save;
nameFVal            = 'allFValAndBestX_%sTrials_model%.3d.txt';
nameModelSAMSpec  	= 'SAM_%sTrials_model%.3d.mat';
nameModelSAMGen    	= 'SAM_%sTrials.mat';
nameUserSpecX       = 'userSpecX_%sTrials_model%.3d.txt';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Simulate the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop over subjects, architectures, and models
warning('off','stats:dataset:ModifiedVarnames')
for iSubject = subject
    for iArchitecture = 1:nArchitecture
        
        
        for iModel = model
            
            % Display progress
            fprintf(1,'Working on subject %d, architecture %s, and model %d. \n',iSubject,architecture{iArchitecture},iModel);
            
            % Load model-specific SAM
            ds = dataset('File',fullfile(sprintf(rootDir,iSubject,dt,trialVarStr,architecture{iArchitecture}),sprintf(nameFVal,optimScope,iModel)));
            
            % Load SAM
            load(fullfile(sprintf(rootDir,iSubject,dt,trialVarStr,architecture{iArchitecture}),ds.FileNameSAM{1}),'SAM');
            
            % Change simulation number to 10000 (default in
            % job_spec_sam_general is 2500
            SAM.sim.n = nSim;
            
            % Extract optimized X
            iBestX = cell2mat(cellfun(@(in1) ~isempty(regexp(in1,'^BestX.*', 'once')),ds.Properties.VarNames,'Uni',0));
            X = double(ds(1,iBestX));
            
            % Get model predictions and costs
            [cost,altCost,prd] = sam_cost(X,SAM);
            
            simData = job_preproc_simulated_data(prd, iSubject, fileStr.preproc);
%             % Save predicted data as observed
%             SAM.optim.obs = prd;
            
            % Make directory if it doesn't exist
            saveDirFull = sprintf(saveDir,iSubject,dt,trialVarStr,architecture{iArchitecture});
            if exist(saveDirFull,'dir') ~= 7
                mkdir(saveDirFull)
            end
            
%             % Change SAM working directory to this simulation date's
%             % directory (used in job_optim_x0.m and job_optim_x.m to save iteration and final log files)
%             SAM.io.workDir = saveDirFull;

            % Save previous best fit parameters as user-specified
            % parameters (see job_spec_x0base.m)
            save(fullfile(saveDirFull,sprintf(nameUserSpecX,optimScope,iModel)),'-ascii','-double','-tabs','X');
            
            % Save the SAM structure with predictions as new observations,
            % both for the geneal SAM structure and the specific model SAM
            save(fullfile(saveDirFull,sprintf(nameModelSAMSpec,optimScope,iModel)),'SAM');
            save(fullfile(saveDirFull,sprintf(nameModelSAMGen,optimScope)),'SAM');
        end
    end
end
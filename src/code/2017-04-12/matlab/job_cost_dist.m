% job_cost_dist
%
% This script is called by ACCRE's SLURM (Simple Linux Utility for Resource
% Management) and aims to simulate data reseeding the random number
% generator to produce a distribution of cost variables (chi-2 and BIC
% values)
%
% DESCRIPTION
% This script contains all the details for the job to run
%
% .........................................................................


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. GET ENVIRONMENT TO DETERMINE PATHS AND SCRIPT-SPECIFICS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
accreRoot = '/gpfs22';
accreHome = '/home/middlepg';
accreScratch = '/scratch/middlepg';
if isdir(fullfile(accreScratch))
    matRoot = fullfile(accreRoot,accreHome,'m-files/');
    modelRoot = fullfile(accreScratch,'perceptualchoice_stop_model/');
    environ = 'accre';
else
    matRoot = '/Volumes/HD-1/Users/paulmiddlebrooks/matlab/';
    modelRoot = '/Volumes/HD-1/Users/paulmiddlebrooks/perceptualchoice_stop_model/';
    environ = 'local';
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. PROCESS INPUTS, DEFINE VARIABLES, ADD ACCESS TO PATHS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1.1. Process inputs
% =========================================================================

switch environ
    case 'accre'
        % Get subject index
        iSubj                 = str2double(getenv('subject'));
        
        % Get architecture
        iModelArch             = getenv('modelArch');
        
        % Get model index
        iModel                = str2double(getenv('iModel'));
        
        % Get starting point index which iteration is this?
        iStartVal 			= str2double(getenv('iStartVal'));
    otherwise
        % otherwise those variables are set by the script
        % exe_job_optim_x0.m
end

% Get path string
pathStr               = getenv('pathStr');

% Get root directory
rootDir               = getenv('rootDir');

% Get save directory
saveDir               = getenv('saveDir');

% Get path string for initial parameters
% pathStrInitParam      = getenv('pathStrInitParam');


% Get time step size
dt                    = str2double(getenv('dt'));

% Trial-to-trial variability
trialVar              = getenv('trialVar');

% Get optimScope
optimScope              = getenv('optimScope');

% Number of starting points
nModleSimulations           = str2double(getenv('nModelSim'));

% Get number of processors per node
nProcessors 	      = str2double(getenv('processorsPerNode'));



% Get job ID and parse digits
jobID                 = getenv('SLURM_JOBID');
jobIDDigits           = regexp(jobID,'\d','match');
jobIDDigits           = [jobIDDigits{:}];




% 1.2. Define dynamic variables
% =========================================================================
timeStr   = datestr(now,'yyyy-mm-dd-THHMMSS');

% 1.3. Add paths
% ========================================================================

addpath(genpath(fullfile(matRoot,'sam')));
addpath(genpath(fullfile(matRoot,'matlab_code_bbz')));
addpath(genpath(fullfile(matRoot,'matlab_file_exchange_tools')));
addpath(genpath(fullfile(matRoot,'cmtb')));


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. SAVE AND RUN JOB
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2.1. Load the SAM file and Get best fitting paramters for this model
nameFVal            = 'allFValAndBestX_%sTrials_model%.3d.txt';
ds = dataset('File',fullfile(sprintf(rootDir,iSubj,dt,trialVar,iModelArch),sprintf(nameFVal,optimScope,iModel)));

% Load SAM
load(fullfile(sprintf(rootDir,iSubj,dt,trialVar,iModelArch),ds.FileNameSAM{1}),'SAM');

% Extract optimized X
iBestX = cell2mat(cellfun(@(in1) ~isempty(regexp(in1,'^BestX.*', 'once')),ds.Properties.VarNames,'Uni',0));
X = double(ds(1,iBestX));



% 2.2. Add details for cluster computing
% =========================================================================
SAM.compCluster.nProcessors   = nProcessors;

% 2.3. Adjust the number of starting points
% =========================================================================
SAM.optim.nStartPoint = nModleSimulations;

% 2.4. Add details for logging
% =========================================================================
fNameIterLog          = sprintf('iterLog_cost_%sTrials_model%.3d_startSet%.3d_started%s.mat',optimScope,iModel,iStartVal,timeStr);
fNameFinalLog         = sprintf('finalLog_cost_%sTrials_model%.3d_startSet%.3d_started%s.mat',optimScope,iModel,iStartVal,timeStr);

% Iteration log file
fitLog.iterLogFile    = fullfile(sprintf(saveDir,iSubj,dt,trialVar,iModelArch),fNameIterLog);

% Iteration lof frequency
fitLog.iterLogFreq    = 50;

% Final log file
fitLog.finalLogFile   = fullfile(sprintf(saveDir,iSubj,dt,trialVar,iModelArch),fNameFinalLog);

SAM.optim.log         = fitLog;

clear fitLog

% 2.5. Save the init-param SAM file
% =========================================================================
% save(sprintf(pathStrInitParam,iSubj,dt,trialVar,iModelArch,optimScope,iModel));

% 2.6. Loop over all candidate starting points and compute cost, and save every iterFreq
% =========================================================================================================================

iterLogFile           = SAM.optim.log.iterLogFile;
iterLogFreq           = SAM.optim.log.iterLogFreq;
finalLogFile          = SAM.optim.log.finalLogFile;

history = nan(nModleSimulations + 1,numel(SAM.optim.x0Base) + 2);

if exist(sprintf(saveDir,iSubj,dt,trialVar,iModelArch),'dir') ~= 7
    mkdir(sprintf(saveDir,iSubj,dt,trialVar,iModelArch))
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% #. START THE PARALLEL POOL
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% In unique directory to prevent collision of parallel jobs
% e.g. see: http://www.mathworks.com/matlabcentral/answers/97141-why-am-i-unable-to-start-a-local-matlabpool-from-multiple-matlab-sessions-that-use-a-shared-preferen
% c = parcluster();
% if isfield(SAM,'compCluster')
%     c.NumWorkers = SAM.compCluster.nProcessors;
% end
% [~,homeDir] = system('echo $HOME');
% homeDir = strtrim(homeDir);
% release = version('-release');
% tempDir = fullfile(homeDir,'.matlab','local_cluster_jobs',release);
% if exist(tempDir) ~= 7
%     mkdir(tempDir)
% end
% t = tempname(tempDir);
% mkdir(t);
% c.JobStorageLocation=t;
% tWait = 1+60*rand();
% pause(tWait);
% myPool = parpool(c);

for iter = 1:nModleSimulations
    tic
    disp(sprintf('This is iter %.3d',iter))
    
    % Re-seed the random number generatory each time:
    SAM.sim.rng.id = rng('shuffle');
    
    [cost,altCost] = sam_cost(X,SAM);
    
    history(iter,1) = cost;
    history(iter,2) = altCost;
    
    if iterLogFreq*round(double(iter)/iterLogFreq) == iter;
        save(iterLogFile,'history');
    end
    toc
end

% Save the final log file
save(finalLogFile,'history');

% Remove iteration log file
delete(iterLogFile);

% Shut down the parallel pool
% =========================================================================
% delete(myPool);

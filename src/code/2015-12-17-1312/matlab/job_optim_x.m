% pbs_job_optim_x
%
% This script is called by ACCRE'sSLURM (Simple Linux Utility for Resource Management) and aims to optimize a set of parameters for a model, using
% settings specified in SAM
% 
%
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

% Get starting point index
iStartVal 			= str2double(getenv('iStartVal'));

    otherwise
        % otherwise those variables are set by the script
        % exe_job_optim_x0.m
end

% Get path string
pathStr               = getenv('pathStr');

% Get time step size
dt                    = str2double(getenv('dt'));

% Trial-to-trial variability
trialVar              = getenv('trialVar');

% Get optimScope
optimScope              = getenv('optimScope');

% Number of starting points
nStartPoint           = str2double(getenv('nStartPoint'));

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
% =========================================================================

addpath(genpath(fullfile(matRoot,'sam')));
addpath(genpath(fullfile(matRoot,'matlab_code_bbz')));
addpath(genpath(fullfile(matRoot,'matlab_file_exchange_tools')));
addpath(genpath(fullfile(matRoot,'cmtb')));


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 2. RUN AND SAVE JOB
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% 2.1. Load the SAM file
% =========================================================================
load(sprintf(pathStr,iSubj,dt,trialVar,iModelArch,optimScope,iModel));

% 2.2. Add details for cluster computing
% =========================================================================
SAM.compCluster.nProcessors   = nProcessors;

SAM.compCluster.jobID         = jobIDDigits;

% 2.3. Add details for logging
% =========================================================================
fNameIterLog                  = sprintf('iterLog_%sTrials_model%.3d_startVal%.3d_started%s.mat',optimScope,iModel,iStartVal,timeStr);
fNameFinalLog                 = sprintf('finalLog_%sTrials_model%.3d_startVal%.3d_started%s.mat',optimScope,iModel,iStartVal,timeStr);

% Iteration log file
fitLog.iterLogFile            = fullfile(SAM.io.workDir,fNameIterLog);

% Iteration lof frequency
fitLog.iterLogFreq            = 50;

% Final log file
fitLog.finalLogFile           = fullfile(SAM.io.workDir,fNameFinalLog);

SAM.optim.log                 = fitLog;

% 2.4. Optimize the fit to the data, starting from parameters corresponding to iStartVal
% =========================================================================
SAM                           = sam_optim(SAM,iStartVal);
fVal                          = SAM.estim.fVal;
X                             = SAM.estim.X;
fValX                         = [fVal,X];

% 2.5. Optimize the fit to the data, starting from parameters corresponding to iStartVal
% =========================================================================
fNameSAM                      = sprintf('SAM_%sTrials_model%.3d_startVal%.3d_exit%d_started%s.mat',optimScope,iModel,iStartVal,SAM.estim.exitFlag,timeStr);
fNameX                        = sprintf('fValAndBestX_%sTrials_model%.3d_startVal%.3d_exit%d_started%s.txt',optimScope,iModel,iStartVal,SAM.estim.exitFlag,timeStr);

save(fullfile(SAM.io.workDir,fNameSAM),'SAM');
save(fullfile(SAM.io.workDir,fNameX),'fValX','-ascii','-double');
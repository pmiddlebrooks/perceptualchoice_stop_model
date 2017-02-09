% pbs_job_optim_x
%
% This script is called by ACCRE's portable batch system (PBS) and aims to optimize a set of parameters for a model, using
% settings specified in SAM
% 
%
%  
% DESCRIPTION 
% This script contains all the details for the job to run
%  
% ......................................................................... 
% Bram Zandbelt, bramzandbelt@gmail.com 
% $Created : Mon 09 Sep 2013 13:07:49 CDT by bram 
% $Modified: Sat 21 Sep 2013 12:24:04 CDT by bram

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 1. PROCESS INPUTS, DEFINE VARIABLES, ADD ACCESS TO PATHS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% 1.1. Process inputs
% =========================================================================

% Get path string
pathStr   			= getenv('pathStr');

% Get subject index
iSubj     			= str2double(getenv('subject'));

% Get time step size
dt        			= str2double(getenv('dt'));

% Trial-to-trial variability
trialVar  			= getenv('trialVar');

% Get architecture
modelArch 			= getenv('modelArch');

% Get optimScope
optimScope  			= getenv('optimScope');

% Get model index
iModel    			= str2double(getenv('iModel'));

% Get starting point index
iStartVal 			= str2double(getenv('iStartVal'));

% Get number of processors per node
nProcessors 		= str2double(getenv('processorsPerNode'));

% 1.2. Define dynamic variables
% =========================================================================
timeStr   = datestr(now,'yyyy-mm-dd-THHMMSS');

% 1.3. Add paths
% =========================================================================
addpath(genpath('/home/middlepg/m-files/sam/'));
addpath(genpath('/home/middlepg/m-files/matlab_code_bbz/'));
addpath(genpath('/home/middlepg/m-files/matlab_file_exchange_tools/'));
addpath(genpath('/home/middlepg/m-files/cmtb/'));

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 2. RUN AND SAVE JOB
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% 2.1. Load the SAM file
% =========================================================================
load(sprintf(pathStr,iSubj,dt,trialVar,modelArch,optimScope,iModel));

% 2.2. Add details for cluster computing
% =========================================================================
SAM.compCluster.nProcessors   = nProcessors;

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
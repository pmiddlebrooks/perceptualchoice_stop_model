% pbs_job_optim_x0
%
% This script is called by ACCRE's portable batch system (PBS) and aims to find optimal initial parameters, before model
% optimization starts
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
pathStr               = getenv('pathStr');

% Get path string for initial parameters
pathStrInitParam      = getenv('pathStrInitParam');

% Get subject index
% iSubj                 = str2double(getenv('subject'));

% Get time step size
dt                    = str2double(getenv('dt'));

% Trial-to-trial variability
trialVar              = getenv('trialVar');

% Get architecture
% modelArch             = getenv('modelArch');
modelArch             = iModelArch;

% Get optimScope
optimScope              = getenv('optimScope');

% Get model index
% iModel                = str2double(getenv('iModel'));

% Number of starting points
nStartPoint           = str2double(getenv('nStartPoint'));

% Get number of processors per node
nProcessors 	      = str2double(getenv('processorsPerNode'));


% 1.2. Define dynamic variables
% =========================================================================
timeStr   = datestr(now,'yyyy-mm-dd-THHMMSS');

% 1.3. Add paths
% =========================================================================
addpath(genpath('/Volumes/HD-1/Users/paulmiddlebrooks/perceptual_choice_stop_model/sam/'));
addpath(genpath('/Volumes/HD-1/Users/paulmiddlebrooks/perceptual_choice_stop_model/matlab_code_bbz/'));
addpath(genpath('/Volumes/HD-1/Users/paulmiddlebrooks/perceptual_choice_stop_model/matlab_file_exchange_tools/'));
addpath(genpath('/Volumes/HD-1/Users/paulmiddlebrooks/perceptual_choice_stop_model/cmtb/'));


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 2. SAVE AND RUN JOB
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% 2.1. Load the general SAM file
% 
sprintf(pathStr,iSubj,dt,trialVar,modelArch,optimScope)
load(sprintf(pathStr,iSubj,dt,trialVar,modelArch,optimScope));

% 2.2. Add details for cluster computing
% =========================================================================
SAM.compCluster.nProcessors   = nProcessors;

% 2.3. Adjust the number of starting points and specify model-specific SAM structure
% =========================================================================
SAM.optim.nStartPoint = nStartPoint;
SAM                   = sam_spec_job_specific(SAM,iModel);

% 2.4. Add details for logging
% =========================================================================
fNameIterLog          = sprintf('iterLog_initParam_%sTrials_model%.3d_started%s.mat',optimScope,iModel,timeStr);
fNameFinalLog         = sprintf('finalLog_initParam_%sTrials_model%.3d_started%s.mat',optimScope,iModel,timeStr);

% Iteration log file
fitLog.iterLogFile    = fullfile(SAM.io.workDir,fNameIterLog);

% Iteration lof frequency
fitLog.iterLogFreq    = 50;

% Final log file
fitLog.finalLogFile   = fullfile(SAM.io.workDir,fNameFinalLog);

SAM.optim.log         = fitLog;

clear fitLog

% 2.5. Save the init-param SAM file
% =========================================================================
save(sprintf(pathStrInitParam,iSubj,dt,trialVar,modelArch,optimScope,iModel));

% 2.6. Loop over all candidate starting points and compute cost, and save every iterFreq
% =========================================================================================================================

iterLogFile           = SAM.optim.log.iterLogFile;
iterLogFreq           = SAM.optim.log.iterLogFreq;
finalLogFile          = SAM.optim.log.finalLogFile;

history = nan(nStartPoint + 1,numel(SAM.optim.x0Base) + 2);


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% #. START THE PARALLEL POOL
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% In unique directory to prevent collision of parallel jobs
% e.g. see: http://www.mathworks.com/matlabcentral/answers/97141-why-am-i-unable-to-start-a-local-matlabpool-from-multiple-matlab-sessions-that-use-a-shared-preferen
c = parcluster();
if isfield(SAM,'compCluster')
	c.NumWorkers = SAM.compCluster.nProcessors;
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

for iter = 1:nStartPoint
  
  disp(sprintf('This is iter %.3d',iter))
  
  [cost,altCost] = sam_cost(SAM.optim.x0(iter,:),SAM);
  
  history(iter,1) = cost;
  history(iter,2) = altCost;
  history(iter,3:end) = SAM.optim.x0(iter,:);
  
  if iterLogFreq*round(double(iter)/iterLogFreq) == iter;
      bestX0 = sortrows(history,1);
      bestX0 = bestX0(1:iterLogFreq,:);
      save(iterLogFile,'history','bestX0');
  end
end

% Save the final log file
save(finalLogFile,'history','bestX0');

% Remove iteration log file
delete(iterLogFile);

% Shut down the parallel pool
% =========================================================================
delete(myPool);

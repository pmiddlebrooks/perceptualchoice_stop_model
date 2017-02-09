function job_spec_sam_general(subj,dt,trialVar,optimScope,choiceMech,stopMech,rngStruct,fileStr)
%
% INPUTS
% subj          vector of subject indices
% dt            time step, in ms
% trialVar      true or false, influences t0, z0, and eta
% optimScope    'go', 'stop', or 'all'
% choiceMech    'race', 'ffi', or 'li'
% stopMech      'race', 'bi', or 'li'
% rngStruct     random number generator structure (e.g. rngStruct = rng('shuffle'))
% fileStr       struct containing 
%
% fileStr.preprocDataDir    = '/scratch/zandbeb/multichoice_stop_model/data/2014-05-23-1515/preproc01/';
% fileStr.rawDataDir        = '/scratch/zandbeb/multichoice_stop_model/data/raw/';
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
  case {'stop','all'}
    modelArch = sprintf('c%s_i%s',choiceMech,stopMech);
end

for iSubj = subj

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. INPUT/OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1.1. Define file/directory strings
% =========================================================================

% File string to root directory
io.rootDir                    = fullfile(fileparts(fileStr.preprocDataDir),sprintf('subj%.2d/',iSubj));

% File string to behavioral data file
io.behavFile                  = fullfile(io.rootDir,sprintf('data_subj%.2d.mat',iSubj));

% File string to working directory
io.workDir                    = fullfile(io.rootDir,sprintf('dt%d',dt),trialVarStr,modelArch);

% Path to directory where raw data reside
io.rawDataDir                 = fullfile(fileparts(fileStr.rawDataDir));

% 1.2. Make directories, if they don't exist
% =========================================================================
if exist(io.rootDir,'dir') ~= 7
  mkdir(io.rootDir)
end

if exist(io.workDir,'dir') ~= 7
  mkdir(io.workDir)
end

% 1.3. Put io into SAM struct
% =========================================================================
SAM.io                        = io;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. EXPERIMENTAL DETAILS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2.1. Numbers
% =========================================================================

% Number of stimuli and responses, across conditions
% -------------------------------------------------------------------------
expt.nStm                     = {[2 1] [2 1] [2 1]};
expt.nRsp                     = {[2 1] [2 1] [2 1]};

% Number of experimental conditions
% -------------------------------------------------------------------------
expt.nCnd                     = numel(expt.nStm);                          

% Number of stop-signal delays
% -------------------------------------------------------------------------
% To find out, load data_subj%d.mat, and evaluate 'numel(unique(data.ssd))'
load(io.behavFile)
expt.nSsd                     = numel(unique(data.ssd));

% 2.2. Timing (in ms)
% =========================================================================

% Trial duration
% -------------------------------------------------------------------------
expt.trialDur                 = 2500;                                             

% Stimulus onset (go-signal, stop-signal)
% -------------------------------------------------------------------------
expt.stmOns                   = [0 nan];

% Stimulus duration
% -------------------------------------------------------------------------
expt.stmDur                   = [2000 100];

% 2.3. Put expt into SAM struct
% =========================================================================
SAM.expt                      = expt;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. SIMULATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 3.1. Goal
% =========================================================================
% The goal of a simulation can be the optimize parameters to account for
% the data ('optimize') or to explore model predictions ('explore')

sim.goal                      = 'optimize';

% 3.2. Scope
% =========================================================================
sim.scope                     = optimScope;
sim.n                         = 2500;

% 3.3. Simulation functions
% =========================================================================
% These functions are used to simulate indivdiual trials and the experiment
% as a whole (consisting of multiple trials)
sim.fun.trial                 = @sam_sim_trial_mex;                         
sim.fun.expt                  = @sam_sim_expt;

% 3.4. Random number generator
% =========================================================================

% 3.4.1. Seeding stage
% -------------------------------------------------------------------------
% The random number generator can be seeded at two points:
% 'sam_optim'     - this results in seeding the rng once
% 'sam_sim_expt'  - this results in seeding the rng every time new model
%                   predictions are simulated
sim.rng.stage                 = 'sam_sim_expt';

% 3.4.2. Identified struct
% -------------------------------------------------------------------------
sim.rng.id                    = rngStruct;

% 3.5. Time windows
% =========================================================================
% These are only relevant when model predictions are explored (sim.goal =
% 'explore') and specify the time windows for which accumulator dynamics
% should be computed.

% 3.5.1. Go-signal-aligned
% -------------------------------------------------------------------------
sim.tWindow.go                = [0 2000];

% 3.5.2. Stop-signal-aligned
% -------------------------------------------------------------------------
sim.tWindow.stop              = [0 2000];

% 3.5.3. Response-aligned
% -------------------------------------------------------------------------
sim.tWindow.resp              = [-500 0];

SAM.sim                       = sim;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 4. MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 4.1. General
% =========================================================================

% Tag specifying model architecture
general.modelCatTag           = modelArch;

% Names of accumulator classes (corresponding to expt.nStm and expt.nRsp)
general.classNames            = {'GO','STOP'};

% 4.1.1. Put general into SAM struct
% -------------------------------------------------------------------------
SAM.model.general             = general;

% 4.2. Accumulation specifics
% =========================================================================

% 4.2.1. Bounds
% -------------------------------------------------------------------------

% Lower bound on activation
accum.zLB                     = 0;

% 4.2.2. Timing
% -------------------------------------------------------------------------
% Time step
accum.dt                      = dt;

% Time scale
accum.tau                     = 1;

% Time window of accumulation
accum.window                  = [0 2000];

% STOP accumulation duration
% Determine if expt.stmDur should be overridden for STOP unit
% 'signal'  - STOP activation accumulates for duration of stop-signal 
% 'trial'   - STOP activation accumulates for duration of entire trial
accum.durSTOP                 = 'trial';

% 4.2.3. Variability
% -------------------------------------------------------------------------
% Non-decision time and starting point may vary between trials (uniform 
% distribution)

if trialVar
  accum.randomT0              = true;
  accum.randomZ0              = true;
else
  accum.randomT0              = false;
  accum.randomZ0              = false;
end

% 4.2.4. Put accum into SAM struct
% -------------------------------------------------------------------------
SAM.model.accum               = accum;

% 4.3. Parameter categories
% =========================================================================

% 4.3.1. Names, numbers, and indices
% -------------------------------------------------------------------------

% Names
XCat.name                     = {'z0', ...    % Starting point
                                 'zc', ...    % Threshold
                                 'v', ...     % Rate, target unit
                                 've', ...    % Rate, non-target unit
                                 'eta', ...   % Variability in rate, between trials
                                 't0', ...    % Non-decision time
                                 'se', ...    % Extrinsic noise
                                 'si', ...    % Intrinsic noise
                                 'k', ...     % Leakage constant
                                 'wliw', ...  % Lateral inhibition weight, within class (e.g. between GO units)
                                 'wlib', ...  % Lateral inhibition weight, between classes (e.g. STOP to GO)
                                 'wffiw'};    % Feed-forward inhibition weight, within class (e.g. between go signals)

% Numbers
XCat.n                        = numel(XCat.name);

% Indices
XCat.i.iZ0                    = find(strncmpi(XCat.name,'z0',2));
XCat.i.iZc                    = find(strncmpi(XCat.name,'zc',2));
XCat.i.iV                     = find(strncmpi(XCat.name,'v',2));
XCat.i.iVe                    = find(strncmpi(XCat.name,'ve',2));
XCat.i.iEta                   = find(strncmpi(XCat.name,'eta',3));
XCat.i.iT0                    = find(strncmpi(XCat.name,'t0',2));
XCat.i.iSe                    = find(strncmpi(XCat.name,'se',2));
XCat.i.iSi                    = find(strncmpi(XCat.name,'si',2));
XCat.i.iK                     = find(strncmpi(XCat.name,'k',2));
XCat.i.iWliw                  = find(strncmpi(XCat.name,'wliw',4));
XCat.i.iWlib                  = find(strncmpi(XCat.name,'wlib',4));
XCat.i.iWffiw                 = find(strncmpi(XCat.name,'wffiw',5));

% 4.3.2. Inclusion/exclusion and class-specificity
% -------------------------------------------------------------------------

switch lower(choiceMech)
  case 'race'
    switch lower(stopMech)
      case {'race','bi'}
        if trialVar
          % 
          XCat.included       = logical([1 1 1 1 1 1 0 1 1 0 0 0]);
          XCat.classSpecific  = logical([1 1 1 1 1 1 0 0 1 0 0 0]);
        else
          XCat.included       = logical([0 1 1 1 0 1 0 1 1 0 0 0]);
          XCat.classSpecific  = logical([0 1 1 1 0 1 0 0 1 0 0 0]);
        end
      case 'li'
        if trialVar
          % 
          XCat.included       = logical([1 1 1 1 1 1 0 1 1 0 1 0]);
          XCat.classSpecific  = logical([1 1 1 1 1 1 0 0 1 0 1 0]);
        else
          XCat.included       = logical([0 1 1 1 0 1 0 1 1 0 1 0]);
          XCat.classSpecific  = logical([0 1 1 1 0 1 0 0 1 0 1 0]);
        end
      otherwise
        % e.g. when only go trials are simulated
        if trialVar
          % 
          XCat.included       = logical([1 1 1 1 1 1 0 1 1 0 0 0]);
          XCat.classSpecific  = logical([1 1 1 1 1 1 0 0 1 0 0 0]);
        else
          XCat.included       = logical([0 1 1 1 0 1 0 1 1 0 0 0]);
          XCat.classSpecific  = logical([0 1 1 1 0 1 0 0 1 0 0 0]);
        end
    end
  case'ffi'
    % NOTE: Consider making noise extrinsic in ffi models and setting se as scale factor
    switch lower(stopMech)
      case {'race','bi'}
        if trialVar
          % 
          XCat.included       = logical([1 1 1 1 1 1 0 1 1 0 0 1]);
          XCat.classSpecific  = logical([1 1 1 1 1 1 0 0 1 0 0 1]);
        else
          XCat.included       = logical([0 1 1 1 0 1 0 1 1 0 0 1]);
          XCat.classSpecific  = logical([0 1 1 1 0 1 0 0 1 0 0 1]);
        end
      case 'li'
        if trialVar
          % 
          XCat.included       = logical([1 1 1 1 1 1 0 1 1 0 1 1]);
          XCat.classSpecific  = logical([1 1 1 1 1 1 0 0 1 0 1 1]);
        else
          XCat.included       = logical([0 1 1 1 0 1 0 1 1 0 1 1]);
          XCat.classSpecific  = logical([0 1 1 1 0 1 0 0 1 0 1 1]);
        end
      otherwise
        % e.g. when only go trials are simulated
        if trialVar
          % 
          XCat.included       = logical([1 1 1 1 1 1 0 1 1 0 0 1]);
          XCat.classSpecific  = logical([1 1 1 1 1 1 0 0 1 0 0 1]);
        else
          XCat.included       = logical([0 1 1 1 0 1 0 1 1 0 0 1]);
          XCat.classSpecific  = logical([0 1 1 1 0 1 0 0 1 0 0 1]);
        end
    end
  case 'li'
    switch lower(stopMech)
      case {'race','bi'}
        if trialVar
          XCat.included       = logical([1 1 1 1 1 1 0 1 1 1 0 0]);
          XCat.classSpecific  = logical([1 1 1 1 1 1 0 0 1 1 0 0]);
        else
          XCat.included       = logical([0 1 1 1 0 1 0 1 1 1 0 0]);
          XCat.classSpecific  = logical([0 1 1 1 0 1 0 0 1 1 0 0]);
        end
      case 'li'
        if trialVar
          XCat.included       = logical([1 1 1 1 1 1 0 1 1 1 1 0]);
          XCat.classSpecific  = logical([1 1 1 1 1 1 0 0 1 1 1 0]);
        else
          XCat.included       = logical([0 1 1 1 0 1 0 1 1 1 1 0]);
          XCat.classSpecific  = logical([0 1 1 1 0 1 0 0 1 1 1 0]);
        end
      otherwise
        % e.g. when only go trials are simulated
        if trialVar
          XCat.included       = logical([1 1 1 1 1 1 0 1 1 1 0 0]);
          XCat.classSpecific  = logical([1 1 1 1 1 1 0 0 1 1 0 0]);
        else
          XCat.included       = logical([0 1 1 1 0 1 0 1 1 1 0 0]);
          XCat.classSpecific  = logical([0 1 1 1 0 1 0 0 1 1 0 0]);
        end
    end
end

% Default parameter values for excluded categories
XCat.valExcluded              = [0 realmax 0 0 0 0 0 0 0 0 0 0];

% 4.3.3. Constraints
% -------------------------------------------------------------------------

% Hard lower and upper bounds on parameter values
XCat.hardLB                   = [SAM.model.accum.zLB 0   0   0   0   0   0   0   -Inf     -Inf -Inf -Inf];
XCat.hardUB                   = [Inf Inf Inf Inf Inf Inf Inf Inf -realmin -realmin -realmin -realmin];

% Specify values for hard lower and upper bounds
% Having specified initials parameter values (x0, i.e. rough estimates),
% parameter bounds are specified depending on x0:
% - Additive bounds:        x0 +/- XCat.additive
% - Multiplicative bounds:  x0 +/- x0.*XCat.multiplicative

XCat.additive                 = [5   0   0   0   0   0   0   0   0.005 0.5 0.5 1];
XCat.multiplicative           = [1   1   1   1   1   1   1   1   0   0 0 0];

% 4.3.4. Scaling
% -------------------------------------------------------------------------
% Determine index and value of the scaling parameter

switch lower(choiceMech)
  case 'race'
    XCat.scale.iX             = find(strncmpi(XCat.name,'si',2));
    XCat.scale.val            = 1;
  case'ffi'
    % NOTE: Consider making noise extrinsic in ffi models and setting se as scale factor
    XCat.scale.iX             = find(strncmpi(XCat.name,'si',2));
    XCat.scale.val            = 1;
  case 'li'
    XCat.scale.iX             = find(strncmpi(XCat.name,'si',2));
    XCat.scale.val            = 1;
end

% 4.3.5. Put accum into SAM struct
% -------------------------------------------------------------------------
SAM.model.XCat                = XCat;

% 4.4. Features
% =========================================================================
% The features matrix indicates how model parameters are supposed to vary
% across experimental factors (first dimension), parameter category (second
% dimension), and accumulator classes (third dimension):
% - experimental factors
%   * stimuli
%   * responses
%   * conditions
% - parameter categories
%   * starting point
%   * threshold
%   * rate, target unit
%   * rate, non-target unit
%   * rate variability
%   * non-decision time
%   * extrinsic noise
%   * intrinsic noise
%   * leakage constant
%   * lateral inhibition weight, within class
%   * lateral inhibition weight, between classes
%   * feed-forward inhibition weight, within class
% - accumulator classes
%   * GO
%   * STOP

% 4.4.1. Define feature matrix
% -------------------------------------------------------------------------
features                      = false(3,XCat.n,2);

% zc_{Go} is allowed to vary between responses; zc_{Go}, v_{Go}, ve_{Go}, 
% t0_{Go}, wli_{Go}, and wffi_{Go} are allowed to vary between conditions
features(2,XCat.i.iZ0,1)      = true; %pgm
features(2,XCat.i.iZc,1)      = true;
features(3,XCat.i.iZc,1)      = true;
features(3,XCat.i.iV,1)       = true;
features(3,XCat.i.iVe,1)      = true;
features(3,XCat.i.iT0,1)      = true;
features(3,XCat.i.iWliw,1)    = true;
features(3,XCat.i.iWffiw,1)   = true;

% 4.4.2. Put features into SAM struct
% -------------------------------------------------------------------------
SAM.model.features            = features;

% 4.5. Matrices
% =========================================================================
% Specify model matrices, including
% - endogenous connectivity (connectivity withing and between accumulators)
% - exogenous connectivity (stimulus-response mappings, feed-forward
%   inhibition)
% - consequences of threshold crossings
%   * trial termination rules
%   * blocked input rules
%   * lateral inhibition rules
%
% Model matrices are defined based on choice and stop mechanisms.

SAM.model.mat                 = sam_spec_conn_mat(SAM,choiceMech,stopMech);

% 4.6. Variants
% =========================================================================
% Specify all possible models, given the feature matrix

SAM.model.variants.tree       = sam_spec_potential_models(SAM);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 6. OPTIMIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 6.1. Solver
% =========================================================================

% Type
optim.solver.type             = 'fminsearchcon';

% Options
optim.solver.opts             = sam_get_solver_opts('fminsearchcon');

% 6.2. Cost function
% =========================================================================

% Function
optim.cost.fun                = @sam_cost;

% Cost statistic
% - 'chisquare', Pearson's chi-squared
% - 'bic', Bayesian Information Criterion
optim.cost.stat.stat          = 'chisquare';

% Cumulative probabilities, used to compute response time bins
optim.cost.stat.cumProb       = [.1 .3 .5 .7 .9];

% Minimum response time bin size
% - for trial categories with fewer trials than the first element, the
% model will try to fit response probability, not response time
% - for trial categories with a greater number of trials than the first
% element, but less trials than the second element will be splitted into
% two bins (i.e. using the median as bin edge)
optim.cost.stat.minBinSize    = [10,30];

% 6.3. Number of independent starting points
% =========================================================================
% This is the number of initial, randomly-selected starting points from
% which the optimization algorithm will start.
optim.nStartPoint             = 20;

% 6.4. Put optim in SAM
% =========================================================================
SAM.optim                     = optim;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 7. COMPUTER CLUSTER VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of processors, only relevant when running job_spec_x0base.m on
% ACCRE

SAM.compCluster.nProcessors   = 4;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 8. SAVE SAM STRUCTURE TO WORKING DIRECTORY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save(fullfile(SAM.io.workDir,sprintf('SAM_%sTrials.mat',SAM.sim.scope)),'SAM');
keep SAM subj dt optimScope modelArch trialVarStr trialVar choiceMech stopMech fileStr rngStruct

cd(SAM.io.workDir)

end

setenv('pathStr', '/Volumes/HD-1/Users/paulmiddlebrooks/perceptualchoice_stop_model/data/2015-11-28-1842/preproc01/subj%.2d/dt%d/%s/%s/SAM_%sTrials.mat');
setenv('pathStrInitParam','/Volumes/HD-1/Users/paulmiddlebrooks/perceptualchoice_stop_model/data/2015-11-28-1842/preproc01/subj%.2d/dt%d/%s/%s/SAM_initParam_%sTrials_model%.3d.mat');

% Optimization scope
setenv('optimScope','go');

% Time step
setenv('dt','20');

% Trial variability tag (for starting point, non-decision time, and rate)
setenv('trialVar','trialvar');

% Model architectures
% setenv('modelArch','crace');
modelArch = {'crace', 'cli', 'cffi'};
modelArch = {'crace'};

% Model variants
% setenv('model','9');
model = [49,99];

% Number of starting points
% nStartPoint=2000;
setenv('nStartPoint','2000');

% Subject
% setenv('subject','1'); %[1 2];
subject = [1 2];

% SBATCH Settings
% WALLTIME="15:00:00";
% PPN="4";
% MEM="12000";
setenv('processorsPerNode','4');

% modelArch             = cellstr(getenv('modelArch'));
% subject                 = str2double(getenv('subject'));

for iModArch = 1 : length(modelArch)
    iModelArch = modelArch{iModArch};
   
    for iMod = 1 : length(model)
        iModel = model(iMod);
    
        for iSub = 1 : length(subject)
            iSubj = subject(iSub);
        
            job_optim_x0;
            
        end
    end
end

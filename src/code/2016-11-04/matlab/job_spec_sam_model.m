function job_spec_sam_model(subj,dt,trialVar,optimScope,choiceMech,stopMech,model,fileStr)
%
% INPUTS
% subj          vector of subject indices
% dt            time step, in ms
% trialVar      true or false, influences t0, z0, and eta
% optimScope      'go','stop', or 'all'
% choiceMech    'race', 'ffi', or 'li'
% stopMech      'race', 'bi', or 'li'
% iModel        vector of model indices

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

rootDir           = fileStr.root;
nameInitLog       = 'iterLog_initParam_%sTrials_model%.3d_started%s.mat';
nameFinalLog      = 'finalLog_initParam_%sTrials_model%.3d_started%s.mat';
nameGeneralSAM    = 'SAM_%sTrials.mat';
nameModelSAM      = 'SAM_%sTrials_model%.3d.mat';

for iSubj = subj
  for iModel = model
    
    % Load the general SAM file
    % =======================================================================================================================
    load(fullfile(sprintf(rootDir,iSubj,dt,trialVarStr,modelArch),sprintf(nameGeneralSAM,optimScope)),'SAM');

    % Specify model-specific SAM structure
    % =======================================================================================================================
    SAM             = sam_spec_job_specific(SAM,iModel);
    
    % Load the set of optimized initial parameters
    % =======================================================================================================================
    finalLogs       = regexpdir(sprintf(rootDir,iSubj,dt,trialVarStr,modelArch),sprintf(nameFinalLog,optimScope,iModel,'.*'));
    iterLogs        = regexpdir(sprintf(rootDir,iSubj,dt,trialVarStr,modelArch),sprintf(nameInitLog,optimScope,iModel,'.*'));

    if ~isempty(finalLogs)
      load(finalLogs{end},'bestX0')
      SAM.optim.x0  = bestX0(1:SAM.optim.nStartPoint,3:end);

    elseif ~isempty(iterLogs)
      load(iterLogs{end},'bestX0')
      SAM.optim.x0  = bestX0(1:SAM.optim.nStartPoint,3:end);
    end

    % Save the model-specific SAM file
    % =======================================================================================================================
    save(fullfile(sprintf(rootDir,iSubj,dt,trialVarStr,modelArch),sprintf(nameModelSAM,optimScope,iModel)),'SAM');

    keep iSubj subj dt trialVarStr modelArch optimScope iModel model rootDir nameInitLog nameFinalLog nameGeneralSAM nameModelSAM

  end
end
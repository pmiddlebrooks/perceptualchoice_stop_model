function job_merge_cost_dist(subj,dt,trialVar,optimScope,choiceMech,stopMech,model,fileStr)
% ACCRE_JOB_MERGE_FVALANDBESTX <Synopsis of what this function does> 
%  
% DESCRIPTION 
% <Describe more extensively what this function does> 
%  
% SYNTAX 
% ACCRE_JOB_MERGE_FVALANDBESTX(subj,dt,trialVar,optimScope,choiceMech,stopMech,model); 
%  
% EXAMPLES 
%  
%  
% REFERENCES 
%  
% ......................................................................... 
% Bram Zandbelt, bramzandbelt@gmail.com 
% $Created : Wed 16 Apr 2014 08:52:13 CDT by bram 
% $Modified: Wed 16 Apr 2014 08:52:13 CDT by bram 

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

% Loop over subjects and models
for iSubj = subj
  for iModel = model
    
    
    % Identify all final log files for this subject and model
    allFinalLog = regexpdir(sprintf(fileStr.root,iSubj,dt,trialVarStr,modelArch),sprintf('finalLog_cost_%sTrials_model%.3d_startSet.*_started.*.mat$',lower(optimScope),iModel));
    nFile       = size(allFinalLog,1);
    
    
chi2 = [];    
BIC = [];    
    for iFile = 1:nFile
       % Load file
      load(allFinalLog{iFile},'history');
     
      chi2 = [chi2; history(:,1)];
      BIC = [BIC; history(:,2)];
      
      
    end
      
    % Sort dataset by function value
    ds = sortrows(ds,'FunctionValue');
    
    % Save as tab-delimited text file (with date/time tag of last 
    export(ds,'file',fullfile(sprintf(fileStr.root,iSubj,dt,trialVarStr,modelArch),sprintf(fileStr.output,optimScope,iModel)));
    
  end
end

function hms = sec2hms(t)
    hours = floor(t / 3600);
    t = t - hours * 3600;
    mins = floor(t / 60);
    secs = t - mins * 60;
    hms = sprintf('%02d:%02d:%05.2f', hours, mins, secs);
end
end
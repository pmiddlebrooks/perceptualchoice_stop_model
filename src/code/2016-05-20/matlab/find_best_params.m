function params = find_best_params(subj,trialVar,optimScope,choiceMech,stopMech,dt,model,fileStr)
%% View the best fitting parameters
accreRoot = '/gpfs22';
accreHome = '/home/middlepg';
accreScratch = '/scratch/middlepg';
if isdir(fullfile(accreScratch))
    matRoot = fullfile(accreRoot,accreHome,'m-files');
    modelRoot = fullfile(accreScratch,'perceptualchoice_stop_model');
    environ = 'accre';
else
    matRoot = '/Volumes/HD-1/Users/paulmiddlebrooks/matlab';
    modelRoot = '/Volumes/HD-1/Users/paulmiddlebrooks/perceptualchoice_stop_model';
    environ = 'local';
end

% Define variables
% subj        = 1;
% dt          = 5;
% trialVar    = true;
% optimScope  = 'stop';
% choiceMech  = 'race'; %{'race', 'li', 'ffi'};
% stopMech    = 'none';
% model       = 191;%[79 191 245 478];

% fileStr.root   = fullfile(modelRoot,'data/2016-05-20/preproc01/subj%.2d/dt%d/%s/%s/');
% fileStr.output = 'allFValAndBestX_%sTrials_model%.3d.txt';
% ===============================================================================================


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
for iSubj = 1 : length(subj)
    
    params = table();
    
    for iMod = 1 : length(model)

        % Identify all final log files for this subject and model
        allFinalLog = regexpdir(sprintf(fileStr.fitRoot,subj(iSubj),dt,trialVarStr,modelArch),sprintf('finalLog_%sTrials_model%.3d_startVal.*_started.*.mat$',lower(optimScope),model(iMod)));
        allSAM      = regexpdir(sprintf(fileStr.fitRoot,subj(iSubj),dt,trialVarStr,modelArch),sprintf('SAM_%sTrials_model%.3d_startVal.*_started.*.mat$',lower(optimScope),model(iMod)));
        nFile       = size(allSAM,1);
        
        
        ds          = dataset({[1:nFile]','FileIndex'}, ...
            {cell(nFile,1),'FileNameLog'}, ...
            {cell(nFile,1),'FileNameSAM'}, ...
            {nan(nFile,1),'StartValueIndex'}, ...
            {nan(nFile,1),'FunctionValue'}, ...
            {cell(nFile,1),'ElapsedTime'}, ...
            {false(nFile,1),'ExitFlag'}, ...
            {cell(nFile,1),'BestX'}, ...
            {cell(nFile,1),'StartX'});
        
        for iFile = 1:nFile
            
            % Get and enter final log filename
            [~,fileNameLog]         = fileparts(allFinalLog{iFile});
            ds.FileNameLog{iFile}   = fileNameLog;
            
            % Get and enter final SAM filename
            [~,fileNameSam]         = fileparts(allSAM{iFile});
            ds.FileNameSAM{iFile}   = fileNameSam;
            
            % Load file
            load(allFinalLog{iFile},'fVal','tElapse','exitFlag','X');
            
            % Enter optimization function value
            ds.FunctionValue(iFile) = fVal;
            
            % Enter elapsed time
            hours = floor(tElapse / 3600);
            tElapse = tElapse - hours * 3600;
            mins = floor(tElapse / 60);
            secs = tElapse - mins * 60;
            hms = sprintf('%02d:%02d:%05.2f', hours, mins, secs);
            ds.ElapsedTime{iFile} = hms;
            
            % Enter exit flag
            ds.ExitFlag(iFile) = exitFlag;
            
            % Enter best-fitting parameter values
            ds.BestX{iFile} = X;
            
            % Load SAM and identify and enter starting value
            settings  = regexp(fileNameLog, 'finalLog_(\w*)Trials_model(\w*)_startVal(\w*)_started(.*)', 'tokens');
            iStartVal = str2num(settings{1}{3});
            
            ds.StartValueIndex(iFile) = iStartVal;
            
            % Enter initial parameter values
            load(allSAM{iFile},'SAM');
            ds.StartX{iFile} = SAM.optim.x0(iStartVal,:);
            
        end
        
        % Sort dataset by function value
        iParams = sortrows(ds,'FunctionValue');
        p{iSubj} = iParams.BestX{1};
        
        
        
        
        
        % Create a table
        paramNames = SAM.model.XCat.name;
        
        
        iParam = table();
        iParam.subject = subj(iSubj);
        iParam.model = model(iMod);
        iParam.scope = {'go'};
        for j = 1 : length(paramNames)
            iParam.(paramNames{j}) = {p{iSubj}(SAM.model.variants.tree(model(iMod)).XSpec.i.all.iCatClass{1,j})};
        end
        
        params = [params; iParam];
        
        iParam = table();
        iParam.subject = subj(iSubj);
        iParam.scope = {'stop'};
        iParam.model = model(iMod);
        for j = 1 : length(paramNames)
            iParam.(paramNames{j}) = {p{iSubj}(SAM.model.variants.tree(model(iMod)).XSpec.i.all.iCatClass{2,j})};
        end
        params = [params; iParam];
        
    end  % Model for loop
    
    
    saveDir             = fullfile(sprintf(fileStr.result,subj(iSubj),dt,trialVarStr,modelArch));
    writetable(params,fullfile(saveDir,'bestFitParams.csv'))
end

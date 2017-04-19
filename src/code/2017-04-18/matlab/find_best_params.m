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

% fileStr.root   = fullfile(modelRoot,'data/2016-08-18/preproc01/subj%.2d/dt%d/%s/%s/');
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

pramSubj = cell(length(subj), 1);
    params = table();
    paramsStart = table();

% Loop over subjects and models
for iSubj = 1 : length(subj)
    
    
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
        allParams = sortrows(ds,'FunctionValue');
        paramsBest{iSubj} = allParams.BestX{1};
        paramsInitial{iSubj} = allParams.StartX{1};
        paramNames = SAM.model.XCat.name;
        
        
        
        
        % Figure out if we're getting go, stop, or all parameters (if the
        % table will have 1 or two rows for each model)
        switch optimScope
            case 'go'
                nRow = 1;
            case {'stop','all'}
                nRow = 2;
        end
        
        iParam = table();
        iParamStart = table();
        
        for k = 1 : nRow
            kParam = table();
            kParamStart = table();
            
            % Define a few descriptive table variables
            kParam.subject =  subj(iSubj);
            kParamStart.subject =  subj(iSubj);
            kParam.model =  model(iMod);
            kParamStart.model =  model(iMod);
            if k == 1
                kParam.unit = {'go'};
                kParamStart.unit = {'go'};
            elseif k == 2
                kParam.unit = {'stop'};
                kParamStart.unit = {'stop'};
            end
            
            % Loop through the parameters and fill in table.
            for j = 1 : length(paramNames)
                jParamInd = SAM.model.variants.tree(model(iMod)).XSpec.i.(optimScope).iCatClass{k,j};
                jParam = paramsBest{iSubj}(jParamInd);
                jParamStart = paramsInitial{iSubj}(jParamInd);
                if isempty(jParam)
                    jParam = nan;
                    jParamStart = nan;
                end
                kParam.(paramNames{j}) = {jParam};
                kParamStart.(paramNames{j}) = {jParamStart};
            end
            
            % Add parameters and desciptive variables to the table
            iParam = [iParam; kParam];
            iParamStart = [iParamStart; kParamStart];
        end
        
        % Add parameter row to full table
        params = [params; iParam];
        paramsStart = [paramsStart; iParamStart];
        
        
    end  % Model for loop
    
    saveDir             = fullfile(sprintf(fileStr.result,subj(iSubj),dt,trialVarStr,modelArch));
    if exist(saveDir,'dir') ~= 7
        mkdir(saveDir)
    end
    writetable(params,fullfile(saveDir,[optimScope,'_bestFitParams.csv']))
    writetable(paramsStart,fullfile(saveDir,[optimScope,'_startParams.csv']))
    
    
end

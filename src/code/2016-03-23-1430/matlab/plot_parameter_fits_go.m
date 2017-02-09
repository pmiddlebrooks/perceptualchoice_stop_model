function plot_parameter_fits_go(subject,model,architecture,dt,trialVar,fileStr,savePlot)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. PROCESS INPUTS & SPECIFY VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

% 1.1. Process inputs
% =========================================================================

if trialVar
    trialVarStr = 'trialvar';
else
    trialVarStr = 'notrialvar';
end

% 1.2. Specify dynamic variables
% =========================================================================

% Numbers
nSubject            = numel(subject);
nModel              = numel(model);
nArchitecture       = numel(architecture);

% 1.3. Specify static variables
% =========================================================================
optimScope          = 'go';
rootDir             = fileStr.root;
nameFVal            = 'allFValAndBestX_%sTrials_model%.3d.txt';
modelStr            = {'zc','t0','v/ve','zc & t0','zc & v/ve','t0 & v/ve','zc, t0 & v/ve'};

nParamSets          = 4;  % How many sets of parameters do we want to visualize?
nRow                = sqrt(nParamSets);
nCol                = sqrt(nParamSets);
sets                = fullfact([nRow nCol]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. Plot the data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop over subjects, architectures, and models
warning('off','stats:dataset:ModifiedVarnames')
for iSubject = 1:nSubject
    for iArchitecture = 1:nArchitecture
        figureHandle = 39 + iSubject * iArchitecture;
        
        % Set up the figure and panels
        %         [axisWidth, axisHeight, xAxesPosition, yAxesPosition] = standard_figure;
        standard_figure(nRow,nCol,'landscape', figureHandle);
        
        %                 set_figure({1024,574,'pixels'},{'USLetter','landscape'},{'Helvetica',18});
        p = panel();
        p.pack(num2cell(repmat(1/nRow,1,nRow)), num2cell(repmat(1/nCol,1,nCol)));
        
        annotation('textbox', [0 0.9 1 0.1], ...
            'String', sprintf('Subject %d, architecture %s',subject(iSubject),architecture{iArchitecture}), ...
            'Interpreter','none', ...
            'EdgeColor', 'none', ...
            'HorizontalAlignment', 'center');
        
        for iModel = 1:nModel
            
            % Display progress
            fprintf(1,'Working on subject %d, architecture %s, and model %d. \n',subject(iSubject),architecture{iArchitecture},model(iModel));
            
            % Load model-specific fits with cost function values
            ds = dataset('File',fullfile(sprintf(rootDir,subject(iSubject),dt,trialVarStr,architecture{iArchitecture}),sprintf(nameFVal,optimScope,model(iModel))));
            
    % Identify all final log files for this subject and model
    allFinalLog = regexpdir(sprintf(fileStr.root,iSubject,dt,trialVarStr,architecture{iArchitecture}),sprintf('finalLog_%sTrials_model%.3d_startVal.*_started.*.mat$',lower(optimScope),model(iModel)));
    allSAM      = regexpdir(sprintf(fileStr.root,iSubject,dt,trialVarStr,architecture{iArchitecture}),sprintf('SAM_%sTrials_model%.3d_startVal.*_started.*.mat$',lower(optimScope),model(iModel)));
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
    ds = sortrows(ds,'FunctionValue');
            
            
    
    
    
    
    
            bestX = cell2num(ds.BestX);
            startX = cell2num(ds.StartX);
            
            
            bestClr= [1 0 0];
            startClr = [0 1 0];
            lineColor = [0 0 0];
            bestObs = 'o';
            markerSize = 50;
            
            % Plot it
            for iSet = 1 : nParamSets

                iRow = sets(iSet, 1);
                iCol = sets(iSet, 2);
            p(iCol,iRow).select();
            p(iCol,iRow).hold('on');
%             p(iCol,iRow).title({sprintf('Model %d',model(iModel)), ...
%                 sprintf('\\chi^2 = %.1f',cost), ...
%                 sprintf('BIC = %.1f',altCost)});
%                 
%                         set(gca,'XLabel',[200 700], ...
%                             'XTick',100:100:700, ...
%                             'YLim',[0 1], ...
%                             'YTick',0:0.2:1)
            iParam1 = iSet * 2 - 1;
            iParam2 = iSet * 2;
            
            scatter(bestX(:,iParam1), bestX(:,iParam2), markerSize, bestClr, 'filled')
            scatter(startX(:,iParam1), startX(:,iParam2), markerSize, startClr, 'filled')
            for iLine = 1 : size(startX,1)
            plot([startX(iLine,iParam1); bestX(iLine,iParam1)], [startX(iLine,iParam2); bestX(iLine,iParam2)], 'color', lineColor)
            end
            
            end
            if savePlot
                saveDir             = fullfile(sprintf(fileStr.result,subject(iSubject),dt,trialVarStr,architecture{iArchitecture}));
                if exist(saveDir,'dir') ~= 7
                    mkdir(saveDir)
                end
                fileName = sprintf('GO_Respond_%s_Accuracy_%s', responseSide, accuracy);
                print(gcf, fullfile(saveDir, fileName),'-dpdf', '-r300')
                print(gcf, fullfile(saveDir, fileName),'-dpng')
            end
        end
    end
end


 end


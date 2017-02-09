function plot_goodness_of_fit(stat,subject,model,architecture,dt,trialVar,optimScope,fileStr,varargin)
%
%
%
%
%
%
%
%
%
%
%
%
%
%   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. PROCESS INPUTS & SPECIFY VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1.1. Process inputs
% =========================================================================

% Any object properties?

if nargin == 0;
elseif nargin == 1
end


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

% Pre-allocate matrices
Y                   = nan(nModel,nArchitecture,nSubject);

% 1.3. Specify static variables
% =========================================================================
rootDir             = fileStr.root;
nameSAMModel        = 'SAM_%sTrials_model%.3d.mat';
nameFVal            = 'allFValAndBestX_%sTrials_model%.3d.txt';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2. EXTRACT RELEVANT DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop over subjects, architectures, and models
warning('off','stats:dataset:ModifiedVarnames')
for iSubject = 1:nSubject
    for iArchitecture = 1:nArchitecture
        for iModel = 1:nModel
            
            % Display progress
            disp(sprintf('Working on subject %d, architecture %s, and model %d.',subject(iSubject),architecture{iArchitecture},model(iModel)));

            
            file = fullfile(sprintf(rootDir,subject(iSubject),dt,trialVarStr,architecture{iArchitecture}),sprintf(nameFVal,optimScope,model(iModel)));
            
            if exist(file,'file')
            
                % Load model-specific SAM
                ds = dataset('File',file);

                % Load SAM
                load(fullfile(sprintf(rootDir,subject(iSubject),dt,trialVarStr,architecture{iArchitecture}),ds.FileNameSAM{1}));

                switch lower(stat)
                    case 'chisquare'
                        switch SAM.optim.cost.stat.stat
                            case 'chisquare'
                                y = SAM.estim.fVal;
                            case 'bic'
                                y = SAM.estim.fValAlt;
                        end
                    case {'bic','dbic','bayesfactor','wbic'}

                        switch SAM.optim.cost.stat.stat
                            case 'chisquare'
                                y = SAM.estim.fValAlt;
                            case 'bic'
                                y = SAM.estim.fVal;
                        end
                end

                % Get BIC value of best-fitting model
                Y(iModel,iArchitecture,iSubject) = y;
            
            else
                
                Y(iModel,iArchitecture,iSubject) = NaN;
                
            end
        end
    end
end
clear y
warning on

assignin('base','Y',Y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 3. PLOT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up the figure and panels
set_figure({982,663,'pixels'},{'USLetter','landscape'},{'Helvetica',18});

% Plot stuff
p = panel();
p.margin = [20 20 10 10];
p.pack(num2cell([.05 .95/nSubject.*ones(1,nSubject)]), 1);
        
if nSubject == 1
annotation('textbox', [0 0.95 1 0.05], ...
           'String', sprintf('%s: Subject %d',stat,subject(iSubject)), ...
           'Interpreter','none', ...
           'EdgeColor', 'none', ...
           'HorizontalAlignment', 'center');
end

for iSubject = 1:nSubject
    
    data = Y(:,:,iSubject);
    ranks = nan(size(data));
    [sData,iData] = sort(data(:));
    ranks(iData) = 1:numel(data);
    
    switch lower(stat)
        case 'chisquare'
            y = data;
            yHeight = mat2cell(50*ones(size(data)),nModel,ones(1,nArchitecture));
        case 'bic'
            y = data;
            yHeight = mat2cell(min(data(:))*ones(size(data)),nModel,ones(1,nArchitecture));
        case 'dbic'
            y = data-min(data(:));
            yHeight = mat2cell(1e-1*ones(size(data)),nModel,ones(1,nArchitecture));
        case 'bayesfactor'
            y = exp(-0.5.*(data-min(data(:))));
            yHeight = mat2cell(min(data(:))*ones(size(data)),nModel,ones(1,nArchitecture));
        case 'wbic'
            BICd = data-min(data(:));
            y = exp(-0.5.*BICd)./sum(sum(exp(-0.5.*BICd)));
            yHeight = mat2cell(0.2*ones(size(data)),nModel,ones(1,nArchitecture));
    end
            
    ha = p(iSubject+1,1).select();
    hb = bar(ha,y,1);
    
    textLabel = mat2cell(ranks,nModel,ones(1,nArchitecture));
    
    xLoc = arrayfun(@(x) mean(get(get(x,'Children'),'XData')),hb,'Uni',0);
    
    
    
    
    % General axes settings
    set(ha,'Box','off', ...
           'XLim',[0.5 nModel+0.5], ...
           'XTickLabel',arrayfun(@(in1) sprintf('model %d ',in1),model,'Uni',0));
%     set(ha,'Box','off', ...
%            'XTick',[], ...
%            'XTickLabel',arrayfun(@(in1) sprintf('model %d ',in1),model,'Uni',0));
       
   % Statistic-specific axes settings
    switch lower(stat)
        case 'chisquare'
            
            axStep = 1000;
            yLimMax = ceil(max(y(:))/axStep)*axStep;
            
            set(ha,'Box','off', ...
                   'YScale','linear', ...
                   'YLim',[0 yLimMax], ...
                   'YTick',0:axStep:yLimMax);
        case 'bic'
            axStep = 1000;
            yLimMax = ceil(max(y(:))/axStep)*axStep;
            yLimMin = floor(min(y(:))/axStep)*axStep;
            set(ha,'Box','off', ...
                   'YScale','linear', ...
                   'YLim',[yLimMin yLimMax], ...
                   'YTick',yLimMin:axStep:yLimMax);
        case 'dbic'
            
            axStep = 1000;
            yLimMax = ceil(max(y(:))/axStep)*axStep;
            
            set(hb,'BaseValue',1e-2);
            line(get(ha,'XLim'),[2 2],'Color',[.8 .8 .8],'LineStyle','--','LineWidth',2)
            line(get(ha,'XLim'),[6 6],'Color',[.4 .4 .4],'LineStyle','--','LineWidth',2)
            line(get(ha,'XLim'),[10 10],'Color',[0 0 0],'LineStyle','--','LineWidth',2)
            
            set(ha,'Box','off', ...
                   'YScale','log', ...
                   'YLim',[1e-2 yLimMax], ...
                   'YTick',[1e-2,1e0,1e2,yLimMax]);
        case 'bayesfactor'
            axStep = 1000;
            yLimMax = ceil(max(y(:))/axStep)*axStep;
            
            set(hb,'BaseValue',1e-2);
            line(get(ha,'XLim'),[2 2],'Color',[.8 .8 .8],'LineStyle','--','LineWidth',2)
            line(get(ha,'XLim'),[6 6],'Color',[.4 .4 .4],'LineStyle','--','LineWidth',2)
            line(get(ha,'XLim'),[10 10],'Color',[0 0 0],'LineStyle','--','LineWidth',2)
            
            set(ha,'Box','off', ...
                   'YScale','log', ...
                   'YLim',[1e-2 yLimMax], ...
                   'YTick',[1e-2,1e0,1e2,yLimMax]);
        case 'wbic'
            set(ha,'Box','off', ...
                   'YScale','linear', ...
                   'YLim',[0 1], ...
                   'YTick',0:0.5:1, ...
                   'YTickLabel','0.0|0.5|1.0');
            
    end
    
    yLim = get(gca,'YLim');
    yText = mat2cell((yLim(1) + 0.05*diff(yLim))*ones(size(data)),nModel,ones(1,nArchitecture));
    cellfun(@(x,y,z) text(x,y,num2str(z),'Color','w','HorizontalAlignment','Center','VerticalAlignment','middle','Rotation',90,'FontSize',12),xLoc(:),yText(:),textLabel(:));
    
    ylabel(stat);
end

% BICs summed across subjects
% 
% data = sum(Y,3);
% ranks = nan(size(data));
% [sData,iData] = sort(data(:));
% ranks(iData) = 1:numel(data);
% 
% switch lower(stat)
%     case 'chisquare'
%         y = data;
%     case 'bic'
%         y = data;
%     case 'dbic'
%         y = data-min(data(:));
%     case 'wbic'
%         y = exp(-0.5.*BICd)./sum(sum(exp(-0.5.*BICd)));
% end
% 
% ha = p(nSubject+1).select();
% hb = bar(ha,y,1);
% set(hb,'BaseValue',1e-2);
% 
% line(get(ha,'XLim'),[2 2],'Color',[.8 .8 .8],'LineStyle','--')
% line(get(ha,'XLim'),[6 6],'Color',[.4 .4 .4],'LineStyle','--')
% line(get(ha,'XLim'),[10 10],'Color',[0 0 0],'LineStyle','--')
% 
% textLabel = mat2cell(ranks,nModel,ones(1,nArchitecture));
% 
% xLoc = arrayfun(@(x) mean(get(get(x,'Children'),'XData')),hb,'Uni',0);
% yHeight = mat2cell(1e-1*ones(size(data)),nModel,ones(1,nArchitecture));
% 
% cellfun(@(x,y,z) text(x,y+0.2,num2str(z),'Color','w','HorizontalAlignment','Center','VerticalAlignment','middle','Rotation',90),xLoc(:),yHeight(:),textLabel(:));
% 
% % General axes settings
% set(ha,'Box','off', ...
%        'XTickLabel',arrayfun(@(in1) sprintf('model %d ',in1),model,'Uni',0));
% 
% % Statistic-specific axes settings
% switch lower(stat)
%     case 'chisquare'
%         set(ha,'Box','off', ...
%                'YScale','linear', ...
%                'YLim',[0 1000]);
%     case 'bic'
%         set(ha,'Box','off', ...
%                'YScale','linear', ...
%                'YLim',[24000 30000]);
%     case 'dbic'
%         set(ha,'Box','off', ...
%                'YScale','log', ...
%                'YLim',[1e-2 1e4]);
%     case 'wbic'
%         set(ha,'Box','off', ...
%                'YScale','linear', ...
%                'YLim',[0 1]);
% 
% end
% 
% 
% switch lower(stat)
%     case 'dbic'
%         set(findobj('String',' 1'),'Color','k');
%     case 'wbic'
%         
% end
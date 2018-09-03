%%
standard_figure(1,1,'landscape');

%                 set_figure({1024,574,'pixels'},{'USLetter','landscape'},{'Helvetica',18});
p = panel();
% nPlotCol = max(3, nModel);
nPlotCol = 2;
if length(conditionArray) == 2
    p.pack({.05, .45, .45, .05}, num2cell(repmat(1/nPlotCol,1,nPlotCol)));
elseif length(conditionArray) == 3
    p.pack({.05, .3 .3 .3, .05}, num2cell(repmat(1/nPlotCol,1,nPlotCol)));
elseif length(conditionArray) == 1
    p.pack({.05, .9, .05}, num2cell(repmat(1/nPlotCol,1,nPlotCol)));
end
annotation('textbox', [0 0.9 1 0.1], ...
    'String', sprintf('Subject %d, architecture %s',subject,architecture), ...
    'Interpreter','none', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');




goCCorrMrkObs           = 'none';
goCCorrLnObs            = '-';
goCCorrLnWidthObs          = 2;
goCCorrMrkPrd           = 'none';
goCCorrLnPrd            = '-';
goCCorrLnWidthPrd          = 5;

goCErrorMrkObs          = '^';
goCErrorLnObs           = 'none';
goCErrorMrkPrd          = 'none';
goCErrorLnPrd           = '-.';
goCErrorLnWidth         = 1;

stopIErrorCCorrMrkObs   = 'none';
stopIErrorCCorrLnObs    = '-';
stopIErrorCCorrLnWidthObs  = 2;
stopIErrorCCorrMrkPrd   = 'none';
stopIErrorCCorrLnPrd    = '-';
stopIErrorCCorrLnWidthPrd  = 5;

stopICorrMrkPrd         = 'none';
stopICorrLnPrd          = ':';
stopICorrLnWidth        = 1;



%%
cohList = sort(unique(obsRts.Coherence))';
rspList = sort(unique(obsRts.Response))';
nSsdPlot = 3;

% Initialize table for Jeff
plotCnd = nan(length(cohList) * length(rspList), 2+nSsdPlot);
plotCnd(:, 1:2) = combvec(cohList,rspList)';
ssdListAll = [];
for iCoh = cohList
    for iRsp = 1 : 2
        
        obsStop = obsRts(~isnan(obsRts.Ssd) & obsRts.Coherence == iCoh & obsRts.Response == iRsp, :);
        
        fprintf('\nCoh: %d\tRsp: %d\n', iCoh, iRsp)
        
        ssdList = unique(obsStop.Ssd)'
        [N,edges] = histcounts(obsStop.Ssd,length(ssdList))
        % N
        [nSort, nInd] = sort(N);
        ssdListPlot = sort(ssdList(nInd(end-2:end)));
stopColorList = flipud(linspace(0.45, 0.85, length(ssdListPlot))');
        
        disp(ssdListPlot)
        ssdListAll = [ssdListAll; ssdListPlot];
        
        % Create a table for excel for Jeff
        
        
        % Set up the plots
        p(iCoh+1, iRsp).select();
        p(iCoh+1, iRsp).title({sprintf('Resp: %d, Coh: %d', iRsp, iCoh)});
        p(iCoh+1, iRsp).hold('on');
        
        % Plot the Go RTs
        % ------------------------
        
        % Get the Go trials
        obsGoRt = sort(obsRts.RT(isnan(obsRts.Ssd) & obsRts.Coherence == iCoh & obsRts.Response == iRsp, :));
        cumPGoCCorrObs   = cmtb_edf(obsGoRt,obsGoRt);
        
        prdGoRt = sort(prdRts.RT(isnan(prdRts.Ssd) & prdRts.Coherence == iCoh & prdRts.Response == iRsp, :));
        cumPGoCCorrPrd   = cmtb_edf(prdGoRt,prdGoRt);

        
        % Plot the Go trials
                        iColor = 'k';
                        plot(obsGoRt,cumPGoCCorrObs,'Color',iColor,'Marker',goCCorrMrkObs,'LineStyle',goCCorrLnObs,'LineWidth',goCCorrLnWidthObs);
                        plot(prdGoRt,cumPGoCCorrPrd,'Color',iColor,'Marker',goCCorrMrkPrd,'LineStyle',goCCorrLnPrd,'LineWidth',goCCorrLnWidthPrd);
        

        
        % Loop through SSDs and plot them
        for j = 1 : length(ssdListPlot)
            jSsd = ssdListPlot(j);
            
            % Get the Stop trials
        obsStopRt = sort(obsStop.RT(obsStop.Ssd == jSsd & obsStop.Coherence == iCoh & obsStop.Response == iRsp, :));
        cumPStopIErrorCCorrObs   = cmtb_edf(obsStopRt,obsStopRt);
            
        prdStopRt = sort(prdRts.RT(prdRts.Ssd == jSsd & prdRts.Coherence == iCoh & prdRts.Response == iRsp, :));
        cumPStopIErrorCCorrPrd   = cmtb_edf(prdStopRt,prdStopRt);
        
        % Plot the Stop trials
                                iColor = repmat(stopColorList(j),1,3);
                                plot(obsStopRt,cumPStopIErrorCCorrObs,'Color',iColor,'Marker',stopIErrorCCorrMrkObs,'LineStyle',stopIErrorCCorrLnObs,'LineWidth',stopIErrorCCorrLnWidthObs);
                                plot(prdStopRt,cumPStopIErrorCCorrPrd,'Color',iColor,'Marker',stopIErrorCCorrMrkPrd,'LineStyle',stopIErrorCCorrLnPrd,'LineWidth',stopIErrorCCorrLnWidthPrd);

        
        end
        
        
        
        
        % Set axes
        switch subject
            case 1
                set(gca,'XLim',[200 600], ...
                    'XTick',100:100:800, ...
                    'YLim',[0.1 0.9], ...
                    'YTick',0.1:0.2:0.9)
            case 2
                set(gca,'XLim',[200 400], ...
                    'XTick',100:100:500, ...
                    'YLim',[0.1 0.9], ...
                    'YTick',0.1:0.2:0.9)
            case 3
                set(gca,'XLim',[200 600], ...
                    'XTick',100:100:700, ...
                    'YLim',[0.1 0.9], ...
                    'YTick',0.1:0.2:0.9)
            otherwise
                error('Need to add axes limits for subject')
        end
        
        
        
    end
end

plotCnd(:,3:5) = ssdListAll;
plotCnd = array2table(plotCnd, 'VariableNames',  {'Coherence', 'Response', 'SSD1', 'SSD2', 'SSD3'});


%%
    saveDir             = '~/perceptualchoice_stop_model/results/2018-05-01';
    if exist(saveDir,'dir') ~= 7
        mkdir(saveDir)
    end
    fileName = sprintf('%s_%s_Cumulative_rts_SelectSsd_Respond_%s_Coherence_%s_Accuracy_%s', subjectName, addData, responseSide{1}, conditionArray, accuracy{1});
    print(gcf, fullfile(saveDir, fileName),'-dpdf', '-r300')

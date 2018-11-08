function check_gazData(DIR,subID)
load( [ DIR.dataRaw '/subInfo.mat' ])

cond = 'IM';%'CT';

%% Find patients who had only one session
ksub = 0;
subID_grp = unique(subInfo(:,1));
for isub = 1:length(subID_grp)
    state = unique(subInfo(subInfo(:,1)==subID_grp(isub),2));
    if length(state)==1 && state == 1
        continue
    elseif length(state)==1 && state ~= 1
        ksub = ksub + 1;
        oneStateSub(ksub,1) = subID_grp(isub);
        %         oneStateSub(ksub,2) = state;
    end
end

time = linspace(0,60,3600);
%%
for isub = 1: length(subID)
    subNum = str2num(subID{isub}(2:end));
    iplot = 0;

%         if subNum ~= 138
%         continue
%         end
    
%     if ismember(subNum, rejSub)||ismember(subNum,oneStateSub)
%         continue
%     end
    
    load( [ DIR.collectedData '/' subID{isub} ])
    load( [ DIR.cfg '/' subID{isub} '_Cfg_a_BR.mat' ] )
    
    figure(1),clf
                    
    for irun = 1:length(Run)
        
        for iblk = 1:length(Run(irun).block)
            %             gazeData = [];
            
            switch cond
                case 'IM'
                    if ~strcmp( Run(irun).block(iblk).condition, 'mixed')
                        continue
                    end
                    
                case 'CT'
                    if ~strcmp( Run(irun).block(iblk).condition, 'br')||...
                         ~strcmp( Run(irun).block(iblk).condition, 'replay')   
                        continue
                    end
            end
            
            
            if ~isempty(Run(irun).block(iblk).gazeData)
                
                stimDes = Run(irun).block(iblk).mixedLabelDesign;
                
                gazeDataB = Run(irun).block(iblk).gazeData( ...
                    Run(irun).block(iblk).gazeData(:,17) >=Run(irun).block(iblk).trial_evTimes(1) &...
                    Run(irun).block(iblk).gazeData(:,17) <=  Run(irun).block(iblk).trial_evTimes(end), 3) ...
                    *  Cfg.windowRect(3);
                
                stimTime = Run(irun).block(iblk).trial_evTimes - Run(irun).block(iblk).trial_evTimes(1);
                stim = [stimTime, stimDes']; %1) stim onset time, 2) stim offset time, 3) stim
                
                iplot = iplot + 1;
                subplot(length(Run)*4, 1, iplot)
                
                for itri = 1:20
                    switch stim(itri,3)
                        case 1
                            c = 'y';
                        case 2
                            c = 'r';
                        case 3
                            c = 'g';
                    end 
                    
                    if sum(isinf(gazeDataB))~=0
                        tmp = gazeDataB;
                        tmp(isinf(tmp)==1) = []; 
                        patchY = [ min(tmp)*1.05,  min(tmp)*1.05,  max(tmp)*1.05, max(tmp)*1.05 ];
                    else
                        patchY = [ min(gazeDataB)*1.05,  min(gazeDataB)*1.05,  max(gazeDataB)*1.05, max(gazeDataB)*1.05 ];
                    end
                    
                    hh = patch( [stim(itri,1) stim(itri,2) stim(itri,2) stim(itri,1)],...
                        patchY, c );
                    set(hh,'edgeColor', 'none')
                    hold on
                end
                
                timeB = linspace(0,stim(end,2)+1,length(gazeDataB)); 
                adj10s = floor(length(timeB)/60*10);
                
                %%% PLOT OKN %%%-------------------------------------------
                h = plot(timeB,gazeDataB, 'k');
                set(h,'lineWidth',2)
                
                % Adjust xaxis because OKN data are usually longer than
                % 18000 for some reason...
                ax = gca;
                ax.XTick = [ timeB(adj10s*1), timeB(adj10s*2), timeB(adj10s*3), timeB(adj10s*4), timeB(adj10s*5), timeB(adj10s*6) ];
                ax.XTickLabel = {'10','20','30','40','50','60'};
                 
                ylim ([ min(gazeDataB)*1.05, max(gazeDataB)*1.05 ])
                xlim( [0, timeB(adj10s*6)] )
                
                title([ subID{isub} '_r' num2str(irun) '_b' num2str(iblk) ], 'interpret', 'none')
                ylabel('Gaze data (pixel)') 
                xlabel('Time (sec)')
                
                
                

                %%% PLOT BP %%%--------------------------------------------
                
                iplot = iplot + 1;
                subplot(length(Run)*4, 1, iplot)
                
                for itri = 1:20
                    switch stim(itri,3)
                        case 1
                            c = 'y';
                        case 2
                            c = 'r';
                        case 3
                            c = 'g';
                    end
                    hh = patch( [ 3*(itri-1) 3*(itri-1)+2 3*(itri-1)+2 3*(itri-1)], [1.1 1.1 -1.1 -1.1], c );
                    set(hh,'edgeColor', 'none')
                    hold on
                end
                
                switch cond
                    case 'IM'
                        
                        rResp = [];
                        lResp = [];
                        for itri = 1:20
                            rResp = cat(2,rResp,Run(irun).block(iblk).rResp(itri,:));
                            lResp = cat(2,lResp,Run(irun).block(iblk).lResp(itri,:));
                        end
                        altResp = rResp - lResp;
                        
                    case 'CT'
                        altResp = Run(irun).block(iblk).rResp - Run(irun).block(iblk).lResp;
                end
                
                h = plot(time, altResp, 'k');
                set(h,'lineWidth',1)
                ylim ([-1.1 1.1])
                ylabel('L-R Button press')
                xlabel('Time (sec)')
                
                title([ subID{isub} '_r' num2str(irun) '_b' num2str(iblk) ], 'interpret', 'none')

            end
        end
    end
    
    set(figure(1),'PaperUnits','inches','PaperPosition',[0 0 30 2*length(Run)*4])
%     set(figure(1),'PaperUnits','inches','PaperPosition',[0 0 60 5])

    print('-dpng', [ DIR.figRaw '/' cond '_' subID{isub} '_rawTimeCourse.png'])
    
end
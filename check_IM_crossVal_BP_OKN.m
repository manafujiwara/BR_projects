function check_IM_crossVal_BP_OKN(DIR,subjects,Cond,C)

time = linspace(-1, 2, 900);
subID_rej = [ 105, 110, 118, 120, 121, 136, 137, 151]; % reject subjects who have too few trials


if exist([DIR.AUC  'alSub_IM_decodability.mat'], 'file') == 2
    
    load([DIR.AUC  'NC_eachSub_IM_decodability.mat'])
    load([DIR.AUC  'PD_eachSub_IM_decodability.mat'])
else
    
    existFile = dir([DIR.slowPhase 'eyeDataAllSub_RunBlkEpcLabelState*.mat'] );
    if ~isempty(existFile)
        loadData = load(existFile(end).name, 'eye');
    else
        error(['file dose not exist'])
    end
    
    % ------------------------------------------------------------------------
    % eye: 1)subID, 2)irun, 3)iblk, 4)iepc, 5)stimLabel, 6)subjestState,     |
    % 7:906)iOKN data                                                        |
    % ------------------------------------------------------------------------
    
    % ------------------------------------------------------------------------
    % RTinfo 1)subID, 2)irun, 3)iblk, 4)iepc, 5)stimLabel, 6)subjestState,    |
    % 7)lRT, 8)rRT, 9)lRT<rRT, 10)lRT>rRT, 11)l or r, 12)RT 13)button press dur
    % 14:913)choosen response
    % ------------------------------------------------------------------------
    % Note that RT and all related things are discribed in sample points, not
    % time (sec or msec)
    
    eyeData = loadData.eye;
    % eyeinfo = [ eye_amb( ismember(eye_amb(:,1:6), RTinfo(:,1:6), 'rows'), 1:6),...
    %     RTinfo(:,[ 11, 16:17 ]), eye_amb( ismember(eye_amb(:,1:6), RTinfo(:,1:6), 'rows'), 7:906 )];
    % ------------------------------------------------------------------------
    % eyeinfo 1)subID, 2)irun, 3)iblk, 4)iepc, 5)stimLabel, 6)subjestState,   |
    % 7)RT, 8)choosenLengthPress, 9)l or r label by resp(l=-1, r=1),         |
    % 10)iOKN data                                                           |
    % ------------------------------------------------------------------------
    % Note that RT and all related things are discribed in sample points, not
    % time (sec or msec)
    
    subID_1     = unique( eyeData( eyeData(:,6) == 1, 1 ) );
    subID_2     = unique( eyeData( eyeData(:,6) == 2, 1 ) );
    subID_3     = unique( eyeData( eyeData(:,6) == 3, 1 ) );
    subID_4     = unique( eyeData( eyeData(:,6) == 4, 1 ) );
    subID_5     = unique( eyeData( eyeData(:,6) == 5, 1 ) );
    subID_com23 = subID_2( ismember(subID_2, subID_3) );
    subID_com45 = subID_4( ismember(subID_4, subID_5) );
    subID_PD = sort([subID_com23; subID_com45]);
    
    subID = [subID_1;subID_PD];
    vsub = 0;
    
    for isub = 1:length(subID)
        
        if ismember( subID(isub), subID_rej )
            continue
        end
        
        vsub = vsub + 1;
        
        if ismember( subID(isub), subID_1 )
            state = 1;
            stateName = 'NC';
            
        elseif ismember( subID(isub), subID_PD)
            state = 2:5;
            stateName = 'PD';
            
        else
            error('Subject state is unknown')
        end
        
        for istim = 1:2 % ambiguous or unambiguous
            
            switch istim
                
                case 1 %unambiguous
                    load( [DIR.epoched 'IM_RTinfo_unambLabel.mat'], 'RTinfo')
                    BPinfo = RTinfo( RTinfo(:,1) == subID(isub) & ismember( RTinfo(:,5), [2,3] ) & ismember( RTinfo(:,6), state ), [1:6, 11:913] );
                    eyeInfo = eyeData( eyeData(:,1) == subID(isub) & ismember( eyeData(:,5), [2,3] ), : ); %pull only unambiguous trials
                    stimName = 'unamb';
                    C.task = 'unamb';
                    
                case 2 % ambiguous
                    load( [DIR.epoched 'IM_RTinfo_AmbLabel02.mat'], 'RTinfo')
                    BPinfo = RTinfo( RTinfo(:,1) == subID(isub) & RTinfo(:,5) == 1 & ismember( RTinfo(:,6), state ), [1:6, 11:913] );
                    eyeInfo = eyeData( eyeData(:,1) == subID(isub) & eyeData(:,5) == 1 & ismember( eyeData(:,6), state ), : ); %pull only ambiguous trials
                    stimName = 'br';
                    C.task = 'br';
                    
            end
            
            % ------------------------------------------------------------------------
            % RTinfo 1)subID, 2)irun, 3)iblk, 4)iepc, 5)stimLabel, 6)subjestState,    |
            % 7)lRT, 8)rRT, 9)lRT<rRT, 10)lRT>rRT, 11)l or r, 12)RT 13)button press dur
            % 14:913)choosen response
            % ------------------------------------------------------------------------
            % Note that RT and all related things are discribed in sample points, not
            % time (sec or msec)
            
            % ------------------------------------------------------------------------
            % BPinfo 1)subID, 2)irun, 3)iblk, 4)iepc, 5)stimLabel, 6)subjestState,    |
            % 7)l or r, 8)RT, 9)button press dur, 10:909)choosen response
            % ------------------------------------------------------------------------
            % Note that RT and all related things are discribed in sample points, not
            % time (sec or msec)
            
            lBPinfo = BPinfo( BPinfo(:,7) == -1, :);
            rBPinfo = BPinfo( BPinfo(:,7) == 1, :);
            
            for i = 1:size(lBPinfo,1)
                lOKN(i,:) = [ lBPinfo(i,1:9), eyeInfo( ismember( eyeInfo(:, 1:6), lBPinfo(i, 1:6), 'rows' ), 7:906 ) ];
            end
            
            for i = 1:size(rBPinfo,1)
                rOKN(i,:) = [ rBPinfo(i,1:9), eyeInfo( ismember( eyeInfo(:, 1:6), rBPinfo(i, 1:6), 'rows' ), 7:906 ) ];
            end
            
            %         labels(1:size(lOKN,1)) = {'l'};
            %         resp(1:size(lOKN,1))   = 0;
            %         labels(size(lOKN,1)+1 : size(lOKN,1)+size(rOKN,1)) = {'r'};
            %         resp(size(lOKN,1)+1 : size(lOKN,1)+size(rOKN,1))   = 1;
            
            lOKN_4dec = lOKN(:,10:909);
            rOKN_4dec = rOKN(:,10:909);
            
            lOKN_4dec(isnan(lOKN_4dec)) = 0;
            rOKN_4dec(isnan(rOKN_4dec)) = 0;
            
            cResp = cell(2,1); %- trial types x 1
            
            %% cross validation, decoding
            for itp = 1:900
                
                cResp{1} = [lOKN_4dec(:,itp)];
                cResp{2} = [rOKN_4dec(:,itp)];
                nTrialPerClass(1) = size(cResp{1},1);
                nTrialPerClass(2) = size(cResp{2},1);
                
                %%
                if min(nTrialPerClass) < 2 ;
                    disp(['# of trials are less than 2 for some class. skip.'])
                    isNotEnoughTrials = 1;
                    continue
                end
                
                C.nFoldValidation = min(nTrialPerClass);
                
                [ind_train,ind_test] = leaveOneOutRLSC_mem(cResp,C);
                [P{itp}] = evaluatePerformance_mem(cResp,C,ind_train,ind_test);
                
                if mod (itp, 300) == 0
                    disp([ 'sub = ' stateName num2str(subID(isub)) ': stim='...
                        num2str(istim) ': time = ' num2str(itp) '/' num2str(900) ' : ' datestr(now) ...
                        ': # of trials = ' num2str(nTrialPerClass)])
                end
                
            end
            
            
            %%
            for itp = 1:900
                if itp == 1
                    decodability(vsub, itp: itp+3, istim)   = [ lOKN(1, [1,5,6]), P{itp}.zmcorrTest ];
                    decodabilitySD(vsub, itp: itp+3, istim) = [ lOKN(1, [1,5,6]), P{itp}.stdcorrTest ];
                    decodeTrain(vsub, itp: itp+3, istim)    = [ lOKN(1, [1,5,6]), P{itp}.zmcorrTrain ];
                    decodeTrainSD(vsub, itp: itp+3, istim)  = [ lOKN(1, [1,5,6]), P{itp}.stdcorrTrain ];
                else
                    decodability(vsub,itp+3,istim)   = P{itp}.zmcorrTest;
                    decodabilitySD(vsub,itp+3,istim) = P{itp}.stdcorrTest;
                    decodeTrain(vsub,itp+3,istim)    = P{itp}.zmcorrTrain;
                    decodeTrainSD(vsub,itp+3,istim)  = P{itp}.stdcorrTrain;
                end
            end
            
            
            mlOKN = nanmean(lOKN(:,10:909),1);
            mrOKN = nanmean(rOKN(:,10:909),1);
            
            lOKN_alSub(vsub,:,istim) = [ lOKN(1, [1,5,6]), mlOKN ]; % - 1) subID, 2)stimLabel, 3)subjectState, 4:903)OKN
            rOKN_alSub(vsub,:,istim) =  [ rOKN(1, [1,5,6]), mrOKN ];
            
            %% Plot eye data and Decodability in each subject
            figure(2)
            
            switch stateName
                case 'NC'
                    plotColor = {'b', 'r'};
                case 'PD'
                    plotColor = {'g', 'm'};
            end
            
            switch istim
                case 1 %Unambiguous
                    subplot(2,2,1)
                    plot(time,lOKN(:,10:909),'b')
                    hold on
                    plot(time,rOKN(:,10:909),'r')
                    xlabel('time from stimulus onset (sec)', 'FontSize', 16)
                    ylabel('velocity (deg/sec)')
                    set(gca, 'FontSize', 16)
                    title([ 'lOKN_and_rOKN_Unamb'], 'interpret', 'none', 'FontSize', 16)
                case 2 %BR
                    subplot(2,2,3)
                    plot(time,lOKN(:,10:909),'b')
                    hold on
                    plot(time,rOKN(:,10:909),'r')
                    xlabel('time from stimulus onset (sec)', 'FontSize', 16)
                    ylabel('velocity (deg/sec)')
                    set(gca, 'FontSize', 16)
                    title([ 'lOKN_and_rOKN_BR'], 'interpret', 'none', 'FontSize', 16)
            end
            plot([0 0],ylim, 'k--')
            plot([-1 2], [.5 .5],'k')
            
            clear rOKN lOKN
        end
        
        subplot(2,2,4),
        h1 = shadedErrorBar(time, decodability(vsub,4:903,1), decodabilitySD(vsub,4:903,1), 'c', 1);
        hold on
        h2 = shadedErrorBar(time, decodability(vsub,4:903,2), decodabilitySD(vsub,4:903,2), 'm', 1);
        legend([h1.mainLine, h2.mainLine],'Unamb', 'BR', 'Location', 'SouthEast')
        legend boxoff
        ylim([0, 1.2])
        plot([0 0],ylim, 'k--')
        plot([-1 2], [.5 .5],'k')
        xlabel('time from stimulus onset (sec)', 'FontSize', 16)
        ylabel('AUC', 'FontSize', 16)
        set(gca,'FontSize', 16)
        title('Decodability_Unamb_and_BR', 'interpret', 'none')
        
        suptitle_m([ stateName num2str( subID(isub) )])
        
        print('-dpng', [DIR.figDecode stateName num2str( subID(isub) ) '_decodability_eachTP'])
        figure(2), clf
        
    end 
    
    save([DIR.AUC  'alSub_IM_decodability.mat'], 'decodability', 'decodabilitySD', 'decodeTrain', 'decodeTrainSD', 'C')
    
    decodability_NC(:,:,1)   = decodability(decodability(:,3,1) == 1, :, 1);
    decodability_NC(:,:,2)   = decodability(decodability(:,3,2) == 1, :, 2);
    decodabilitySD_NC(:,:,1) = decodabilitySD(decodabilitySD(:,3,1) == 1, :, 1);
    decodabilitySD_NC(:,:,2) = decodabilitySD(decodabilitySD(:,3,2) == 1, :, 2);
    decodeTrain_NC(:,:,1)    = decodeTrain(decodeTrain(:,3,1) == 1, :, 1);
    decodeTrain_NC(:,:,2)    = decodeTrain(decodeTrain(:,3,2) == 1, :, 2);
    decodeTrainSD_NC(:,:,1)  = decodeTrainSD(decodeTrainSD(:,3,1) == 1, :, 1);
    decodeTrainSD_NC(:,:,2)  = decodeTrainSD(decodeTrainSD(:,3,2) == 1, :, 2);
    
    decodability_PD(:,:,1)   = decodability( ismember( decodability(:,3,1), 2:5 ), :, 1);
    decodability_PD(:,:,2)   = decodability( ismember( decodability(:,3,2), 2:5 ), :, 2);
    decodabilitySD_PD(:,:,1) = decodabilitySD( ismember( decodabilitySD(:,3,1), 2:5 ), :, 1);
    decodabilitySD_PD(:,:,2) = decodabilitySD( ismember( decodabilitySD(:,3,2), 2:5 ), :, 2);
    decodeTrain_PD(:,:,1)    = decodeTrain( ismember( decodeTrain(:,3,1), 2:5 ), :, 1);
    decodeTrain_PD(:,:,2)    = decodeTrain( ismember( decodeTrain(:,3,2), 2:5 ), :, 2);
    decodeTrainSD_PD(:,:,1)  = decodeTrainSD( ismember( decodeTrainSD(:,3,1), 2:5 ), :, 1);
    decodeTrainSD_PD(:,:,2)  = decodeTrainSD( ismember( decodeTrainSD(:,3,2), 2:5 ), :, 2);
    
    save([DIR.AUC  'NC_eachSub_IM_decodability.mat'], 'decodability_NC', 'decodabilitySD_NC', 'decodeTrain_NC', 'decodeTrainSD_NC', 'C')
    save([DIR.AUC  'PD_eachSub_IM_decodability.mat'], 'decodability_PD', 'decodabilitySD_PD', 'decodeTrain_PD', 'decodeTrainSD_PD', 'C')
    
end
%% plot decodability

m_decodability_NC   = squeeze( mean(decodability_NC(:,4:903,:), 1) );
sem_decodability_NC = squeeze( std(decodability_NC(:,4:903,:), 1) ) /sqrt( size(decodability_NC,1) );
m_decodability_PD   = squeeze( mean(decodability_PD(:,4:903,:), 1) );
sem_decodability_PD = squeeze( std(decodability_PD(:,4:903,:), 1) ) /sqrt( size(decodability_PD,1) );

for icond = 1:4
    
    switch icond
        case 1 % NC-Unamb
            data = squeeze(decodability_NC(:,4:903,1));
        case 2 % NC-BR
            data = squeeze(decodability_NC(:,4:903,2));
        case 3 % PD-Unamb
            data = squeeze(decodability_PD(:,4:903,1));
        case 4 % PD-BR
            data = squeeze(decodability_PD(:,4:903,2));
    end
    
    [h01(icond,:),p01(icond,:)] = ttest(data,.5, 'Alpha', .01, 'tail', 'right');
    [h001(icond,:),p001(icond,:)] = ttest(data,.5, 'Alpha', .001, 'tail', 'right');
    
end

h01(h01==0) = nan;
h001(h001==0) = nan;


figure(3),clf
subplot(3,1,[1:2])

h(1) = shadedErrorBar(time, m_decodability_NC(:,1), sem_decodability_NC(:,1), 'b', 1);
hold on
h(2) = shadedErrorBar(time, m_decodability_NC(:,2), sem_decodability_NC(:,2), 'r', 1);
h(3) = shadedErrorBar(time, m_decodability_PD(:,1), sem_decodability_PD(:,1), 'c', 1);
h(4) = shadedErrorBar(time, m_decodability_PD(:,2), sem_decodability_PD(:,2), 'm', 1);

for i = 1:4
    set(h(i).mainLine, 'LineWidth', 2)
end

legend([ h(1).mainLine, h(2).mainLine, h(3).mainLine, h(4).mainLine ],...
    'NC-Unamb','NC-BR','PD-Unamb','PD-BR', 'Location', 'NorthWest')
legend boxoff


%%% significant poit
plot(time, h01(1,:)*.45, '.b', 'MarkerSize', 5)
plot(time, h01(2,:)*.43, '.r', 'MarkerSize', 5)
plot(time, h01(3,:)*.41, '.c', 'MarkerSize', 5)
plot(time, h01(4,:)*.39, '.m', 'MarkerSize', 5)

plot(time, h001(1,:)*.45, '.b', 'MarkerSize', 20)
plot(time, h001(2,:)*.43, '.r', 'MarkerSize', 20)
plot(time, h001(3,:)*.41, '.c', 'MarkerSize', 20)
plot(time, h001(4,:)*.39, '.m', 'MarkerSize', 20)

plot([0 0],ylim, 'k--')
plot([-1 2], [.5 .5],'k')

set(gca, 'FontSize', 16)
xlabel('time from stimulus onset (sec)  ', 'FontSize', 16)
ylabel('decodability ', 'FontSize', 16)

title('mean_decodability_accSub_eachSubGroup_eachCond  ', 'interpret', 'none')
set(gcf,'PaperUnits','centimeters','PaperPosition', [0 0 20 15])
print('-dpng', [DIR.figDecode  'alSub_mdecodability_accSub_eachState_eachCond'])


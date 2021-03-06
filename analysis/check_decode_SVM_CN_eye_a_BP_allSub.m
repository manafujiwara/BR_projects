function check_decode_SVM_CN_eye_a_BP_allSub(DIR, subID)
rejCriOKN = 'NanZero';
%1) NaN or 2)NanZero; Reject episodes that contain NaN (or NaN and zero)
%more than hald of it %%%
F.subRej = 0;
F.median = 1;
if F.median
    filtWin = 30;
end
F.biResp = 1;
F.cleanCriteria = 'All'; % 1) OKN or 2)All (both OKN and button press resp)
timeWin  = 500; % ms

switch F.cleanCriteria
    case 'OKN'
        clnName = 'rejEpByOKN';
    case 'All'
        clnName = 'rejEpByAll';
end

if F.median
    load([ DIR.allSub '/CN_subDecodePerCorrectness_m10crossVal_' clnName rejCriOKN '_biResp_timeWin' num2str(timeWin) ...
        '_medFilt' num2str(filtWin) 'pnt.mat' ], 'mperCorrectnessStore_resp', 'mperCorrectnessStore_borderPos', 'perCorrectness'  )
else
    load([ DIR.allSub '/CN_subDecodePerCorrectness_m10crossVal_' clnName rejCriOKN '_biResp_timeWin' num2str(timeWin)...
        '.mat' ], 'mperCorrectnessStore_resp', 'mperCorrectnessStore_borderPos', 'perCorrectness' )
end

%% Check OKN of subject who had low decoding accuracy %%%
stimName = {[], [], [], 'Replay'};
stateName = {'C','Med-On', 'Med-Off', 'DBS-On', 'DBS-Off'};
shresholdAcc = 75; % shreshold for high and low
binVel = -15:15;

for iqualAcc = 1:2 % high and low
    
    rawOKNStore = [];
    medOKNStore = [];
    
    switch iqualAcc
        case 1 % group of subject who had LOW decoding accuracy
            qualAcc = 'Low';
            accAllSubInfo = mperCorrectnessStore_resp( mperCorrectnessStore_resp(:,4) < shresholdAcc, : );
            
        case 2 % group of subject who had HIGH decoding accuracy
            qualAcc = 'High';
            accAllSubInfo = mperCorrectnessStore_resp( mperCorrectnessStore_resp(:,4) > shresholdAcc, : );
    end
    
    accAllSubID = unique( accAllSubInfo(:,1) );
    jsub = 0;
    
    for isub = 1:length(accAllSubID)
        subName = num2str(accAllSubID(isub));
        %     if F.median
        %         D = dir( [ DIR.clnData '*' num2str(lowAccAllSubID(isub)) '_60Hz_medFilt' num2str(filtWin) '_cln*.mat']);
        %         load( [DIR.clnData D(2).name ],...
        %             'clnOKN_OKN', 'clnResp_OKN', 'clnBiResp_OKN', 'clnBorderPos_OKN')
        %         load( [DIR.clnData D(1).name ],...
        %             'clnOKN_all', 'clnResp_all', 'clnBiResp_all', 'clnBorderPos_all')
        %     else
        %         D = dir( [ DIR.clnData '*' num2str(lowAccAllSubID(isub)) '_60Hz_cln*.mat']);
        %         load( [DIR.clnData D(2).name ],...
        %             'clnOKN_OKN', 'clnResp_OKN', 'clnBiResp_OKN', 'clnBorderPos_OKN')
        %         load( [DIR.clnData D(1).name ],...
        %             'clnOKN_all', 'clnResp_all', 'clnBiResp_all', 'clnBorderPos_all')
        %     end
        
        D = dir( [ DIR.excel '/*' num2str(accAllSubID(isub)) '_60Hz.mat' ] );
        load( [ DIR.excel '/' D.name ], 'd_CT_mat' )
        % d_CT_mat: 1)subNum, 2)state, 3)time, 4)stimCondition, 5)alternative resp(just subtracted), 6)smothed OKN vel, ...
        % 7)borderPos_lb , 8:9)left/right pupil diameter ] );
        
        accSubInfo = accAllSubInfo( accAllSubInfo(:,1)==accAllSubID(isub), : );
        % 1)subNum, 2)state, 3)stimCond, 4)decoding accuracy
        
        state    = unique( accSubInfo(:,2) );
        stimCond = unique( accSubInfo(:,3) );
        
        for ises = 1:size(accSubInfo,1)
            
            state = accSubInfo(ises,2);
            stimCond = accSubInfo(ises,3);
            
            if state == 1
                stateCat2 = 'c';
            else
                stateCat2 = 'p';
            end
            
            
            jsub = jsub + 1;
            
            rawOKN = d_CT_mat( d_CT_mat(:,2) == state & d_CT_mat(:,4) == stimCond, 6 );
            medOKN = medfilt1(rawOKN,filtWin);
            
            kMedOKN = kurtosis(medOKN);
            
            perNanRaw  = sum( isnan(rawOKN) ) / length(rawOKN) * 100;
            perZeroRaw = sum( rawOKN == 0 ) / length(rawOKN) * 100;
            
            
            accAllSubInfo( accAllSubInfo(:,1)==accAllSubID(isub) &...
                accAllSubInfo(:,2)==state &...
                accAllSubInfo(:,3)==stimCond, 5:7 ) = [ kMedOKN, perNanRaw, perZeroRaw ];
            
            
            histRawOKN = hist( rawOKN, binVel);
            histMedOKN = hist( medOKN, binVel);
            
            
            rawOKNStore( jsub, : ) = [ accAllSubID(isub), state, stimCond, histRawOKN ];
            medOKNStore( jsub, : ) = [ accAllSubID(isub), state, stimCond, histMedOKN ];
            
            subPerCorr = perCorrectness( perCorrectness(:,1)==accAllSubID(isub) & ...
                perCorrectness(:,2) == state & ...
                perCorrectness(:,3) == stimCond, :);
            
            %                 figure(1),clf
            %                 subplot(2,3,1:3)
            %                 plot(rawOKN)
            %                 hold on
            %                 plot(medOKN)
            %                 ylim([-15 15])
            %
            %                 subplot(2,3,4)
            %                 hist(medOKN)
            %                 title('medOKN vel')
            %                 ylabel('nPnt')
            %                 xlabel('Velocity(deg/sec)')
            %
            %                 subplot(2,3,5:6)
            %                 text(.05, .9, [ 'state = ' stateName{state(istate)} ], 'FontSize', 13, 'interpret', 'none')
            %                 text(.05, .8, [ 'stim = ' stimName{stimCond(istimCond)} ], 'FontSize', 13)
            %                 text(.05, .7, [ 'NaN: OKN = ' num2str( perNanRaw ) '%' ], 'FontSize', 13)
            %                 text(.05, .6, [ 'Zero: OKN = ' num2str( perZeroRaw ) '%' ], 'FontSize', 13)
            %                 text(.05, .5, [ 'rejEpisode: OKN = ' num2str( subPerCorr(4)*100 ) '%' ], 'FontSize', 13)
            %                 text(.05, .4, [ 'rejEpisode: Resp = ' num2str( subPerCorr(5)*100 ) '%' ], 'FontSize', 13)
            %                 text(.05, .3, [ 'rejEpisode: Total = ' num2str( subPerCorr(6)*100 ) '%' ], 'FontSize', 13)
            %                 text(.05, .2, [ 'decodAcc: Resp = ' num2str( subPerCorr(7) ) '%' ], 'FontSize', 13)
            %                 text(.05, .1, [ 'decodeAcc: borderPos = ' num2str( subPerCorr(8) ) '%' ], 'FontSize', 13)
            %                 set(gca, 'XTick', [], 'YTick', [])
            %
            %
            %                 suptitle_m([ stateCat2 subName '_' qualAcc ])
            %
            %                 print( [ DIR.figDecodeQual '/CN_OKNQualCmp_borderPos_decAcc' qualAcc num2str(shresholdAcc)...
            %                     '_' stateCat2 subName '_state' num2str( state(istate) ) '_m10crossVal_' ...
            %                     clnName rejCriOKN '_biResp_timeWin' num2str(timeWin) ], '-dpng' )
            %
            
        end
        
    end
    
    switch iqualAcc
        case 1 % group of subject who had LOW decoding accuracy
            rawOKNStore_l   = rawOKNStore;
            medOKNStore_l   = medOKNStore;
            accAllSubInfo_l = accAllSubInfo;
            
        case 2 % group of subject who had HIGH decoding accuracy
            rawOKNStore_h   = rawOKNStore;
            medOKNStore_h   = medOKNStore;
            accAllSubInfo_h = accAllSubInfo;
    end
    
    clear accAllSubInfo rawOKNStore
end

% %% Compute rate of total missing point %%%
% accAllSubInfo_l(:,8) =  accAllSubInfo_l(:,6) + accAllSubInfo_l(:,7);
% accAllSubInfo_h(:,8) =  accAllSubInfo_h(:,6) + accAllSubInfo_h(:,7);
% % 1)subNum, 2)state, 3)stimCond, 4)decoding accuracy, 5)perNan, 6)perZero,
% % 7)perMissPnt
% 
% [H PH] = ttest2( accAllSubInfo_l(:,5), accAllSubInfo_h(:,5) )
% 
% %%% Correlation coefficient between rate od total missing point and
% %%% decoding accuracy %%%
% [R,PR] = corrcoef( [ accAllSubInfo_l(:,4); accAllSubInfo_h(:,4) ], [ accAllSubInfo_l(:,8); accAllSubInfo_h(:,8) ] )
% 
% 
% perHistRawOKN_l = rawOKNStore_l( :, 1:3 );
% perHistRawOKN_h = rawOKNStore_h( :, 1:3 );
% perHistRawOKN_l = [ perHistRawOKN_l, rawOKNStore_l(:,4:end) ./ repmat( sum(rawOKNStore_l(:,4:end),2), 1, length(binVel) ) * 100 ];
% perHistRawOKN_h = [ perHistRawOKN_h, rawOKNStore_h(:,4:end) ./ repmat( sum(rawOKNStore_h(:,4:end),2), 1, length(binVel) ) * 100 ];
% 
% accAllSubInfo_l(:,9) = perHistRawOKN_l( :, find( abs(binVel)==0 )+ 3 );
% accAllSubInfo_h(:,9) = perHistRawOKN_h( :, find( abs(binVel)==0 )+ 3 );
% 
% accAllSubInfo_l(:,10) = sum( perHistRawOKN_l( :, find( abs(binVel)>=2 & abs(binVel)<=5 ) + 3 ), 2);
% accAllSubInfo_h(:,10) = sum( perHistRawOKN_h( :, find( abs(binVel)>=2 & abs(binVel)<=5 ) + 3 ), 2);
% 
% accAllSubInfo_l(:,11) = sum( perHistRawOKN_l( :, find( abs(binVel) > 10 ) + 3 ), 2);
% accAllSubInfo_h(:,11) = sum( perHistRawOKN_h( :, find( abs(binVel) > 10 ) + 3 ), 2);
% 
% accAllSubInfo_all = cat(1, accAllSubInfo_l, accAllSubInfo_h);
% accAllSubInfo_all = sortrows( accAllSubInfo_all, 4);
% 
% accAllSubInfo_all(:,12) = sum( [ accAllSubInfo_all(:,5)>4,...
%     accAllSubInfo_all(:,8)>20,...
%     accAllSubInfo_all(:,9)>25,...
%     accAllSubInfo_all(:,10)<30, ...
%     accAllSubInfo_all(:,11)>2 ], 2 );
% 
% accAllSubInfo_all(:,12) = accAllSubInfo_all(:,12) .* (1 - [accAllSubInfo_all(:,11)<1] );
% 
% badSub = accAllSubInfo_all( accAllSubInfo_all(:,12)>=2, 1 );
% %%% Correlation between decoding accuracy and rate of OKN velocity in a
% %%% certain range
% corrcoef([ accAllSubInfo_l(:,8); accAllSubInfo_h(:,8) ], [ accAllSubInfo_l(:,4); accAllSubInfo_h(:,4) ])
% 
% nsub_l = size( perHistRawOKN_l, 1 );
% nsub_h = size( perHistRawOKN_h, 1 );
% 
% mperHistRawOKN_l = mean(perHistRawOKN_l(:,4:end),1);
% mperHistRawOKN_h = mean(perHistRawOKN_h(:,4:end),1);
% 
% semperHistRawOKN_l = std(perHistRawOKN_l(:,4:end),1)/sqrt(nsub_l);
% semperHistRawOKN_h = std(perHistRawOKN_h(:,4:end),1)/sqrt(nsub_h);
% 
% figure(2),clf
% hl = shadedErrorBar( binVel, mperHistRawOKN_l, semperHistRawOKN_l, 'b', 1 );
% hold on
% hh = shadedErrorBar( binVel, mperHistRawOKN_h, semperHistRawOKN_h, 'r', 1 );
% hl.mainLine.LineWidth = 2;
% hh.mainLine.LineWidth = 2;
% legend( [ hl.mainLine, hh.mainLine], ...
%     {['LowAcc n=' num2str( nsub_l ) ],...
%     ['HighAcc n=' num2str( nsub_h ) ]} )
% ylabel('Rate (%)')
% xlabel('Velocity (deg/sec)')
% title( ['Rate of OKN vel sample in high/low dec acc sub (p=' num2str(PH) ')'] )
% 
% 
% print( gcf, [ DIR.figDecodeQual '/CN_OKNQualCmpHist_decRespAcc' num2str(shresholdAcc)...
%     '_m10crossVal_' clnName rejCriOKN '_biResp_timeWin' num2str(timeWin) ], '-dpng' )


%%
for iResp = 1:2 % 1)button press or 2)border position (only in replay)
    switch iResp
        case 1
            respName = 'resp';
            mperCorrectnessStore = mperCorrectnessStore_resp;
        case 2
            respName = 'borderPos';
            mperCorrectnessStore = mperCorrectnessStore_borderPos;
    end
    
    if F.subRej
        respName = [ respName '_subRej'];
        mperCorrectnessStore( ismember( mperCorrectnessStore(:,1), badSub ), : ) = [];
    else
        respName = [ respName '_noSubRej'];
    end
    
    for icmp = 1:4
        
        switch icmp
            case 1 % C vs PD
                grp1 = 1;
                grp2 = 2:5;
                grp1Name  = 'C';
                grp2Name  = 'PD';
                titleName = 'C_vs_PD';
            case 2 % med on vs med off
                grp1 = 2;
                grp2 = 3;
                grp1Name  = 'med-on';
                grp2Name  = 'med-off';
                titleName = 'med-on_vs_med-off';
            case 3 % dbs on vs dbs off
                grp1 = 4;
                grp2 = 5;
                grp1Name = 'DBS-on';
                grp2Name = 'DBS-off';
                titleName = 'DBS-on_vs_DBS-off';
            case 4 % C vs PD on
                grp1 = 1;
                grp2 = [2,4];
                grp1Name = 'C';
                grp2Name = 'PD-on';
                titleName = 'C_vs_PD-on';
        end
        
        stimCond = unique(mperCorrectnessStore(:,3));
        
        for istimCond = 1:2 %BR and Replay
            
            if iResp == 2 && istimCond == 1 
                % if response is detected based on border position, which
                % is avalable only in replay not BR condition,
                % but it is stimulus condition is BR and 
                continue
            end
            
            switch stimCond(istimCond)
                case 1
                    stimName = 'br';
                case 4
                    stimName = 'replay';
            end
            
            perCorr1 = mperCorrectnessStore( ismember( mperCorrectnessStore(:,2), grp1 )...
                & mperCorrectnessStore(:,3)==stimCond(istimCond), : );
            perCorr2 = mperCorrectnessStore( ismember( mperCorrectnessStore(:,2), grp2 )...
                 & mperCorrectnessStore(:,3)==stimCond(istimCond), : );
            
            
            switch icmp
                case {1,4} % unpaired t-test
                    [H,P,CI,STATS] = ttest2( perCorr1(:,4), perCorr2(:,4) );
                case {2,3} % paired t-test
                    perCorr1( ~ismember( perCorr1(:,1), perCorr2(:,1) ), : ) = [];
                    perCorr2( ~ismember( perCorr2(:,1), perCorr1(:,1) ), : ) = [];
                    [H,P,CI,STATS] = ttest(perCorr1(:,4),perCorr2(:,4));
            end
            
            mPerCorr1   = mean(perCorr1(:,4));
            mPerCorr2   = mean(perCorr2(:,4));
            medPerCorr1 = median(perCorr1(:,4));
            medPerCorr2 = median(perCorr2(:,4));
            sdPerCorr1  = std(perCorr1(:,4));
            sdPerCorr2  = std(perCorr2(:,4));
            nPerCorr1   = length(perCorr1(:,4));
            nPerCorr2   = length(perCorr2(:,4));
            
            %%% -----------------------------------------------------------------
            figure(1),clf
            
            subplot(2,2,1)
            his = hist(perCorr1(:,4), 15);
            hist(perCorr1(:,4), 15)
            hold on
            % average
            h = plot([mPerCorr1 mPerCorr1], [0 max(his)+1]);
            h.LineWidth = 2;
            h.Color = [.8 0 .8];
            % median
            h = plot([medPerCorr1 medPerCorr1], [0 max(his)+1]);
            h.LineWidth = 2;
            h.Color = [0 .8 .8];
            
            legend('Histogram', 'Average', 'Median', 'Location', 'NorthWest')
            legend boxoff
            
            ylim([0 max(his)+1])
            xlim([25 100])
            title([ grp1Name ' (N=' num2str(nPerCorr1) ')'],'interpret', 'none')
            xlabel('Decoding Accuracy (%)')
            ylabel('nSubject')
            
            
            
            subplot(2,2,2)
            text(.05,.8,titleName, 'FontSize', 13, 'interpret', 'none')
            text(.05, .7, [ 'p = ' num2str(P) ], 'FontSize', 13)
            text(.05, .6, [ 'tstat = ' num2str(STATS.tstat) ], 'FontSize', 13)
            text(.05, .5, [ 'df = ' num2str(STATS.df)  ], 'FontSize', 13)
            text(.05, .4, [ 'sd = ' num2str(STATS.sd)  ], 'FontSize', 13)
            set(gca, 'XTick', [], 'YTick', [])
            
            
            subplot(2,2,3)
            his = hist(perCorr2(:,4), 15);
            hist(perCorr2(:,4), 15)
            hold on
            % average
            h = plot([mPerCorr2 mPerCorr2], [0 max(his)+1]);
            h.LineWidth = 2;
            h.Color = [.8 0 .8];
            
            % median
            h = plot([medPerCorr2 medPerCorr2], [0 max(his)+1]);
            h.LineWidth = 2;
            h.Color = [0 .8 .8];
            
            legend('Histogram', 'Average', 'Median', 'Location', 'NorthWest')
            legend boxoff
            
            ylim([0 max(his)+1])
            xlim([25 100])
            title([ grp2Name ' (N=' num2str(nPerCorr2) ')'], 'interpret', 'none')
            xlabel('Decoding Accuracy (%)')
            ylabel('nSubject')
            
            
            %         subplot(2,2,4)
            %         h = bar( [mPerCorr1, medPerCorr1; medPerCorr2, mPerCorr2] );
            %         hold on
            %         d = .15;
            %         errorbar( [ 1-d, 2-d ], [ mPerCorr1, mPerCorr2 ], [ sdPerCorr1, sdPerCorr2 ], '.k', 'LineWidth', 2);
            %         xticklabels( { grp1Name, grp2Name } )
            
            subplot(2,2,4)
            h1 = bar( 1, mPerCorr1);
            hold on
            h2 = bar( 2, mPerCorr2);
            h1.FaceColor = [0 0 .8];
            h1.BarWidth = .5;
            h2.FaceColor = [.8 0 0];
            h2.BarWidth = .5;
            d = 0;
            errorbar( [ 1-d, 2-d ], [ mPerCorr1, mPerCorr2 ], [ sdPerCorr1, sdPerCorr2 ], '.k', 'LineWidth', 2);
            xlim([0 3])
            ylim([ 50 100 ])
            xticklabels( { grp1Name, grp2Name } )
            
            
            %         legend(h, { 'average', 'median' } )
            %         legend boxoff
            
            suptitle_m([ titleName '_' respName '_' stimName ])
            
            if F.median
                print( [ DIR.figCTGrpCmp '/CN_decodingAccuracyCmp_m10crossVal_' clnName rejCriOKN '_biResp_timeWin' num2str(timeWin) ...
                    '_medFilt' num2str(filtWin) 'pnt_' titleName '_' respName '_' stimName], '-dpng' )
            else
                print( [ DIR.figCTGrpCmp '/CN_decodingAccuracyCmp_m10crossVal_' clnName rejCriOKN '_biResp_timeWin' num2str(timeWin) ...
                    '_' titleName '_' respName '_' stimName ], '-dpng' )
            end
            
            %
            %         if F.median
            %             print( [ DIR.figDecode '/CN_decodingAccuracyCmp_m10crossVal_' clnName rejCriOKN '_biResp_timeWin' num2str(timeWin) ...
            %                 '_medFilt' num2str(filtWin) 'pnt_' titleName '_' ansName ], '-dpng' )
            %         else
            %             print( [ DIR.figDecode '/CN_decodingAccuracyCmp_m10crossVal_' clnName rejCriOKN '_biResp_timeWin' num2str(timeWin) ...
            %                 '_' titleName '_' ansName ], '-dpng' )
            %         end
        end
    end
end

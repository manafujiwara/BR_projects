function decode_SVM_CN_eye_a_BP_allSub_cmpEyeBP(DIR, subID)
rejCriOKN = 'NanZero';
%1) NaN or 2)NanZero; Reject episodes that contain NaN (or NaN and zero)
%more than hald of it %%%

F.median = 1;
if F.median
    filtWin = 30;
end
F.biResp = 1;
F.cleanCriteria = 'All'; % 1) OKN or 2)All (both OKN and button press resp)

Plot.OKN_a_resp = 1;

timeWin  = 500; % ms
nSampPnt = 60/1000*timeWin;
sampWin  = -ceil(nSampPnt/2) : ceil(nSampPnt/2) ; % 30 samples (=500 ms) are used for each observation
sampWin  = sampWin(1:nSampPnt);
ises     = 0;
load( [ DIR.clnData  '/allSubInfo_badEp' rejCriOKN '.mat'], 'subInfo', 'badEp' )

for isub = 4:length(subID)
    %% Prepare data %%
    switch F.cleanCriteria
        % clnOKN_OKN: 1D=episode, 2D= 1)subNum, 2)state, 3)stimCondition,
        %             4:30)data
        case 'OKN'
            clnName = 'rejEpByOKN';
            
            if F.median
                load( [DIR.clnData '/' subID{isub} '_60Hz_medFilt' num2str(filtWin) '_clnOKN' rejCriOKN '.mat'],...
                    'clnOKN_OKN', 'clnResp_OKN', 'clnBiResp_OKN', 'clnBorderPos_OKN')
            else
                load( [DIR.clnData  '/' subID{isub} '_60Hz_clnOKN' rejCriOKN '.mat'],...
                    'clnOKN_OKN', 'clnResp_OKN', 'clnBiResp_OKN', 'clnBorderPos_OKN')
            end
            
            OKN = clnOKN_OKN;
            switch F.biResp
                case 1
                    resp = clnBiResp_OKN;
                case 0
                    resp = clnResp_OKN;
            end
            borderPos = clnBorderPos_OKN;
            
        case 'All'
            clnName = 'rejEpByAll';
            
            if F.median
                load( [DIR.clnData '/' subID{isub} '_60Hz_medFilt' num2str(filtWin) '_clnAll' rejCriOKN '.mat'],...
                    'clnOKN_all', 'clnResp_all', 'clnBiResp_all', 'clnBorderPos_all')
            else
                load( [DIR.clnData  '/' subID{isub} '_60Hz_clnAll' rejCriOKN '.mat'],...
                    'clnOKN_all', 'clnResp_all', 'clnBiResp_all', 'clnBorderPos_all')
            end
            
            OKN = clnOKN_all;
            switch F.biResp
                case 1
                    resp = clnBiResp_all;
                case 0
                    resp = clnResp_all;
            end
            borderPos = clnBorderPos_all;
    end
    
    subNum   = str2num( subID{isub}(2:end) );
    state    = unique( OKN(:,2) );
    stimCond = unique( OKN(:,3) );
    
    for istate = 1:length(state)
        
        for istimCond = 1:length(stimCond)
            ises = ises + 1;
            
            switch stimCond(istimCond)
                case 1 %BR condition
                    respCond = 1; %
                case 4 %Replay condition
                    respCond = [1,2];
            end
            
            perCorrStore_resp_OKN( ises, 1 ) = subNum;
            perCorrStore_resp_OKN( ises, 2 ) = state(istate);
            perCorrStore_resp_OKN( ises, 3 ) = stimCond(istimCond);
            perCorrStore_resp_resp( ises, 1 ) = subNum;
            perCorrStore_resp_resp( ises, 2 ) = state(istate);
            perCorrStore_resp_resp( ises, 3 ) = stimCond(istimCond);
            
            perCorrStore_borderPos_OKN( ises, 1 ) = subNum;
            perCorrStore_borderPos_OKN( ises, 2 ) = state(istate);
            perCorrStore_borderPos_OKN( ises, 3 ) = stimCond(istimCond);
            perCorrStore_borderPos_resp( ises, 1 ) = subNum;
            perCorrStore_borderPos_resp( ises, 2 ) = state(istate);
            perCorrStore_borderPos_resp( ises, 3 ) = stimCond(istimCond);
            
            %%% Pick up data that suit for sitmulus condition %%%
            % 1)subNum, 2)state, 3)stimCondition, 4:30)data
            pOKN       = OKN( OKN(:,2) == state(istate) & OKN(:,3) == stimCond(istimCond), 4:end );
            pResp      = resp( resp(:,2) == state(istate) & resp(:,3) == stimCond(istimCond), 4:end );
            pBorderPos = borderPos( borderPos(:,2) == state(istate) & borderPos(:,3) == stimCond(istimCond), 4:end );
            
            % randomly pick up 70% of the whole block for training, 30% for test
            nEpisode    = size( pOKN, 1 );
            nSamp_train = floor( nEpisode*.7 );
            nSamp_test  = nEpisode - nSamp_train;
            
            for iResp = 1:length(respCond) % 1)button response, 2)border position (replay only)
                
                for icrossVal = 1:10 % cross validation
                    randIdx = randperm(nEpisode);
                    t_train = randIdx( 1 : nSamp_train )';
                    t_test  = randIdx( nSamp_train + 1 : nSamp_train + nSamp_test )';
                    
                    OKN_train       = pOKN( t_train, : )' ;
                    resp_train      = pResp( t_train, : )' ;
                    borderPos_train = pBorderPos( t_train, : )' ;
                    
                    OKN_test       = pOKN( t_test, : )' ;
                    resp_test      = pResp( t_test, : )' ;
                    borderPos_test = pBorderPos( t_test, : )' ;
                    
                    
                    switch iResp
                        case 1 %predict button press response
                            ans_train = resp_train;
                            ans_test  = resp_test;
                        case 2 %predict border position
                            ans_train = borderPos_train;
                            ans_test  = borderPos_test;
                    end
                    
                    
                    [Mdl_OKN,FitInfo]= fitclinear( OKN_train, ans_train(nSampPnt/2,:),...
                        'ObservationsIn', 'columns');
                    ans_pred_OKN = predict(Mdl_OKN, OKN_test');
                    ans_real_OKN = ans_test(nSampPnt/2,:);
                    
                    [Mdl_resp,FitInfo]= fitclinear( resp_train, ans_train(nSampPnt/2,:),...
                        'ObservationsIn', 'columns');
                    
                    ans_pred_resp = predict(Mdl_resp, resp_test');
                    ans_real_resp = ans_test(nSampPnt/2,:);
                    
%                     % For an exemplar plot for the poster in NIPS
%                     load( [DIR.excel '/' subID{isub} '_60Hz.mat'], 'd_CT_mat' )
%                     t = 1:1/60:240;
%                     time = [t(min(t_test)):1/60:t(max(t_test))] - t(min(t_test));
%                     figure(3),clf
%                     plot( time, d_CT_mat(t_test,6),'LineWidth',2)
%                     hold on
%                     plot( time, medfilt1(d_CT_mat(t_test,6),30),'LineWidth',2)
%                     ylim([-15 15])
%                     
%                     plot( time, ans_pred*13,'.k')
%                     plot( time, ans_real*15,'.k')
%                     xlim([time(1) time(end)])
                    
                    
                    corr_OKN = ans_pred_OKN - ans_real_OKN';
                    corr_OKN( corr_OKN ~= 0 ) = nan;
                    corr_OKN( ~isnan(corr_OKN) ) = 1;
                    perCorr_OKN = nansum(corr_OKN) / length(~isnan(corr_OKN)) * 100;
                    
                    corr_resp = ans_pred_resp - ans_real_resp';
                    corr_resp( corr_resp ~= 0 ) = nan;
                    corr_resp( ~isnan(corr_resp) ) = 1;
                    perCorr_resp = nansum(corr_resp) / length(~isnan(corr_resp)) * 100;
                    
                    switch iResp
                        case 1
                            perCorrStore_resp_OKN( ises, icrossVal + 3 ) = perCorr_OKN;
                            perCorrStore_resp_resp( ises, icrossVal + 3 ) = perCorr_resp;
                        case 2
                            perCorrStore_borderPos_OKN( ises, icrossVal + 3 ) = perCorr_OKN;
                            perCorrStore_borderPos_resp( ises, icrossVal + 3 ) = perCorr_resp;
                    end
                    
                    clear Mdl correctness perCorrectness
                end
                
                clear t_train t_test
                
            end
        end
        
        
        figure(2),clf
        t = 1:1/60:240;
        time = [t(min(t_test)):1/60:t(max(t_test))] - t(min(t_test));
        plot( time( 1:length(OKN_train) ), OKN_train(15,:))
        hold on
        plot( time( 1:length(OKN_train) ), resp_train*10)
        
        
        % Plot examplar sample for check
        if Plot.OKN_a_resp
            figure(1), clf
            time = linspace(0, length(OKN)/60, length(OKN));
            plot(time,OKN, 'LineWidth', 2) %smooted OKN vel
            hold on
            plot(time,resp*10) %alternative response
            plot(time,zeros(1,length(OKN))) % line for vel 0
            title([ subID{isub} ' OKN & binary response (%correct=' num2str(perCorrectness) ')' ],'FontSize',16)
            ylim([-30 30])
            leg = legend( 'OKN', 'resp', 'Location', 'EastOutside');
            leg.FontSize = 16;
            xlabel('time (sec)')
            print([ DIR.figAppearance subID{isub} '_OKN_a_resp.png'], '-dpng')
        end
        
        
        figure(3),clf
        plot( time( 1:length(t_test) ), resp_real )
        hold on
        plot( time( 1:length(t_test) ), resp_pred-3)
        plot( time( 1:length(t_test) ), correctness*-1.5, 'k.')
        leg = legend('real resp', 'pred resp', 'correctness',...
            'Location', 'EastOutside');
        leg.FontSize = 16;
        axis tight
        xlabel('time (sec)')
        title( [ subID{isub} ' decode resp from OKN (%correct=' num2str(perCorrectness) ')' ], 'FontSize', 16 )
        
        
    end
    disp([ 'Done ' subID{isub} ])
end

mperCorrectnessStore_resp      = [ perCorrStore_resp_OKN(:,1:3), mean(perCorrStore_resp_OKN(:,4:end),2) ];
mperCorrectnessStore_borderPos = [ perCorrStore_borderPos_OKN(:,1:3), mean(perCorrStore_borderPos_OKN(:,4:end),2) ];

perCorr_OKN = [ subInfo, mperCorrectnessStore_resp(:,end), mperCorrectnessStore_borderPos(:,4) ];

if F.median
    save([ DIR.allSub '/CN_subDecodePerCorrectness_m10crossVal_' clnName rejCriOKN '_biResp_timeWin' num2str(timeWin) ...
        '_medFilt' num2str(filtWin) 'pnt.mat' ], 'mperCorrectnessStore_resp', 'mperCorrectnessStore_borderPos', 'perCorrectness'  )
else
    save([ DIR.allSub '/CN_subDecodePerCorrectness_m10crossVal_' clnName rejCriOKN '_biResp_timeWin' num2str(timeWin)...
        '.mat' ], 'mperCorrectnessStore_resp', 'mperCorrectnessStore_borderPos', 'perCorrectness' )
end




function decode_SVM_CN_eye_a_BP_allSub(DIR, subID)
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
load( [ DIR.clnData '/allSubInfo_badEp' rejCriOKN '.mat'], 'subInfo', 'badEp' )

for isub = 2:length(subID)
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
                load( [DIR.clnData '/' subID{isub} '_60Hz_clnOKN' rejCriOKN '.mat'],...
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
                load( [DIR.clnData '/' subID{isub} '_60Hz_clnAll' rejCriOKN '.mat'],...
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
            
            perCorrectnessStore_resp( ises, 1 ) = subNum;
            perCorrectnessStore_resp( ises, 2 ) = state(istate);
            perCorrectnessStore_resp( ises, 3 ) = stimCond(istimCond);
            perCorrectnessStore_borderPos( ises, 1 ) = subNum;
            perCorrectnessStore_borderPos( ises, 2 ) = state(istate);
            perCorrectnessStore_borderPos( ises, 3 ) = stimCond(istimCond); 
            
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
                    %                 randIdx = randperm(nEpisode);
                    randIdx = 1:nEpisode;
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
                    
                    
                    [Mdl,FitInfo]= fitclinear( OKN_train, ans_train(nSampPnt/2,:),...
                        'ObservationsIn', 'columns');
                    
                    ans_pred = predict(Mdl, OKN_test');
                    ans_real = ans_test(nSampPnt/2,:);
                    
                    % For a exemplar plot for the poster in NIPS
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
                    
                    
                    correctness = ans_pred - ans_real';
                    correctness( correctness ~= 0 ) = nan;
                    correctness( ~isnan(correctness) ) = 1;
                    perCorrectness = nansum(correctness) / length(~isnan(correctness)) * 100;
                    
                    switch iResp
                        case 1
                            perCorrectnessStore_resp( ises, icrossVal + 3 ) = perCorrectness;
                        case 2
                            perCorrectnessStore_borderPos( ises, icrossVal + 3 ) = perCorrectness;
                    end
                    
                    clear Mdl correctness perCorrectness
                end
                
                clear t_train t_test
                
            end
        end
        %     figure(2),clf
        %     plot( time( 1:length(OKN_train) ), OKN_train(15,:))
        %     hold on
        %     plot( time( 1:length(OKN_train) ), resp_train*10)
        
        
        %% Plot examplar sample for check
        %     if Plot.OKN_a_resp
        %         figure(1), clf
        %         time = linspace(0, length(OKN)/60, length(OKN));
        %         plot(time,OKN, 'LineWidth', 2) %smooted OKN vel
        %         hold on
        %         plot(time,resp*10) %alternative response
        %         plot(time,zeros(1,length(OKN))) % line for vel 0
        %         title([ subID{isub} ' OKN & binary response (%correct=' num2str(perCorrectness) ')' ],'FontSize',16)
        %         ylim([-30 30])
        %         leg = legend( 'OKN', 'resp', 'Location', 'EastOutside');
        %         leg.FontSize = 16;
        %         xlabel('time (sec)')
        %         print([ DIR.figAppearance subID{isub} '_OKN_a_resp.png'], '-dpng')
        %     end
        
        
        %     figure(3),clf
        %     plot( time( 1:length(t_test) ), resp_real )
        %     hold on
        %     plot( time( 1:length(t_test) ), resp_pred-3)
        %     plot( time( 1:length(t_test) ), correctness*-1.5, 'k.')
        %     leg = legend('real resp', 'pred resp', 'correctness',...
        %         'Location', 'EastOutside');
        %     leg.FontSize = 16;
        %     axis tight
        %     xlabel('time (sec)')
        %     title( [ subID{isub} ' decode resp from OKN (%correct=' num2str(perCorrectness) ')' ], 'FontSize', 16 )
        %
        %
    end
    disp([ 'Done ' subID{isub} ])
end

mperCorrectnessStore_resp      = [ perCorrectnessStore_resp(:,1:3), mean(perCorrectnessStore_resp(:,4:end),2) ];
mperCorrectnessStore_borderPos = [ perCorrectnessStore_borderPos(:,1:3), mean(perCorrectnessStore_borderPos(:,4:end),2) ];

perCorrectness = [ subInfo, mperCorrectnessStore_resp(:,end), mperCorrectnessStore_borderPos(:,4) ];

if F.median
    save([ DIR.allSub '/CN_subDecodePerCorrectness_m10crossVal_' clnName rejCriOKN '_biResp_timeWin' num2str(timeWin) ...
        '_medFilt' num2str(filtWin) 'pnt.mat' ], 'mperCorrectnessStore_resp', 'mperCorrectnessStore_borderPos', 'perCorrectness'  )
else
    save([ DIR.allSub '/CN_subDecodePerCorrectness_m10crossVal_' clnName rejCriOKN '_biResp_timeWin' num2str(timeWin)...
        '.mat' ], 'mperCorrectnessStore_resp', 'mperCorrectnessStore_borderPos', 'perCorrectness' )
end






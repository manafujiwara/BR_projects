function reject_smooth_CN(DIR,subID)
rejCriOKN = 'NanZero'; 
% =========================================================================
% 1) NaN or 2)NanZero; Reject episodes that contain NaN (or NaN and zero)
% more than half of it 
% =========================================================================

F.median = 1; % Median filtering on or off
if F.median
    filtWin = 30;
end

timeWin  = 500; % ms
nSampPnt = 60/1000*timeWin;
sampWin  = -ceil(nSampPnt/2) : ceil(nSampPnt/2) ; % 30 samples (=500 ms) are used for each observation
sampWin  = sampWin(1:nSampPnt);
ises     = 0;

for isub = 2:length(subID) % start from isub=2 because the first data is of an experiemnter
    clnOKN_OKN       = [];
    clnResp_OKN      = [];
    clnBiResp_OKN    = [];
    clnBorderPos_OKN = [];
    clnOKN_all       = [];
    clnResp_all      = [];
    clnBiResp_all    = [];
    clnBorderPos_all = [];
    
    subNum = str2num( subID{isub}(2:end) );
    %% Prepare data
    % load( [DIR.excel '/' subID{isub} '.mat'], 'CT_mat')
    % CT_mat: 1)subNum, 2)state, 3)time, 4)stimCondition, 5)lResp, 6)rResp, 7:8)iOKN, ...
    % 9:10)vel_smtiOKN, 11:12)smtvel_smtiOKN, 13)borderPos, 14)borderPos_lb , 15:16)left/right pupil diameter ] );
    
    load( [DIR.excel '/' subID{isub} '_60Hz.mat'], 'd_CT_mat' )
    % d_CT_mat: 1)subNum, 2)state, 3)time, 4)stimCondition, 5)alternative resp(just subtracted), 6)smothed OKN vel, ...
    % 7)borderPos_lb , 8:9)left/right pupil diameter ] );
    state    = unique( d_CT_mat(:,2) );
    stimCond = unique( d_CT_mat(:,4) );
    
    for istate = 1:length(state)
        
        for istimCond = 1:length(stimCond)
            ises = ises + 1;
            
            %%% pick up data that suit for sitmulus condition %%%
            d_CT_cnd = d_CT_mat( d_CT_mat(:,2) == state(istate) & d_CT_mat(:,4) == stimCond(istimCond), : );
            
            OKN       = d_CT_cnd( :, 6);
            resp      = d_CT_cnd( :, 5);
            borderPos = d_CT_cnd( :, 7);
            
            %% Set binary response
            
            zeroIdx = find(resp==0);
            tmp = diff(zeroIdx);
            endZeroIdx = find(diff(zeroIdx)~=1);
            swapZeroIdx = [ [ 1;endZeroIdx(1:end-1)+1 ], endZeroIdx ];
            swapZeroIdx = [ swapZeroIdx; [ swapZeroIdx(end)+1, length(zeroIdx) ]];
            
            biResp = resp; % preallocation
            for iswap = 1:size(swapZeroIdx,1)
                
                if zeroIdx( swapZeroIdx(iswap,2) ) < length(biResp)
                    biResp( zeroIdx( swapZeroIdx(iswap,1) ): zeroIdx( swapZeroIdx(iswap,2) ) ) = ...
                        biResp( zeroIdx( swapZeroIdx(iswap,2) ) + 1 ); % copy resp from after the end of zero time
                else
                    biResp( zeroIdx( swapZeroIdx(iswap,1) ): zeroIdx( swapZeroIdx(iswap,2) ) ) = nan;
                    % if BP in very end of the trial is missing, put nan
                end
            end
            
            %%
            sampIdx = [ 1:size( d_CT_cnd, 1 ) ]'; %index for pickIdx
            obsIdx  = [ repmat( sampIdx, [ 1, length(sampWin) ] ) + repmat( sampWin, [ length(sampIdx), 1 ] ) ];
            avlPickIdx = obsIdx( obsIdx(:,1)>0 & obsIdx(:,end)<=length(sampIdx), : ) ;
            % remove negative index and pick up only available time points
            % : this happens in the begining and ending of a block
            % because we use samples before and after a certain time point
            
            if F.median
                OKN = medfilt1(OKN,filtWin);
            end
            
            matOKN       = reshape( OKN( avlPickIdx ), [size(avlPickIdx)] );
            matResp      = reshape( resp( avlPickIdx ), [size(avlPickIdx)] );
            matBiResp    = reshape( biResp( avlPickIdx ), [size(avlPickIdx)] );
            matBorderPos = reshape( borderPos( avlPickIdx ), [size(avlPickIdx)] );
            
            %%% episodes whose data contain NaN/zero more than half of the
            %%% episode
            
            switch rejCriOKN
                case 'Nan'
                    badEpOKN = find( sum( isnan(matOKN), 2 ) > nSampPnt/2 );
                case 'NanZero'
                    badEpOKN = find( sum( isnan(matOKN) | matOKN == 0, 2 ) > nSampPnt/2 );
            end
            
            badEpResp = find( sum( isnan(matResp), 2 ) > nSampPnt/2 | sum( matResp==0, 2 ) > nSampPnt/2 );
            badEpAll  = unique( [ badEpOKN; badEpResp ] );
            nEpisode  = size( matOKN, 1 );
            
            subInfo( ises, 1 ) = subNum;
            subInfo( ises, 2 ) = state(istate);
            subInfo( ises, 3 ) = stimCond(istimCond);
            subInfo( ises, 4 ) = length(badEpOKN)/nEpisode;
            subInfo( ises, 5 ) = length(badEpResp)/nEpisode; % rate of bad trials
            subInfo( ises, 6 ) = length(badEpAll)/nEpisode;
                        
            badEp(ises).subID     = subNum;
            badEp(ises).state     = state(istate);
            badEp(ises).stimCond  = stimCond(istimCond);
            badEp(ises).nEpisode  = nEpisode;
            badEp(ises).rejEpOKN  = badEpOKN;
            badEp(ises).rejEpResp = badEpResp;
            badEp(ises).rejEpAll  = badEpAll;
            
            

            %%
            % Reject eposodes according to OKN (this can be used for replay
            % trial that dose not require button press response)
            clnOKN_OKN_tmp       = matOKN( ~ismember( 1:size(matOKN,1), badEpOKN ), : );
            clnResp_OKN_tmp      = matResp( ~ismember( 1:size(matOKN,1), badEpOKN ), : );
            clnBiResp_OKN_tmp    = matBiResp( ~ismember( 1:size(matOKN,1), badEpOKN ), : );
            clnBorderPos_OKN_tmp = matBorderPos( ~ismember( 1:size(matOKN,1), badEpOKN ), : );
            
            % Reject eposodes according to OKN AND button press response;
            % accept only eposodes that satisfy both OKN and button press
            % quality
            clnOKN_all_tmp       = matOKN( ~ismember( 1:size(matOKN,1), badEpAll ), : );
            clnResp_all_tmp      = matResp( ~ismember( 1:size(matOKN,1), badEpAll ), : );
            clnBiResp_all_tmp    = matBiResp( ~ismember( 1:size(matOKN,1), badEpAll ), : );
            clnBorderPos_all_tmp = matBorderPos( ~ismember( 1:size(matOKN,1), badEpAll ), : );
            
            
            % Store data each subject
            clnOKN_OKN       = cat(1, clnOKN_OKN, [ repmat( subInfo( ises, 1:3 ), nEpisode - length(badEpOKN), 1 ), clnOKN_OKN_tmp ] );
            clnResp_OKN      = cat(1, clnResp_OKN, [ repmat( subInfo( ises, 1:3 ), nEpisode - length(badEpOKN), 1 ), clnResp_OKN_tmp ]);
            clnBiResp_OKN    = cat(1, clnBiResp_OKN, [ repmat( subInfo( ises, 1:3 ), nEpisode - length(badEpOKN), 1 ), clnBiResp_OKN_tmp ]);
            clnBorderPos_OKN = cat(1, clnBorderPos_OKN, [ repmat( subInfo( ises, 1:3 ), nEpisode - length(badEpOKN), 1 ), clnBorderPos_OKN_tmp ]);
    
            clnOKN_all       = cat(1, clnOKN_all, [ repmat( subInfo( ises, 1:3 ), nEpisode - length(badEpAll), 1 ), clnOKN_all_tmp ]);
            clnResp_all      = cat(1, clnResp_all, [ repmat( subInfo( ises, 1:3 ), nEpisode - length(badEpAll), 1 ), clnResp_all_tmp ]);
            clnBiResp_all    = cat(1, clnBiResp_all, [ repmat( subInfo( ises, 1:3 ), nEpisode - length(badEpAll), 1 ), clnBiResp_all_tmp ]);
            clnBorderPos_all = cat(1, clnBorderPos_all, [ repmat( subInfo( ises, 1:3 ), nEpisode - length(badEpAll), 1 ), clnBorderPos_all_tmp ]);

        end
    end
    
    if F.median
        save( [DIR.clnData '/' subID{isub} '_60Hz_medFilt' num2str(filtWin) '_clnOKN' rejCriOKN '.mat'],...
            'clnOKN_OKN', 'clnResp_OKN', 'clnBiResp_OKN', 'clnBorderPos_OKN')
        save( [DIR.clnData '/' subID{isub} '_60Hz_medFilt'  num2str(filtWin) '_clnAll' rejCriOKN '.mat'],...
            'clnOKN_all', 'clnResp_all', 'clnBiResp_all', 'clnBorderPos_all')
    else
        save( [DIR.clnData '/' subID{isub} '_60Hz_clnOKN' rejCriOKN '.mat'],...
            'clnOKN_OKN', 'clnResp_OKN', 'clnBiResp_OKN', 'clnBorderPos_OKN')
        save( [DIR.clnData '/' subID{isub} '_60Hz_clnAll' rejCriOKN '.mat'],...
            'clnOKN_all', 'clnResp_all', 'clnBiResp_all', 'clnBorderPos_all')
    end
    
    disp([ 'Done ' subID{isub} ])
end

save( [ DIR.clnData '/allSubInfo_badEp' rejCriOKN '.mat'], 'subInfo', 'badEp' )


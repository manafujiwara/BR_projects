function rej_upsmp_smth_IM(DIR, subID)

load([DIR.allSub '/IM_BPinfo.mat'], 'BPInfo', 'respInfo')
load([DIR.allSub '/IM_eyeinfo.mat'], 'slwPhsBoxInfo', 'iOKNInfo', 'velInfo')

eye30info  = slwPhsBoxInfo;
% eye150info = [velInfo(:,1:10), boxcar(velInfo(:,11:910),150) ];

% -----------------------------------------------------------------------
% BPCln: 1)subID, 2)irun, 3)iblk, 4)iepc, 5)stimLabel,                 |
% 6)subjestState, 7)BPP (-1~1), 8)label, 9)fBP, 10)fcBP, 11:910)BP(t)   |
% -----------------------------------------------------------------------
% eyeinfo: 1)subID, 2)irun, 3)iblk, 4)iepc, 5)stimLabel, 6)subjestState |
% 7)BPP (-1~1), 8)label, 9)fBP, 10)fcBP, 11:910)eyeinfo(t)              |
% -----------------------------------------------------------------------

BPCln = []; respCln = []; slwPhsBoxCln = []; slwPhsBox150Cln = []; iOKNCln = []; velCln = [];
jsub = 0;
for isub = 1: length(subID)
    subNum = str2num(subID{isub}(2:4));
    
    BP_sub = BPInfo( BPInfo(:,1) == subNum, : );
    
    state = unique(BP_sub(:,6));
    nstate = length(state);
    
    for istate = 1:nstate
        jsub = jsub + 1;
        BP = BP_sub(BP_sub(:,6)==state(istate),:);
        %         nallTri = size(BP_sub,1);
        nallTri = size(BP,1);
        
        br_BP = BP( BP(:,5) == 1, :);
        nonbr_BP = BP( ismember(BP(:,5),[2,3]), :);
        
        %%% check balance of right and left
        br_rBP = find( BP(:,8)== 1 & BP(:,5) == 1 );
        br_lBP = find( BP(:,8)== -1 & BP(:,5) == 1 ) ;
        
        %%% check num of bad trials
        eye30 = eye30info( eye30info(:,1) == subNum & eye30info(:,6)==state(istate), : );
        
        badTriBP    = find( sum( BP(:,311:910) == 0 | isnan( BP(:,311:910) ), 2) > 300 | isnan( BP(:,7)) | isnan( BP(:,8) ) );
        badTrieye30 = find( sum( eye30(:,311:910) == 0 ,2) > 300 | sum( isnan( eye30(:,311:910) ), 2 ) > 300 | isnan( eye30(:,7) ) );
        rejBlockTri = find( isnan( eye30(:,5) ) );
%         badTriBP( ismember(badTriBP, rejBlockTri) ) = []; 
%         % remove trials in removed BLOCK because they are rejected due to
%         % gaze data not BP
        
        badTri(jsub).BPeye30 = unique( [badTriBP', badTrieye30'] ); %this includes trial rejected by block rejection
        
        rejBPeye30 = length(badTri(jsub).BPeye30)/nallTri ; % percentage of trial rejection
        
        val_br_rBP = sum( ~ismember( br_rBP, badTri(jsub).BPeye30 ) );
        val_br_lBP = sum( ~ismember( br_lBP, badTri(jsub).BPeye30 ) );
        
        if isempty( val_br_rBP )
            val_br_rBP = nan;
        end
        
        if isempty( val_br_lBP )
            val_br_lBP = nan;
        end
        
        subInfo(jsub,:) = [ subNum, state(istate), nallTri, length(badTri(jsub).BPeye30),...
            rejBPeye30, val_br_rBP, val_br_lBP, length(badTriBP), length(badTrieye30),...
            length(br_rBP), length(br_lBP), min( [ length(br_lBP),length(br_rBP)] )/max([length(br_lBP), length(br_rBP)]) ];
        
        
        %%% Remove trials %%%
        BP(badTri(jsub).BPeye30,:) = [];
        eye30(badTri(jsub).BPeye30,:) = [];
        
%         iOKN = iOKNInfo(iOKNInfo(:,1)==subNum , : );
%         iOKN(badTri(jsub).BPeye30,:)=[];
%         
%         vel = velInfo(velInfo(:,1)==subNum , : );
%         vel(badTri(jsub).BPeye30,:)=[];
        
        BPCln           = cat(1, BPCln, BP);
        slwPhsBoxCln    = cat(1, slwPhsBoxCln, eye30);
%         iOKNCln         = cat(1, iOKNCln, iOKN);
%         velCln          = cat(1, velCln, vel);
    end
end


rejSub = unique( subInfo( find( subInfo (:,5) > .5 | subInfo (:,6) < 3 | subInfo (:,7)< 3), 1));
accSub = unique( subInfo( ~ismember( subInfo(:,1), rejSub(:,1)), 1 ));

BPCln(ismember(BPCln(:,1),rejSub),:)               = [];
slwPhsBoxCln(ismember(slwPhsBoxCln(:,1),rejSub),:) = [];
% iOKNCln(ismember(iOKNCln(:,1),rejSub),:)           = [];
% velCln(ismember(velCln(:,1),rejSub),:)             = [];

save( [ DIR.allSub '/subInfo.mat' ],'rejSub', 'accSub', 'subInfo', 'badTri')
save( [DIR.allSub '/IM_BPCln.mat'], 'BPCln')
save( [DIR.allSub '/IM_eyeCln.mat'], 'slwPhsBoxCln', 'iOKNCln', 'velCln', 'slwPhsBox150Cln')


%% Upsmple, Smooth %%%
tw = 31; % time window (time points, in this case trials)
ms = 1; % moving step
time = linspace(-1, 2, 900);
colormap jet
jsub = 0;

for isub = 1: length(subID)
    subNum = str2num(subID{isub}(2:4));
    
    if ~ismember(subNum, accSub)
        continue
    end
    
    jsub = jsub + 1;
    
    BP     = BPCln( BPCln(:,1) == subNum, : );
    eye30  = slwPhsBoxCln( slwPhsBoxCln(:,1) == subNum, : );
    
    % ----------------------------------------------------------------------
    % BPCln: 1)subID, 2)irun, 3)iblk, 4)iepc, 5)stimLabel,                 |
    % 6)subjestState, 7)pBP (-1~1), 8)label, 9)fBP, 10)fcBP, 11:910)BP(t)  |
    % ----------------------------------------------------------------------
    
    % Flip non-rivalrous leftwards and rivalrous left-labeled
    BP( BP(:,5)==3, [7,11:910] ) = -1*BP( BP(:,5)==3, [7,11:910] );
    eye30( eye30(:,5)==3, [7,11:910] ) = -1*eye30( eye30(:,5)==3, [7,11:910] );
    
    BP( BP(:,5)==1 & BP(:,8)==-1, [7,11:910] ) = -1*BP( BP(:,5)==1 & BP(:,8)==-1, [7,11:910] );
    eye30( eye30(:,5)==1 & eye30(:,8)==-1, [7,11:910] ) = -1*eye30( eye30(:,5)==1 & eye30(:,8)==-1, [7,11:910] );
    
    % sort based on pBP (Proportion of Button Press)
    BPval     = flipud(sortrows( BP, 7 ));
    eye30val  = flipud(sortrows( eye30, 7 ));
    
    pBP = [ min( eye30val( ~isnan( eye30val( :, 7 ) ), 7 ) ), max( eye30val( ~isnan( eye30val( :,7 ) ),7 ) ) ];
    
    for idata = 1:2
        
        switch idata
            case 1
                data = BPval;
                titleD = '_BP';
            case 2
                data = eye30val;
                titleD = '_eye30';
        end
        
        for istim = 1:2
            
            switch istim
                case 1 % non-rivalrous
                    data_put = data( ismember( data(:,5), [2,3] ), 11:910 );
                    titleStim = '_non';
                case 2 % br
                    data_put = data( data(:,5) == 1,  11:910 );
                    titleStim = '_be';
            end
            
            switch idata
                case 1 % data = BP
                    yLim = [0 1];
                case 2 % data = eye
                    switch istim
                        case 1 % stim = non-br
                            yLim = [0 4];
                        case 2 % stim = br
                            yLim = [0 2];
                    end
            end
            
            %%% Upsampling %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            idx = sum( isnan( data_put ),2) == 0;
            nvalid = sum(idx);
            X  = repmat( time,[nvalid, 1]);
            Y  = repmat( linspace( 0, 1, nvalid )',[ 1, 900 ]);
            Z  = data_put( idx, : );
            Xq = repmat( time, 1000, 1);
            Yq = repmat( linspace( 0, 1, 1000 )',[ 1, 900 ]);
            Zq = interp2( X, Y, Z, Xq, Yq, 'cubic');
            
            
            %%% Plot BEFORE smoothing %%% (for methods)
            figure(4),clf, colormap jet
            subplot(1,2,1)
            imagesc(time,[],data_put, [yLim])
            title('BeforeSmooth')
            colorbar
            
            %%% Smoothing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            ntri = 1000;
            tmp = [ 1 : ntri - (tw - 1) ; tw : ntri]';
            vSmtIdx = [ ones( (tw-1)/2, 1) , ( tw - (tw-1)/2 : tw - 1)' ; tmp;
                ( ntri - (tw - 2): ntri - (tw - 1) + (tw-1)/2 )', ntri * ones( (tw-1)/2, 1) ];
            
            for itri = 1:ntri
                if mod(itri,200)==0
                    msg = [ 'Upsampling: sub= ' subID{isub} ', data= ' num2str(idata) ', stim=' num2str(istim) ', tri=' num2str(itri) ' / 1000' ];
                    disp(msg)
                end
                
                for itime = 1:size(time,2)
                    data_smth(itri,itime) = mean( Zq( vSmtIdx(itri,1):vSmtIdx(itri,2), itime ) );
                end
            end
            
            
            switch istim
                case 1
                    data_smth_non = data_smth;
                case 2
                    data_smth_br = data_smth;
            end
            
            figure(4)
            subplot(1,2,2)
            imagesc(time,[],data_smth, [yLim])
            title('AfterSmooth')
            colorbar
            set(figure(4),'PaperUnits','inches','PaperPosition',[0 0 6 3])
            print(figure(4), '-dpng', [DIR.figIMBP '/smooth_' num2str(subNum) titleD titleStim '.png'])
        end
        
        switch idata
            case 1
                blb_BPsmth_non(:,:,jsub) = data_smth_non;
                blb_BPsmth_br(:,:,jsub)  = data_smth_br;
            case 2
                blb_eye30smth_non(:,:,jsub) = data_smth_non;
                blb_eye30smth_br(:,:,jsub)    = data_smth_br;
        end
    end
    
    figure(1), clf, colormap(jet)
    
    for iplot = 1:4
        
        switch iplot
            case 1
                data_plot = blb_BPsmth_non(:,:,jsub) ;
                title_plot = 'non-br / BP';
                yLim = [0 1];
            case 2
                data_plot = blb_BPsmth_br(:,:,jsub) ;
                title_plot = 'br / BP';
                yLim = [0 1];
            case 3
                data_plot = blb_eye30smth_non(:,:,jsub);
                title_plot = 'non-br/ OKN30tp';
                yLim = [0 4];
            case 4
                data_plot = blb_eye30smth_br(:,:,jsub);
                title_plot = 'br/ OKN30tp';
                yLim = [0 4];
        end
        
        subplot(3,2,iplot)
        imagesc(time, [], data_plot, yLim)
        colorbar
        set(gca, 'YTick',   pBP , 'YTickLabel', fliplr(pBP))
        ylabel('pBP')
        xlabel('time (s)')
        title(title_plot)
        
    end
    suptitle_m(subID{isub})
    
    set(figure(1),'PaperUnits','inches','PaperPosition',[0 0 20 10])
    print('-dpng', [DIR.figIMSortBPSlwPhs '/' subID{isub} '_sort_BP_a_slwPhs_lb.png'] )
    
    
    timeClock = clock;
    msg = sprintf('rej_upsmp_smth_IM. Done with sub: %s, %d-%d-%d. %d:%d:%d', subID{isub}, round(timeClock));
    disp(msg);
    
end

BPsmth.non.blb     = blb_BPsmth_non;
BPsmth.br.blb      = blb_BPsmth_br;
eye30smth.non.blb  = blb_eye30smth_non;
eye30smth.br.blb   = blb_eye30smth_br;

save([DIR.allSub '/allSub_IM_rej_upsamp_BPsmth.mat'],'BPsmth', '-v7.3')
save([DIR.allSub '/allSub_IM_rej_upsamp_eye30smth.mat'],'eye30smth', '-v7.3')

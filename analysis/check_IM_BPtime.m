function check_IM_BPtime(DIR,subID)

% Check button press time %================================================
% FBP: First button press - literally the first button press
% FCBP: First continuous button press - set to avoid mistakenly pressed
% first button press
% =========================================================================

nSubGrp = 5; % Choose number of subject groups from [2, 4, 5].
% 2: Control vs patients on treatment
% 4: Control vs all patients, with- vs without-electrodes patients
% 5: Control, Med-on, Med-off, DBS-on, DBS-off

load( [DIR.allSub '/IM_BPCln.mat'])
load( [ DIR.allSub '/subInfo.mat'  ] ,'rejSub', 'accSub', 'subInfo', 'badTri')

time = linspace(-1, 2, 900);

%% HERE, FIX WRONG SUBJECT INFO
BPCln( ismember( BPCln(:,1),[162, 166] ), 6 ) = 3;

%% REMOVE SUBJECT WHO HAVE DONE ONLY ON- OR OFF- TREATMENT
% Remove those patients only when compare within subject
ksub = 0;

subID_grp = unique(BPCln(:,1));
for isub = 1:length(subID_grp)
    state = unique(BPCln(BPCln(:,1)==subID_grp(isub),6));
    if length(state)==1 && state == 1
        continue
    elseif length(state)==1 && state ~= 1
        ksub = ksub + 1;
        oneStateSub(ksub) = subID_grp(isub);
    end
end

BPCln( ismember( BPCln(:,1), oneStateSub ), : ) = [];

%%
% -------------------------------------------------------------------------
% subInfo: 1) subNum, 2)state(istate), 3)nallTri, 3)nBadTri(eye30))
%         4)rejBPeye30(rate), 5)rejBPeye30, 6)val_br_rBP, 7)val_br_lBP;
% -------------------------------------------------------------------------
% BPCln: 1)subID, 2)irun, 3)iblk, 4)iepc, 5)stimLabel,
% 6)subjestState, 7)BPP (-1~1), 8)label, 9)fBP, 10)fcBP, 11:910)BP(t)
% -------------------------------------------------------------------------

jsub = 0; ksub = 0;
FBP_sub = []; % First button press (possibly mistakenly pressed)
FCBP_sub = []; % First continuous button press
BPAcc_sub = [];

for isub = 1:length(subID)
    
    subNum = str2num(subID{isub}(2:4));
    
    if ~ismember(subNum, accSub) || ismember(subNum, oneStateSub)
        continue
    end
    
    d = BPCln;
    d = d( d(:,1) == subNum, :);
    
    for istim = 1:2 % non-rivalry or br
        
        switch istim
            
            case 1 % stim = non-br
                dd = d( ismember( d(:,5), [2,3] ) , :);
                titleStim = '_non';
                stimLabel = dd(:,5);
                stimLabel(stimLabel==2) = 1;
                stimLabel(stimLabel==3) = -1;
            case 2 % stim = br
                dd = d( d(:,5) == 1 , :);
                titleStim = '_br';
        end
        %%
        switch nSubGrp
            case 4
                FBP_sub  = [FBP_sub; nan(1,4)];
                FCBP_sub = [FCBP_sub; nan(1,4)];
                
                if istim == 1
                    BPAcc_sub = [BPAcc_sub; nan(1,4)];
                    ksub = ksub + 1;
                end
                
                jsub = jsub + 1;
                state = str2num(sprintf('%d',unique(dd(:,6))));
                
                ddd = dd;
                
                nTri = size(ddd,1);
                label = ddd(:,8);
                
                for iTri = 1:nTri
                    
                    if isnan(label(iTri))
                        FBP(iTri) = nan;
                        FCBP(iTri) = nan;
                        
                        if istim == 1
                            BPAcc(iTri) = nan;
                        end
                        
                    else
                        FBP(iTri) = find( abs( ddd(iTri,311:910) ) >= 0.5, 1, 'first')/300;
                        switch istim
                            case 1 % Unambiguous stimulus
                                iFBP = FBP(iTri)*300;
                                
                                tmp = find( ddd(iTri,311:910) * stimLabel(iTri) >= 0.5, 1, 'first')/300;
                                if isempty(tmp)
                                    FCBP(iTri) = nan;
                                    BPAcc(iTri) = 0;
                                else
                                    FCBP(iTri) = tmp;
                                    BPAcc(iTri) = length( find( ddd(iTri, iFBP+311 : 760) == stimLabel(iTri) ) )/ length(iFBP+311:760);
                                end
                                
                                % Accuracy was calculated as a rate of correctness of
                                % button press that made in between first button press and
                                % 1.5 s from the stimulus onset.
                                
                            case 2 % Umbiguous stimulus
                                FCBP(iTri) = find( ddd(iTri,311:910) * label(iTri) >= 0.5, 1, 'first')/300;
                        end
                        
                    end
                end
                
                FBP_sub(jsub,:)  = [ subNum, istim, state, nanmean(FBP)  ];
                FCBP_sub(jsub,:) = [ subNum, istim, state, nanmean(FCBP) ];
                
                if istim == 1
                    BPAcc_sub(ksub,:) = [ subNum, istim, state, nanmean(BPAcc)];
                end
                clear FBP FCBP BPAcc
                
                
            case {5,2}
                state = unique(dd(:,6));
                nstate = length(state);
                
                for istate = 1:nstate
                    FBP_sub = [FBP_sub; nan(1,4)];
                    FCBP_sub = [FCBP_sub; nan(1,4)];
                    
                    if istim == 1
                        BPAcc_sub = [BPAcc_sub; nan(1,4)];
                        ksub = ksub + 1;
                    end
                    
                    jsub = jsub + 1;
                    
                    ddd = dd( dd( :, 6 ) == state( istate ), : );
                    
                    nTri = size(ddd,1);
                    label = ddd(:,8);
                    
                    for iTri = 1:nTri
                        
                        if isnan(label(iTri))
                            FBP(iTri) = nan; % First button press
                            FCBP(iTri) = nan; % First continuous button press
                            
                            if istim == 1
                                BPAcc(iTri) = nan;
                            end
                            
                        else
                            FBP(iTri) = find( abs( ddd(iTri,311:910) ) >= 0.5, 1, 'first')/300;
                            switch istim
                                case 1
                                    iFBP = FBP(iTri)*300;
                                    
                                    tmp = find( ddd(iTri,311:910) * stimLabel(iTri) >= 0.5, 1, 'first')/300;
                                    if isempty(tmp)
                                        FCBP(iTri) = nan;
                                        BPAcc(iTri) = 0;
                                    else
                                        FCBP(iTri) = tmp;
                                        BPAcc(iTri) = length( find( ddd(iTri, iFBP+311 : 760) == stimLabel(iTri) ) )/ length(iFBP+311:760);
                                    end
                                    
                                    % Accuracy was calculated as a rate of correctness of
                                    % button press that made in between first button press and
                                    % 1.5 s from the stimulus onset.
                                    
                                case 2
                                    FCBP(iTri) = find( ddd(iTri,311:910) * label(iTri) >= 0.5, 1, 'first')/300;
                            end
                            
                        end
                    end
                    
                    FBP_sub(jsub,:) = [ subNum, istim,  state( istate ), nanmean(FBP)];
                    FCBP_sub(jsub,:) = [ subNum, istim,  state( istate ), nanmean(FCBP)];
                    
                    if istim == 1
                        BPAcc_sub(ksub,:) = [ subNum, istim,  state( istate ), nanmean(BPAcc)];
                    end
                    clear FBP FCBP BPAcc
                end
        end
    end
end
save([DIR.figIMBP '/IM_BPtime_eachSub_Grp' num2str(nSubGrp) '.mat'], 'FBP_sub', 'FCBP_sub')



switch nSubGrp
    case 4
        comGrpName = {'WO-Electrode vs W-','Conreol vs All PD'};
        
        %% STATS FOR ACCURACY
        acc_c   = BPAcc_sub(BPAcc_sub(:,3)==1, 4);
        acc_p   = BPAcc_sub(BPAcc_sub(:,3)~=1, 4);
        acc_w   = BPAcc_sub(ismember(BPAcc_sub(:,3),[4,5,45]), 4);
        acc_wo  = BPAcc_sub(ismember(BPAcc_sub(:,3),[2,3,23]), 4);
        
        [H1,P1,CI1,STATS1] = ttest2(acc_c, acc_p);
        [H2,P2,CI2,STATS2] = ttest2(acc_w, acc_wo);
        
        
        %% STATS FOR FBP and FCBP
        j = 0; iloc = 0; istats  = 0;
        figure(1),clf
        for idata = 1:2
            
            switch idata
                case 1
                    D = FBP_sub;
                    dataType = 'FBP';
                case 2
                    D = FCBP_sub;
                    dataType = 'FCBP';
            end
            
            for icmp = 1:2 % we are going to compare 2 sets of subject groups
                iloc = iloc + 1;
                
                switch icmp
                    case 1
                        grp1 = [2,3,23]; % ACTUAL NUMBER OF GROUP
                        grp2 = [4,5,45];
                        grp = [1 2]; %Number for comparison (same in "check_cmp_inR.m")
                        cmpName = {'woEle', 'wEle'};
                    case 2
                        grp1 = 1;
                        grp2 = [2,3,4,5,23,45];
                        grp = [3, 4];
                        cmpName = {'Control', 'PD'};
                end
                
                %% t-test: Compare subject group
                %lat 1)subID, 2)stim, 3)isubgrp(NOT STATE), 4)nTri, 5)data type,
                % 6)latenct
                for istim = 1:2
                    istats = istats + 1;
                    alLat_1 = D( ismember( D(:,2), istim) & ismember(D(:,3), grp1), 4 );
                    alLat_2 = D( ismember( D(:,2), istim) & ismember(D(:,3), grp2), 4 );
                    
                    [H(istats),P(istats),CI(:,istats),STATS] = ttest2(alLat_1, alLat_2) ;
                    tstat(istats) = STATS.tstat;
                    df(istats) = STATS.df;
                    sd(istats) = STATS.sd;
                    statSum(istats,:) = [idata, icmp, istim, P(istats), tstat(istats), df(istats), sd(istats)];
                    
                    j = j + 1;
                    mLat1(j,:) = [ idata, icmp, istim, grp(1), mean(alLat_1), std(alLat_1)/sqrt(length(alLat_1)) ];
                    j = j + 1;
                    mLat1(j,:) = [ idata, icmp, istim, grp(2), mean(alLat_2), std(alLat_2)/sqrt(length(alLat_2)) ];
                end
                
                lat_1 = mLat1( mLat1(:,1)==idata & mLat1(:,2)==icmp & mLat1(:,4)==grp(1),:);
                lat_2 = mLat1( mLat1(:,1)==idata & mLat1(:,2)==icmp & mLat1(:,4)==grp(2),:);
                
                subplot(3,2,iloc),
                
                bar([lat_1(1,5), lat_2(1,5); lat_1(2,5), lat_2(2,5) ])
                k = .15;
                hold on
                h = errorbar([1-k, 1+k, 2-k, 2+k], [lat_1(1,5), lat_2(1,5), lat_1(2,5), lat_2(2,5) ], ...
                    [lat_1(1,6),  lat_2(1,6), lat_1(2,6), lat_2(2,6) ],...
                    [lat_1(1,6), lat_2(1,6), lat_1(2,6), lat_2(2,6) ], '.k');
                set(h, 'LineWidth', 3)
                set(gca, 'XTickLabel',  {'Non-Riv','Riv'})
                
                legend(cmpName,'Location', 'SouthOutside')
                ylabel('Latency(s)')
                
                ylim([0 .9])
                title(dataType, 'interpret', 'none')
            end
            
        end
        
    case 5 % control, med-on, med-off, DBS-on, DBS-off patients.
        comGrpName = {'Med-ON vs -OFF','DBS-ON vs -OFF'};
        
        BPAcc_sub(ismember(BPAcc_sub(:,1),oneStateSub),:) = [];
        acc_mon = BPAcc_sub(ismember(BPAcc_sub(:,3),[2]), 4);
        acc_mof = BPAcc_sub(ismember(BPAcc_sub(:,3),[3]), 4);
        acc_don = BPAcc_sub(ismember(BPAcc_sub(:,3),[4]), 4);
        acc_dof = BPAcc_sub(ismember(BPAcc_sub(:,3),[5]), 4);
        
        [H1,P1,CI1,STATS1] = ttest(acc_mon, acc_mof);
        [H2,P2,CI2,STATS2] = ttest(acc_don, acc_dof);
        
        
        j = 0; iloc = 0; istats  = 0;
        figure(1),clf
        for idata = 1:2
            
            switch idata
                case 1
                    D = FBP_sub;
                    dataType = 'FBP';
                case 2
                    D = FCBP_sub;
                    dataType = 'FCBP';
            end
            
            for icmp = 1:2 % we are going to compare 2 sets of subject groups
                iloc = iloc + 1;
                
                switch icmp
                    case 1
                        grp1 = 2;
                        grp2 = 3;
                        grp = [2, 3];
                        cmpName = {'MON', 'MOF'};
                    case 2
                        grp1 = 4;
                        grp2 = 5;
                        grp = [4, 5];
                        cmpName = {'DON', 'DOF'};
                end
                
                %% t-test: Compare subject group
                %lat 1)subID, 2)stim, 3)isubgrp(NOT STATE), 4)nTri, 5)data type,
                % 6)latenct
                for istim = 1:2
                    istats = istats + 1;
                    alLat_1 = D( ismember( D(:,2), istim) & ismember(D(:,3), grp1), 4 );
                    alLat_2 = D( ismember( D(:,2), istim) & ismember(D(:,3), grp2), 4 );
                    
                    [H(istats),P(istats),CI(:,istats),STATS] = ttest(alLat_1, alLat_2) ;
                    tstat(istats) = STATS.tstat;
                    df(istats) = STATS.df;
                    sd(istats) = STATS.sd;
                    statSum(istats,:) = [idata, icmp, istim, P(istats), tstat(istats), df(istats), sd(istats)];
                    
                    j = j + 1;
                    mLat1(j,:) = [ idata, icmp, istim, grp(1), mean(alLat_1), std(alLat_1)/sqrt(length(alLat_1)) ];
                    j = j + 1;
                    mLat1(j,:) = [ idata, icmp, istim, grp(2), mean(alLat_2), std(alLat_2)/sqrt(length(alLat_2)) ];
                end
                %%
                
                lat_1 = mLat1( mLat1(:,1)==idata & mLat1(:,2)==icmp & mLat1(:,4)==grp(1),:);
                lat_2 = mLat1( mLat1(:,1)==idata & mLat1(:,2)==icmp & mLat1(:,4)==grp(2),:);
                
                subplot(3,2,iloc),
                
                bar([lat_1(1,5), lat_2(1,5); lat_1(2,5), lat_2(2,5) ])
                k = .15;
                hold on
                h = errorbar([1-k, 1+k, 2-k, 2+k], [lat_1(1,5), lat_2(1,5), lat_1(2,5), lat_2(2,5) ], ...
                    [lat_1(1,6),  lat_2(1,6), lat_1(2,6), lat_2(2,6) ],...
                    [lat_1(1,6), lat_2(1,6), lat_1(2,6), lat_2(2,6) ], '.k');
                set(h, 'LineWidth', 3)
                set(gca, 'XTickLabel',  {'Non-Riv','Riv'})
                
                legend(cmpName,'Location', 'SouthOutside')
                ylabel('Latency(s)')
                
                ylim([0 .9])
                %         xlim([0.3 2.7])
                title(dataType, 'interpret', 'none')
                
            end
        end
        
    case 2 % Control vs treatment-on patients
        comGrpName = {'Control vs ON treatment PD'};
        
        icmp = 1;
        acc_c   = BPAcc_sub(ismember(BPAcc_sub(:,3),[1]), 4);
        acc_pon = BPAcc_sub(ismember(BPAcc_sub(:,3),[2,4]), 4);
        
        [H1,P1,CI1,STATS1] = ttest2(acc_c, acc_pon);
        
        j = 0; iloc = 0; istats  = 0;
        figure(1),clf
        for idata = 1:2
            
            switch idata
                case 1
                    D = FBP_sub;
                    dataType = 'FBP';
                case 2
                    D = FCBP_sub;
                    dataType = 'FCBP';
            end
            
            iloc = iloc + 1;
            
            grp1 = 1;
            grp2 = [3,5];
            grp  = [1, 35];
            cmpName = {'C', 'PON'};
            
            %% t-test: Compare subject group
            %lat 1)subID, 2)stim, 3)isubgrp(NOT STATE), 4)nTri, 5)data type,
            % 6)latenct
            for istim = 1:2
                istats = istats + 1;
                alLat_1 = D( ismember( D(:,2), istim) & ismember(D(:,3), grp1), 4 );
                alLat_2 = D( ismember( D(:,2), istim) & ismember(D(:,3), grp2), 4 );
                
                [H(istats),P(istats),CI(:,istats),STATS] = ttest2(alLat_1, alLat_2) ;
                tstat(istats) = STATS.tstat;
                df(istats) = STATS.df;
                sd(istats) = STATS.sd;
                statSum(istats,:) = [idata, icmp, istim, P(istats), tstat(istats), df(istats), sd(istats)];
                
                j = j + 1;
                mLat1(j,:) = [ idata, icmp, istim, grp(1), mean(alLat_1), std(alLat_1)/sqrt(length(alLat_1)) ];
                j = j + 1;
                mLat1(j,:) = [ idata, icmp, istim, grp(2), mean(alLat_2), std(alLat_2)/sqrt(length(alLat_2)) ];
            end
            %%
            
            lat_1 = mLat1( mLat1(:,1)==idata & mLat1(:,2)==icmp & mLat1(:,4)==grp(1),:);
            lat_2 = mLat1( mLat1(:,1)==idata & mLat1(:,2)==icmp & mLat1(:,4)==grp(2),:);
            
            subplot(3,2,iloc),
            
            bar([lat_1(1,5), lat_2(1,5); lat_1(2,5), lat_2(2,5) ])
            k = .15;
            hold on
            h = errorbar([1-k, 1+k, 2-k, 2+k], [lat_1(1,5), lat_2(1,5), lat_1(2,5), lat_2(2,5) ], ...
                [lat_1(1,6),  lat_2(1,6), lat_1(2,6), lat_2(2,6) ],...
                [lat_1(1,6), lat_2(1,6), lat_1(2,6), lat_2(2,6) ], '.k');
            set(h, 'LineWidth', 3)
            set(gca, 'XTickLabel',  {'Non-Riv','Riv'})
            
            legend(cmpName,'Location', 'SouthOutside')
            ylabel('Latency(s)')
            
            ylim([0 .9])
            %         xlim([0.3 2.7])
            title(dataType, 'interpret', 'none')
            
        end
              
end


% Insert table for stats results %
f=figure(1);

% Create the column name in cell arrays
cnames = {'Data', 'Compared groups','Stim','p','t-stat','df','sem'};


% Data translation
statData = cell(size(statSum,1), 7);
statData(statSum(:,1)==1,1) = {'FBP'};
statData(statSum(:,1)==2,1) = {'FCBP'};
statData(statSum(:,2)==1,2) = comGrpName(1);
if nSubGrp ~= 2, statData(statSum(:,2)==2,2) = comGrpName(2); end
statData(statSum(:,3)==1,3) = {'Non-Riv'};
statData(statSum(:,3)==2,3) = {'Riv'};
statData(:,4:end) = num2cell(statSum(:,4:end));


% Create the uitable
t = uitable(f,'Data', statData,...
    'ColumnName',cnames,...
    'ColumnWidth',{50, 200, 60, 60, 60, 50, 60});

subplot(3,2,iloc+1:iloc+2),
pos = get(subplot(2,2,3:4),'position');
title('stats results')
set(subplot(3,2,iloc+1:iloc+2),'yTick',[])
set(subplot(3,2,iloc+1:iloc+2),'xTick',[])

set(t,'units','normalized')
set(t,'position',pos)
set(t,'ColumnName',cnames)


set(figure(1),'PaperUnits','inches','PaperPosition',[0 0 8 8])
print('-dpng', [DIR.figIMBP '/IM_BPtime2_nGrp' num2str(nSubGrp) '.png'])


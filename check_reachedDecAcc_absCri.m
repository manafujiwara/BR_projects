function check_reachedDecAcc_absCri(DIR)
% =========================================================================
% Check time of reachde discriminability here. Gray dots above each figure
% represent statistically significant points between groups, while black
% dots indicate such points between data kinds. Choose nSubjectGroup for
% comparison. 2= C vs PD on treatment, 4= WO- vs W-Electrode, 5= ON- vs OFF
% medication and DBS. Color is set as "lineColor" (Ctrl+F please for which
% color is which in each nSubjectGroup setting). Resukts of BR and non-BR
% trials are shown with dotted and solid lines respectively.
% =========================================================================


%%% Chose one of the two options %%%
Plot.reachedAccEachSubThenAverage = 1;
Plot.reachedAccOfAveragedAcc = 0;

%% Prepare data
load( [ DIR.allSub '/subInfo.mat'  ] ,'rejSub', 'accSub', 'subInfo', 'badTri')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nSubGrp          = 4; % 4 or 5. 4)with electrode, without electrode, control, pd. 5) c, mon, mof, don, dof
ntp              = 20;
tRange           = [271 900];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time = linspace(-1, 2, 900);
figure(1), clf, colormap jet
figure(2), clf
figure(4), clf
iLocPlot = 0;

% -------------------------------------------------------------------------
% BPinfo: 1)subID, 2)irun, 3)iblk, 4)iepc, 5)stimLabel,
% 6)subjestState, 7)BPP (-1~1), 8)label, 9)fBP, 10)fcBP, 11:910)BP(t)
% -------------------------------------------------------------------------
% eyeinfo: 1)subID, 2)irun, 3)iblk, 4)iepc, 5)stimLabel, 6)subjestState
% 7)BPP (-1~1), 8)label, 9)fBP, 10)fcBP, 11:910)eyeinfo(t)
% -------------------------------------------------------------------------
% subInfo  1)subNum, 2)nallTri, 3)length(badTri(isub).BPeye30),
% 4)length(badTri(isub).BPeye150), 5)rejBPeye30, 6)val_br_rBP, 7)val_br_lBP;
% -------------------------------------------------------------------------

load([DIR.SVM '/IM_discriminability_nGrp5_tp' num2str(ntp) '_BP_balTri.mat'])
load([DIR.SVM '/IM_discriminability_nGrp5_tp' num2str(ntp) '_OKN30_balTri.mat'])

for istim = 1:2 % non-rivalry or br
    
    switch istim
        case 1 % stim = non-br
            titleStim = '_non';
            
        case 2 % stim = br
            titleStim = '_br';
    end
    
    time4reachedAccStore         = [];
    time4reachedAccStoreSubGrp   = [];
    time4reachedAccStoreDataType = [];
    
    for idata = 1:2 % BP or eye
        
        switch idata
            case 1
                D = BP;
                lineStyle = '-';
                dataType = 'BP';
            case 2
                D = OKN30;
                lineStyle = '--';
                dataType = 'OKN';
        end
        
        %% Compute mean
        %%% mean across trials then subject : for shadedError
        
        for isubgrp = 1:nSubGrp % control, medOn medOff dbsOn dbsOff
            
            switch nSubGrp
                case 5
                    if isubgrp == 1
                        continue
                    end
                    nSubYLim = [0 18];
                    
                    switch istim
                        case 1
                            oriYLim = 1.2;
                            adjYLim = 1;
                        case 2
                            oriYLim = 1.7;
                            adjYLim = 1.5;
                    end
                    
                    subgrp = isubgrp;
                    switch isubgrp
                        case 1
                            lineColor = [0 0 .9];
                            cmpName   = 'control';
                        case 2
                            nSubYLim = [0 18];
                            iLocPlot  = iLocPlot + 1;
                            lineColor = [0 .8 0];
                            cmpName   = 'mon-mof';
                        case 3
                            nSubYLim = [0 18];
                            lineColor = [0 .2 0];
                            cmpName   = 'mon-mof';
                        case 4
                            nSubYLim = [0 8];
                            iLocPlot  = iLocPlot + 1;
                            lineColor = [0 0 1];
                            cmpName   = 'don-dof';
                        case 5
                            nSubYLim = [0 8];
                            lineColor = [0 0 .2];
                            cmpName   = 'don-dof';
                    end
                    
                case 4
                    nSubYLim = [0 32];
                    oriYLim = 1.4;
                    adjYLim = 1.2;
                    switch isubgrp
                        case 1 % without electrode
                            iLocPlot  = iLocPlot + 1;
                            subgrp = [2,3];
                            lineColor = [0 .5 0];
                            cmpName = 'w-woElectrode';
                        case 2 % with electrode
                            subgrp = [4,5];
                            lineColor = [0 0 .6];
                            cmpName = 'w-woElectrode';
                        case 3 % control
                            iLocPlot  = iLocPlot + 1;
                            subgrp = 1;
                            lineColor = [0 0 .9];
                            cmpName = 'c-pd';
                        case 4 % PD patients
                            subgrp = [2,3,4,5];
                            lineColor = [.8 0 0];
                            cmpName = 'c-pd';
                    end
                    
                case 2
                    nSubYLim = [0 25];
                    oriYLim = 1.2;
                    adjYLim = 1;
                    switch isubgrp
                        case 1 % control
                            iLocPlot  = iLocPlot + 1;
                            subgrp = 1;
                            lineColor = [0 0 .9];
                            cmpName = 'c-pd';
                        case 2 % PD patients
                            subgrp = [3,5];
                            lineColor = [.8 0 0];
                            cmpName = 'c-pd';
                    end
            end
            
            %BP/OKN 1)subID, 2)stim, 3)isubgrp(NOT STATE), 4)nTri, 5)data
            % type, 6:905)mTriAmplitude
            
            tmpD = D( ismember(D(:,2),istim) & ismember(D(:,3), subgrp ) & D(:,5)==idata,:);
            nSub = size(tmpD,1);
            
            if Plot.reachedAccEachSubThenAverage
                reachedAcc = tmpD(:,306);
                for isamp = 2:600
                    reachedAcc(:,isamp) = max( [ reachedAcc(:,isamp-1), tmpD(:,305+isamp) ], [], 2);
                end
                
                reachedAccLinSpace = (.5:.001:1)';
                
                for isub = 1:nSub
                    
                    uniqAcc = unique(reachedAcc(isub,:));
                    nAcc = length(uniqAcc);
                    time_acc = [ time(301:900) ; reachedAcc(isub,:) ]';
                    
                    for iAcc = 1 : nAcc-1
                        timeBetween = time_acc( time_acc(:,2) >= uniqAcc(iAcc) & time_acc(:,2) < uniqAcc(iAcc+1), : );
                        reachedAccLinSpace( reachedAccLinSpace(:,1) >= uniqAcc(iAcc) & reachedAccLinSpace(:,1) < uniqAcc(iAcc+1), 1+isub ) = max(timeBetween(:,1));
                    end
                    
                    reachedAccLinSpace( reachedAccLinSpace(:,1) == uniqAcc(nAcc), 1+isub ) = min( time_acc( time_acc(:,2) == uniqAcc(nAcc), 1 ) );
                    reachedAccLinSpace( reachedAccLinSpace(:,1) > uniqAcc(nAcc), 1+isub ) = nan; % accuracy doesn't reach the level
                    
                    
                end
                
                time4reachedAcc = reachedAccLinSpace(:,2:end);
                mtime4reachedAcc = nanmean(reachedAccLinSpace(:,2:end), 2);
                semtime4reachedAcc = nanstd(reachedAccLinSpace(:,2:end), [], 2) / sqrt(nSub);
                
                %%% Plot reached accuracy on x axis
                time4reachedAccStore         = cat( 1, time4reachedAccStore, time4reachedAcc' ); % Collect data for t-test
                time4reachedAccStoreSubGrp   = cat( 1, time4reachedAccStoreSubGrp, isubgrp * ones( size(time4reachedAcc,2), 1 ) );
                time4reachedAccStoreDataType = cat( 1, time4reachedAccStoreDataType, idata * ones( size(time4reachedAcc,2), 1 ) );
                
                nSunSamp( isubgrp, : ) = sum( ~isnan( time4reachedAcc ), 2 )';
                
                
                %%% Plot average %%%
                figure(2),
                switch nSubGrp
                    case 5
                        subplot( 5, 4, iLocPlot + 2*(istim-1) )
                    case 4
                        subplot( 5, 2, iLocPlot + 2*(istim-1) )
                    case 2
                        subplot( 5, 2, iLocPlot * istim )
                end
                hh = shadedErrorBar( reachedAccLinSpace(:,1), mtime4reachedAcc, semtime4reachedAcc, [], 1 );
                
                ylim([0 oriYLim])
                xlim([.5 1])
                xlabel('Discriminability')
                ylabel('Time (sec)')
                
                hold on
                set(hh.mainLine, 'LineWidth', 2, 'LineStyle', lineStyle , 'Color', lineColor)
                set(hh.edge,'Color', 'none')
                set(hh.patch,'FaceColor', lineColor)
                
                %%% Plot nSubject at each discriminability level point %%%
                figure(4),
                switch nSubGrp
                    case 5
                        subplot( 5, 4, iLocPlot + 2*(istim-1) )
                    case 4
                        subplot( 5, 2, iLocPlot + 2*(istim-1) )
                    case 2
                        subplot( 5, 2, iLocPlot * istim )
                end
                
                plot( reachedAccLinSpace(:,1), nSunSamp( isubgrp, : ) ,...
                    'LineStyle', lineStyle, 'Color', lineColor, 'LineWidth', 2);
                ylim( nSubYLim )
                xlim([.5 1])
                xlabel('Discriminability reached')
                ylabel('n subjects')
                hold on
            end
            
        end
        iLocPlot = 0;
    end
    
    if Plot.reachedAccEachSubThenAverage
        
        ySig = 1.5;
        
        switch nSubGrp
            case 2 % ANOVA
                for isamp = 1:501
                    statsData     = time4reachedAccStore(:,isamp);
                    statsSubGrp   = time4reachedAccStoreSubGrp(~isnan(statsData));
                    statsDataType = time4reachedAccStoreDataType(~isnan(statsData));
                    statsData     = statsData(~isnan(statsData));
                    
                    [ p(:,isamp), ~, stats(:,isamp) ] = anovan( statsData,...
                        { statsSubGrp, statsDataType },'display','off');
                end
                %         arrayfun( @(n) fdr( p( :, n ) ), 1:size(p,2) )
                
                [pthr,pcor,padj_subGrp] = fdr( p(1,:), 0.05 ) ;
                [pthr,pcor,padj_dataType] = fdr( p(2,:), 0.05 ) ;
                
                sig_subGrp = zeros(1,length(p));
                sig_subGrp( padj_subGrp < 0.05) = 1;
                sig_subGrp( padj_subGrp >= 0.05 ) = nan;
                
                sig_dataType = zeros(1,length(p));
                sig_dataType( padj_dataType < 0.05) = 1;
                sig_dataType( padj_dataType >= 0.05 ) = nan;
                
                figure(2),
                a_subGrp = reachedAccLinSpace(:,1)' .* sig_subGrp;
                a_dataType = reachedAccLinSpace(:,1)' .* sig_dataType;
                plot( a_subGrp, oriYLim, '.k', 'MarkerSize', 10, 'Color', [.5 .5 .5]);
                plot( a_dataType, 1.1, '.k', 'MarkerSize', 10 );
                
                %%% TO PLOT OUTSIDE OF THE BOX %%%
                %# create a second axis as copy of first (without its content),
                %# reduce its size, and set limits accordingly
                hAx1 = gca;
                hAx2 = copyobj(hAx1,gcf);
                set(hAx2, 'Position', get(hAx1,'Position').*[1 1 1 1/oriYLim], ...
                    'XLimMode', 'manual', 'YLimMode', 'manual', ...
                    'YLim', get(hAx1,'YLim').*[1 1/oriYLim])
                delete(get(hAx2,'Children'))
                
                %# hide first axis, and adjust Z-order
                axis(hAx1,'off')
                uistack(hAx1,'top')
                
                
            case 4
                
                for icmp = 1:2
                    
                    switch icmp
                        case 1 % compare med on vs off
                            state = [ 1, 2 ];
                        case 2 % compare dbs on vs off
                            state = [ 3, 4 ];
                    end
                    
                    for isamp = 1:501
                        time4reachedAccStoreCmp         = time4reachedAccStore( ismember(time4reachedAccStoreSubGrp, state ), : );
                        time4reachedAccStoreSubGrpCmp   = time4reachedAccStoreSubGrp( ismember(time4reachedAccStoreSubGrp, state ), : );
                        time4reachedAccStoreDataTypeCmp = time4reachedAccStoreDataType( ismember(time4reachedAccStoreSubGrp, state ), : );
                        
                        statsData     = time4reachedAccStoreCmp( :, isamp );
                        statsSubGrp   = time4reachedAccStoreSubGrpCmp( ~isnan(statsData) );
                        statsDataType = time4reachedAccStoreDataTypeCmp( ~isnan(statsData) );
                        statsData     = statsData( ~isnan(statsData) );
                        
                        [ p(:,isamp), ~, stats(:,isamp) ] = anovan( statsData,...
                            { statsSubGrp, statsDataType }, 'display', 'off');
                    end
                    
                    %         arrayfun( @(n) fdr( p( :, n ) ), 1:size(p,2) )
                    
                    [pthr,pcor,padj_subGrp]   = fdr( p(1,:), 0.05 ) ;
                    [pthr,pcor,padj_dataType] = fdr( p(2,:), 0.05 ) ;
                    
                    sig_subGrp = zeros(1,length(p));
                    sig_subGrp( padj_subGrp < 0.05) = 1;
                    sig_subGrp( padj_subGrp >= 0.05 ) = nan;
                    
                    sig_dataType = zeros(1,length(p));
                    sig_dataType( padj_dataType < 0.05) = 1;
                    sig_dataType( padj_dataType >= 0.05 ) = nan;
                    
                    figure(2),
                    subplot( 5, 2, icmp + 2*(istim-1)  )
                    hold on
                    a_subGrp = reachedAccLinSpace(:,1)' .* sig_subGrp;
                    a_dataType = reachedAccLinSpace(:,1)' .* sig_dataType;
                    plot( a_subGrp, oriYLim, '.k', 'MarkerSize', 10, 'Color', [.5 .5 .5]);
                    plot( a_dataType, oriYLim - .1, '.k', 'MarkerSize', 10 );
                    
                    %%% TO PLOT OUTSIDE OF THE BOX %%%
                    %# create a second axis as copy of first (without its content),
                    %# reduce its size, and set limits accordingly
                    hAx1 = gca;
                    hAx2 = copyobj(hAx1,gcf);
                    set(hAx2, 'Position', get(hAx1,'Position').*[1 1 1 adjYLim/oriYLim], ...
                        'XLimMode', 'manual', 'YLimMode', 'manual', ...
                        'YLim', get(hAx1,'YLim').*[1 adjYLim/oriYLim])
                    delete(get(hAx2,'Children'))
                    
                    %# hide first axis, and adjust Z-order
                    axis(hAx1,'off')
                    uistack(hAx1,'top')
                    
                end
            case 5 % PAIRED t-test
                
                for icmp = 1:2
                    
                    switch icmp
                        case 1 % compare med on vs off
                            state = [ 2, 3 ];
                        case 2 % compare dbs on vs off
                            state = [ 4, 5 ];
                    end
                    
                    for isamp = 1:501
                        time4reachedAccStoreCmp         = time4reachedAccStore( ismember(time4reachedAccStoreSubGrp, state ), : );
                        time4reachedAccStoreSubGrpCmp   = time4reachedAccStoreSubGrp( ismember(time4reachedAccStoreSubGrp, state ), : );
                        time4reachedAccStoreDataTypeCmp = time4reachedAccStoreDataType( ismember(time4reachedAccStoreSubGrp, state ), : );
                        
                        statsData     = time4reachedAccStoreCmp( :, isamp );
                        statsSubGrp   = time4reachedAccStoreSubGrpCmp( ~isnan(statsData) );
                        statsDataType = time4reachedAccStoreDataTypeCmp( ~isnan(statsData) );
                        statsData     = statsData( ~isnan(statsData) );
                        
                        [ p(:,isamp), ~, stats(:,isamp) ] = anovan( statsData,...
                            { statsSubGrp, statsDataType }, 'display', 'off');
                    end
                    
                    [pthr,pcor,padj_subGrp]   = fdr( p(1,:), 0.05 ) ;
                    [pthr,pcor,padj_dataType] = fdr( p(2,:), 0.05 ) ;
                    
                    sig_subGrp = zeros(1,length(p));
                    sig_subGrp( padj_subGrp < 0.05) = 1;
                    sig_subGrp( padj_subGrp >= 0.05 ) = nan;
                    
                    sig_dataType = zeros(1,length(p));
                    sig_dataType( padj_dataType < 0.05) = 1;
                    sig_dataType( padj_dataType >= 0.05 ) = nan;
                    
                    figure(2),
                    subplot( 5, 4, icmp + 2*(istim-1)  )
                    hold on
                    a_subGrp = reachedAccLinSpace(:,1)' .* sig_subGrp;
                    a_dataType = reachedAccLinSpace(:,1)' .* sig_dataType;
                    plot( a_subGrp, oriYLim, '.k', 'MarkerSize', 10, 'Color', [.5 .5 .5]);
                    plot( a_dataType, oriYLim - .1, '.k', 'MarkerSize', 10 );
                    
                    %%% TO PLOT OUTSIDE OF THE BOX %%%
                    %# create a second axis as copy of first (without its content),
                    %# reduce its size, and set limits accordingly
                    hAx1 = gca;
                    hAx2 = copyobj(hAx1,gcf);
                    set(hAx2, 'Position', get(hAx1,'Position').*[1 1 1 adjYLim/oriYLim], ...
                        'XLimMode', 'manual', 'YLimMode', 'manual', ...
                        'YLim', get(hAx1,'YLim').*[1 adjYLim/oriYLim])
                    delete(get(hAx2,'Children'))
                    
                    %# hide first axis, and adjust Z-order
                    axis(hAx1,'off')
                    uistack(hAx1,'top')
                    
                end
        end
    end
end


switch nSubGrp
    case 5
        set(figure(2),'PaperUnits','inches','PaperPosition',[0 0 14 9])
        set(figure(4),'PaperUnits','inches','PaperPosition',[0 0 14 9])
    case {2,4}
        set(figure(2),'PaperUnits','inches','PaperPosition',[0 0 7 9])
        set(figure(4),'PaperUnits','inches','PaperPosition',[0 0 7 9])
end


if Plot.reachedAccEachSubThenAverage
    figure(2)
    print('-dpng', [DIR.figRevision '/IM_svm_absReachedAcc_nGrp' num2str(nSubGrp) '_1_balTri_farDots.png'] )
elseif Plot.reachedAccOfAveragedAcc
    figure(2)
    print('-dpng', [DIR.figRevision '/IM_svm_absReachedAcc_nGrp' num2str(nSubGrp) '_2.png'] )
end

if Plot.reachedAccEachSubThenAverage
    figure(4)
    print('-dpng', [DIR.figRevision '/IM_svm_nSub_absReachedAcc_nGrp' num2str(nSubGrp) '_1_balTri_farDots.png'] )
elseif Plot.reachedAccOfAveragedAcc
    figure(4)
    print('-dpng', [DIR.figRevision '/IM_svm_nSub_absReachedAcc_nGrp' num2str(nSubGrp) '_2.png'] )
end




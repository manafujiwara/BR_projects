function check_cmp_subGrp(DIR)
%% Changed from check_cmp_inR
% Eventually no need to run R, so changed the name and some more

%% Prepare data
load( [ DIR.allSub '/subInfo.mat'  ] ,'rejSub', 'accSub', 'subInfo', 'badTri')
% load( [ DIR.allSub '/allSub_IM_rej_upsamp_smth_sort_BP_a_slwPhs_lb.mat']) % smoothed data for imagesc, but not for plot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nSubGrp          = 2; % 4 or 5. 4)with electrode, without electrode, control, pd. 5) c, mon, mof, don, dof
discriminability = 1;
triBalance       = 1;
ntp              = 20;
tRange           = [271 900];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

time = linspace(-1, 2, 900);
figure(1), clf, colormap jet
figure(2), clf
iLocPlot = 0;
iLocImage = 0;
jsubGrp = 0;
yLimDot = [ -.4 1.2];
jsub = 0;
lat = [];
j = 0;
maxOfMean = [];
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

switch discriminability
    case 0
        load([DIR.allSub '/IM_meanAmp_latency_nGrp5_latTp' num2str(ntp) '.mat'])
    case 1
        if triBalance
            load([DIR.SVM '/IM_discriminability_nGrp5_tp' num2str(ntp) '_BP_balTri.mat'])
            load([DIR.SVM '/IM_discriminability_nGrp5_tp' num2str(ntp) '_OKN30_balTri.mat'])
        else
            load([DIR.SVM '/IM_discriminability_nGrp5_tp' num2str(ntp) '_BP.mat'])
            load([DIR.SVM '/IM_discriminability_nGrp5_tp' num2str(ntp) '_OKN30.mat'])
        end
end

for idata = 1:2 % BP or eye gaze
    
    switch idata
        case 1
            D = BP;
            dataType  = 'BP';
            titleData = 'BP'; 
        case 2
            D = OKN30;
            dataType = 'velSmooth100ms';
            titleData = 'OKN30';
    end
    
    for istim = 1:2 % non-rivalry or br
        
        switch istim
            case 1 % stim = non-br
                titleStim = '_nonBR';
            case 2 % stim = br
                titleStim = '_BR';
        end
        
        if discriminability
            yLim = [ .3 1.1];
            
        else
            switch idata
                case 1 % BP
                    yLim = [ -.4 1.2];
                    
                case 2 % eye gaze
                    switch istim
                        case 1 % stim = non-br
                            yLim = [-2 5];
                        case 2 % stim = br
                            yLim = [-1 2];
                    end
            end
        end
        
        switch discriminability
            case 0
                switch istim
                    case 1
                        stim = [2,3];
                    case 2
                        stim = 1;
                end
            case 1
                stim = istim;
        end
        
        
        %% Compute mean
        %%% mean across trials then subject : for shadedError

        
        for isubgrp = 1:nSubGrp % control, medOn medOff dbsOn dbsOff
            
            iLocImage = iLocImage + 1;
            jsubGrp = jsubGrp + 1;
            switch nSubGrp
                case 5
                    comGrpName = {'Med-ON vs -OFF','DBS-ON vs -OFF'};
                    subgrp = isubgrp;
                    switch isubgrp
                        case 1
                            iLocPlot = iLocPlot+ 1;
                            lineStyle = '-';
                            lineColor = [0 0 .9];
                            cmpName = 'control';
                        case 2
                            iLocPlot = iLocPlot+ 1;
                            lineStyle = '-';
                            lineColor = [0 .8 0];
                            cmpName = 'mon-mof';
                        case 3
                            lineStyle = '-';
                            lineColor = [0 .2 0];
                            cmpName = 'mon-mof';
                        case 4
                            iLocPlot = iLocPlot+ 1;
                            lineStyle = '-';
                            lineColor = [0 0 1];
                            cmpName = 'don-dof';
                        case 5
                            lineStyle = '-';
                            lineColor = [0 0 .2];
                            cmpName = 'don-dof';
                    end
                    
                case 4
                    comGrpName = {'WO-Electrode vs W-','Conreol vs All PD'};
                    
                    switch isubgrp
                        case 1 % without electrode
                            subgrp = [2,3];
                            iLocPlot = iLocPlot+ 1;
                            lineStyle = '-';
                            lineColor = [0 .5 0];
                            cmpName = 'w-woElectrode';
                        case 2 % with electrode
                            subgrp = [4,5];
                            lineStyle = '-';
                            lineColor = [0 0 .6];
                            cmpName = 'w-woElectrode';
                        case 3 % control
                            subgrp = 1;
                            iLocPlot = iLocPlot+ 1;
                            lineStyle = '-';
                            lineColor = [0 0 .9];
                            cmpName = 'c-pd';
                        case 4 % PD patients
                            subgrp = [2,3,4,5];
                            lineStyle = '-';
                            lineColor = [.8 0 0];
                            cmpName = 'c-pd';
                    end
                    
                case 2
                    comGrpName = {'Control vs ON treatment PD'};
                    switch isubgrp
                        case 1 % control
                            subgrp = 1;
                            iLocPlot = iLocPlot+ 1;
                            lineStyle = '-';
                            lineColor = [0 0 .9];
                            cmpName = 'c-pon';
                        case 2 % PD patients on treatment
                            subgrp = [3,5];
                            lineStyle = '-';
                            lineColor = [.8 0 0];
                            cmpName = 'c-pon';
                    end
            end
            
            %BP/OKN 1)subID, 2)stim, 3)state, 4)nTri, 5)data type, 6:905)mTriAmplitude
            
            tmpD = D( ismember( D(:,2), stim ) & ismember( D(:,3), subgrp ) & D(:,5)==idata, : );
            nSub = size(tmpD,1);
            
            mTrimSub_DD = nanmean(tmpD(:,6:905),1);
            semSub_DD   = nanstd(tmpD(:,6:905),1)/sqrt(nSub);
            
            tmp3(isubgrp).d = tmpD; % Collect data for t-test
            
            %% Compute half maximum
            maxVel = max(tmpD(:,306:905), [], 2);
            lat = [ lat;  tmpD(:,1:5), nan(nSub,1) ];
            
            switch discriminability
                case 0
                    threshold = mean( [zeros(nSub,1), maxVel], 2 );
                case 1
                    threshold = mean( [ones(nSub,1)*0.5, maxVel], 2 );
            end
            
            for isub = 1:nSub
                jsub = jsub + 1;
                lat(jsub, 6) = time( 300 + min( find( tmpD( isub, 306:905 ) >= threshold (isub) ) ) );
            end
            
            mLat2(jsubGrp,:) = [ idata, istim, isubgrp, nanmean(lat(end-nSub+1:end,6)), std(lat(end-nSub+1:end,6))/sqrt(nSub) ];
            
            %% Let's plot
            figure(2), % SHADEDERRORBAR ===================================
            switch nSubGrp
                case 5
                    subplot(5,3,iLocPlot)
                case {2,4}
                    subplot(5,2,iLocPlot)
            end
           
            
            % PLOT TIME COURSE ===================================
            hh = shadedErrorBar( time( tRange(1):tRange(2) ), mTrimSub_DD( tRange(1):tRange(2) ), semSub_DD( tRange(1):tRange(2) ), [], 1);
            hold on
            set(hh.mainLine, 'LineWidth', 2, 'LineStyle', lineStyle , 'Color', lineColor)
            set(hh.edge,'Color', 'none')
            set(hh.patch,'FaceColor', lineColor)

            maxOfMean(iLocPlot,isubgrp) = max(mTrimSub_DD(300:900));
            
            % PLOT LATENCY ===================================
            hh = area([ mLat2(jsubGrp,4) - mLat2(jsubGrp,5), mLat2(jsubGrp,4) + mLat2(jsubGrp,5) ], [yLim(2), yLim(2)], yLim(1));
            set(hh, 'LineWidth', 2, 'LineStyle', lineStyle , 'FaceColor', lineColor, 'EdgeColor', 'none')
            alpha(.3)
            
            hh = plot([ mLat2(jsubGrp,4), mLat2(jsubGrp,4)], yLim);
            set(hh, 'Color', lineColor, 'LineWidth', 1.5)
            
            plot([ 0 0 ], [yLim], 'k--');
            
            ylim(yLim)
            xlim( [time(tRange(1)),time(tRange(2))] )
            titlePlot = cat(2, titleData, titleStim);
            title(titlePlot, 'interpret', 'none')
            
            switch discriminability
                case 0
                    plot([time(tRange(1)), time(tRange(2))], [ 0 0 ], 'k--');
                    ySig = yLim(1) + diff(yLim)*0.0625;
                case 1
                    plot([time(tRange(1)), time(tRange(2))], [ .5 .5 ], 'k--');
                    ySig = 0.35;
            end
            
            if [nSubGrp == 4 && ismember(isubgrp, [2,4])] || [nSubGrp == 5 && ismember(isubgrp, [3,5])] ||...
                   [nSubGrp == 2 && ismember(isubgrp, 2)] 
                j = j + 1;
                
                switch nSubGrp
                    case {4,2} % UNpaired t-test
                        [H, P, CI, STATS] = ttest2( tmp3(isubgrp-1).d(:,6:905), tmp3(isubgrp).d(:,6:905) );
                        
                    case 5 % PAIRED t-test
                        [H, P, CI, STATS] = ttest( tmp3(isubgrp-1).d(:,6:905), tmp3(isubgrp).d(:,6:905) );
                end
                [RHO, P_corr] = ttest2( tmp3(isubgrp-1).d(:,6:905), tmp3(isubgrp).d(:,6:905) );

                
                [pthr,pcor,padj] = fdr(P,0.05); 
                sig = zeros(1,length(P));
                sig(padj < 0.05) = 1;
                sig( padj >= 0.05 ) = nan;
                
                a = time .* sig;
                h = plot( a( tRange(1):tRange(2) ), ySig , '.k');
                set(h,'MarkerSize', 10)
                
                if discriminability
                    %%% max discriminability %%%
                    
                    maxDisc1 = max( tmp3(isubgrp-1).d(:,6:905) , [], 2 );
                    maxDisc2 = max( tmp3(isubgrp).d(:,6:905) , [], 2 );
                    
                    [H2, P2, CI, STATS2] = ttest2( maxDisc1, maxDisc2 );
                    
                    cmpMaxAmp(j,:) = [idata, istim, mean( maxDisc1 ), mean( maxDisc2 ),...
                        P2, STATS2.tstat, STATS2.df, STATS2.sd];
                    
                    
%                     cmpMaxAmp(j,:) = [idata, istim, mean( maxDisc1 ), mean( maxDisc2 ),...
%                         std( maxDisc1 )/sqrt(length( maxDisc1 )), std( maxDisc2 )/sqrt(length( maxDisc2 )),...
%                         P2, STATS2.tstat, STATS2.df, STATS2.sd];
                    
                    
%                     maxAmp2 (2*j-1, :) = [idata, istim, isubgrp-1, ...
%                         mean(maxAmpSub(2*j-1).d(:,4)), std(maxAmpSub(2*j-1).d(:,4))/sqrt(length(maxAmpSub(2*j-1).d(:,4))), ...
%                         mean(maxAmpSub(2*j-1).d(:,5)), std(maxAmpSub(2*j-1).d(:,5))/sqrt(length(maxAmpSub(2*j-1).d(:,5))) ];
%                     
%                     maxAmp2 (2*j, :)   = [idata, istim, isubgrp,...
%                         mean(maxAmpSub(2*j).d(:,4)), std(maxAmpSub(2*j).d(:,4))/sqrt(length(maxAmpSub(2*j).d(:,4))) ...
%                         mean(maxAmpSub(2*j).d(:,5)), std(maxAmpSub(2*j).d(:,5))/sqrt(length(maxAmpSub(2*j).d(:,5)))];
                    
%                     %%% 75% discriminability %%%
%                     [H3, P3, CI, STATS3] = ttest2( maxAmpSub(2*j-1).pct75(:,5), maxAmpSub(2*j).pct75(:,5));
%                     
%                     cmpAmp75(j,:) = [idata, istim, isubgrp, H3, P3, STATS3.tstat, STATS3.df, STATS3.sd];
%                     
%                     amp75 (2*j-1, :) = [idata, istim, isubgrp-1, ...
%                         nanmean(maxAmpSub(2*j-1).pct75(:,4)), nanstd(maxAmpSub(2*j-1).pct75(:,4))/sqrt(length(maxAmpSub(2*j-1).pct75(:,4))), ...
%                         nanmean(maxAmpSub(2*j-1).pct75(:,5)), nanstd(maxAmpSub(2*j-1).pct75(:,5))/sqrt(length(maxAmpSub(2*j-1).pct75(:,5))) ];
%                     
%                     amp75 (2*j, :)   = [idata, istim, isubgrp,...
%                         nanmean(maxAmpSub(2*j).pct75(:,4)), nanstd(maxAmpSub(2*j).pct75(:,4))/sqrt(length(maxAmpSub(2*j).pct75(:,4))) ...
%                         nanmean(maxAmpSub(2*j).pct75(:,5)), nanstd(maxAmpSub(2*j).pct75(:,5))/sqrt(length(maxAmpSub(2*j).pct75(:,5)))];
                    
                end
            end
            
            
        end
    end
end

statSum = cmpMaxAmp;

% Insert a table for stats results %
f=figure(2);
iLocPlot = iLocPlot+1;

% Create the column name in cell arrays
cnames = {'dataTrained', 'Compared groups', 'stim', 'maxDiscriminability1', ...
    'maxDiscriminability2', 'p','t-stat','df','sd'};



% Data translation
statData = cell(size(statSum,1), length(cnames));
statData(statSum(:,1)==1,1) = {'BP'};
statData(statSum(:,1)==2,1) = {'eyeGaze'};
statData(:,2) = comGrpName(1);
statData(statSum(:,2)==1,3) = {'Non-Riv'};
statData(statSum(:,2)==2,3) = {'Riv'};
statData(:,4:end) = num2cell(statSum(:,3:end));

% Create the uitable
t = uitable(f,'Data', statData,...
    'ColumnName',cnames,...
    'ColumnWidth',{60, 170, 60, 100, 100, 50, 60, 30, 60});

switch nSubGrp
    case 5
        pos = get( subplot(5,3,iLocPlot:iLocPlot+1),'position');
        set( subplot(5,3,iLocPlot:iLocPlot+1),'yTick',[])
        set( subplot(5,3,iLocPlot:iLocPlot+1),'xTick',[])
    case {2,4}
        pos = get( subplot(5,2,iLocPlot:iLocPlot+1),'position');
        set( subplot(5,2,iLocPlot:iLocPlot+1),'yTick',[])
        set( subplot(5,2,iLocPlot:iLocPlot+1),'xTick',[])
end

title('stats results')


set(t,'units','normalized')
set(t,'position',pos)
set(t,'ColumnName',cnames)


fid = fopen( [ DIR.figRevision 'maxDiscStats_betwnGrp_nGrp'  num2str(nSubGrp) '.csv' ], 'wt' );
fprintf( fid, '%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n', ...
    'data', 'istim', 'meanGrp1', 'meanGrp2', 'semGrp1', 'semGrp2', 'significance','p', 'tstat', 'df', 'sd');
fprintf( fid, '%d,%d,%f,%f,%f,%f,%d,%f,%f,%f,%f\n', cmpMaxAmp' );
fclose(fid);

% if discriminability
%     maxAmp2
%     cmpMaxAmp
%     amp75
%     cmpAmp75
% end

switch nSubGrp
    case 5
        set(figure(2),'PaperUnits','inches','PaperPosition',[0 0 13 9])
    case {4,2}
        set(figure(2),'PaperUnits','inches','PaperPosition',[0 0 13 9])
end

switch discriminability
    case 0
        print('-dpng', [DIR.figRevision '/IM_dynamics_nGrp' num2str(nSubGrp) '_vertBar.png'] )
    case 1
        if triBalance
            print('-dpng', [DIR.figRevision '/IM_svm_disc_nGrp' num2str(nSubGrp) '_10crsVal_balTri_vertBar.png'] )
        else
            print('-dpng', [DIR.figRevision '/IM_svm_disc_nGrp' num2str(nSubGrp) '_10crsVal.png'] )
        end
end
clear statSum

%% Plot mean latencies

% mLat -----------------------------------------------------
% 1)idata, 2)istim, 3)isubGroup, 4)mean:50% of max, 5)SEM  |
% ----------------------------------------------------------
figure(1),clf
iloc = 0;
j = 0;
istats = 0;
clear H P CI STATS
for idata = 1:2
    switch idata
        case 1
            dataType  = 'BP';
        case 2
            dataType = 'OKN';
    end
    
    switch nSubGrp
        case {4, 5}
                nCmp = 2;
        case 2
                nCmp = 1;
    end

    for icmp = 1:nCmp % we are going to compare 2 sets of subject groups
        iloc = iloc + 1;
        
        switch icmp
            case 1
                switch nSubGrp
                    case 4
                        grp1 = [2,3];
                        grp2 = [4,5];
                        cmpName = {'woEle', 'wEle'};
                        barColor1 = [0 .5 0];
                        barColor2 = [0 0 .6];
                    case 5
                        grp1 = 2;
                        grp2 = 3; %MAYBE [1, 2]
                        cmpName = {'MON', 'MOF'};
                        barColor1 = [0 .8 0];
                        barColor2 = [0 .2 0];
                    case 2
                        grp1 = 1; %MAYBE [1, 2]
                        grp2 = [3,5];
                        cmpName = {'C', 'PON'};
                        barColor1 = [0 0 .9];
                        barColor2 = [.8 0 0];
                end
                
            case 2
                switch nSubGrp
                    case 4
                        grp1 = 1;
                        grp2 = [2,3,4,5];
                        cmpName = {'Control', 'PD'};
                        barColor1 = [0 0 .9];
                        barColor2 = [.8 0 0];
                        
                    case 5
                        grp1 = 4;
                        grp2 = 5;
                        cmpName = {'DON', 'DOF'};
                        barColor1 = [0 0 1];
                        barColor2 = [0 0 .2];
                end
        end
        
        if idata == 1
            %% 3WAY-ANOVA
            % lat: 1)subID, 2)stim, 3)state, 4)nTri, 5)data type,
            % 6)latency
            cmp_lat = lat( ismember( lat(:,3), [grp1,grp2]), : );
            tmp = cmp_lat;
            tmp( ismember( cmp_lat(:,3), grp1 ), 3 ) = 1;
            tmp( ismember( cmp_lat(:,3), grp2 ), 3 ) = 2;
            cmp_lat(:,3) = tmp(:,3);
            
            Y = cmp_lat(:, 6);
            [P_anova,T,STATS_anova,TERMS] = anovan(Y, { cmp_lat(:,2), cmp_lat(:,3), cmp_lat(:,5) },...
                'model', 3, 'varnames', {'Stim','SubGrp','DataType'} ); %Stim, SubGrp, dataType
            
            anovaSum(icmp).T = T;
        end
        
        %% t-test: Compare subject group %%
        %lat 1)subID, 2)stim, 3)isubgrp(NOT STATE), 4)nTri, 5)data type,
        % 6)latency
        for istim = 1:2
            
            switch discriminability
                case 0
                    switch istim
                        case 1
                            stim = [2,3];
                        case 2
                            stim = 1;
                    end
                case 1
                    stim = istim;
            end
            
            istats = istats + 1;
            alLat_1 = lat( ismember( lat(:,2), stim ) & lat(:,5) == idata & ismember( lat(:,3), grp1 ), 6 );
            alLat_2 = lat( ismember( lat(:,2), stim ) & lat(:,5) == idata & ismember( lat(:,3), grp2 ), 6 );
            
          
            switch nSubGrp
                case {4,2}
                    [H(istats),P(istats),CI(istats,:),STATS] = ttest2(alLat_1, alLat_2) ;
                case 5
                    [H(istats),P(istats),CI(istats,:),STATS] = ttest(alLat_1, alLat_2) ;
            end
            
            tstat(istats) = STATS.tstat;
            df(istats)    = STATS.df;
            sd(istats)    = STATS.sd;
            statSum(istats,:) = [idata, icmp, istim, P(istats), tstat(istats), df(istats), sd(istats)];
            
            j = j + 1;
            mLat1(j,:) = [ idata, icmp, istim, 1, mean(alLat_1), std(alLat_1)/sqrt(length(alLat_1)) ];
            j = j + 1;
            mLat1(j,:) = [ idata, icmp, istim, 2, mean(alLat_2), std(alLat_2)/sqrt(length(alLat_2)) ];
            
        end
        %%
        
        lat_1 = mLat1( mLat1(:,1)==idata & mLat1(:,2)==icmp & mLat1(:,4)==1, : );
        lat_2 = mLat1( mLat1(:,1)==idata & mLat1(:,2)==icmp & mLat1(:,4)==2, : );
        
        subplot(2,2,iloc),
        
        hh = bar([lat_1(1,5), lat_2(1,5); lat_1(2,5), lat_2(2,5) ]);
        hh(1).FaceColor = barColor1;
        hh(2).FaceColor = barColor2;
        
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

fid = fopen( [ DIR.figRevision 'barStats_betwnGrp_nGrp'  num2str(nSubGrp) '.csv' ], 'wt' );
fprintf( fid, '%s, %s, %s, %s, %s, %s,%s\n', 'data', 'icmp', 'istim', 'p', 'tstat', 'df', 'sd');
fprintf( fid, '%d, %d,%d, %f, %f, %f, %f\n', statSum' );
fclose(fid);

% set(figure(1),'PaperUnits','inches','PaperPosition',[0 0 6 8])
% 
% switch discriminability
%     case 0
%         print('-dpng', [DIR.figRevision '/IM_halfMax_nGrp' num2str(nSubGrp) '_bar.png'] )
%         save( [DIR.figRevision '/IM_halfMax_nGrp' num2str(nSubGrp) '.mat'], 'statSum', 'mLat1', 'lat' )
%         
%     case 1
%         print('-dpng', [DIR.figRevision '/IM_disc_halfMax_nGrp' num2str(nSubGrp) '_bar.png'] )
%         save( [DIR.figRevision '/IM_disc_halfMax_nGrp' num2str(nSubGrp) '.mat'], 'statSum', 'mLat1', 'lat' )
%         
% end

%% t-test: Compare Button press vs OKN %%% no plotting
% lat: 1)subID, 2)stim, 3)state, 4)nTri, 5)data type, 6)latency

switch nSubGrp
    case {4, 5}
        nCmp = 2;
    case 2
        nCmp = 1;
end

istats = 0;
for icmp = 1:nCmp % we are going to compare 2 sets of subject groups
    for istim = 1:2
        
        switch discriminability
            case 0
                switch istim
                    case 1
                        stim = [2,3];
                    case 2
                        stim = 1;
                end
            case 1
                stim = istim;
        end
        
        switch icmp
            case 1
                switch nSubGrp
                    case 4
                        grp1 = [2, 3];
                        grp2 = [4, 5];
                        cmpName = {'woEle', 'wEle'};
                    case 5
                        grp1 = 2;
                        grp2 = 3;
                        cmpName = {'MON', 'MOF'};
                    case 2
                        grp1 = 1;
                        grp2 = [3, 5];
                        cmpName = {'C', 'PON'};
                        
                end
                
            case 2
                switch nSubGrp
                    case 4
                        grp1 = 1;
                        grp2 = [2, 3, 4, 5];
                        cmpName = {'Control', 'PD'};
                    case 5
                        grp1 = 4;
                        grp2 = 5;
                        cmpName = {'DON', 'DOF'};
                end
        end
        
        for iGrp = 1:2
            switch iGrp
                case 1
                    grp = grp1;
                case 2
                    grp = grp2;
            end
            
            istats = istats + 1;
            alLat_3 = lat( ismember( lat(:,2), stim) & ismember( lat(:,3), grp ) & lat(:,5) == 1, 6 );
            alLat_4 = lat( ismember( lat(:,2), stim) & ismember( lat(:,3), grp ) & lat(:,5) == 2, 6 );
            
            [H_dType,P_dType,CI_dType,STATS_dType] = ttest( alLat_3, alLat_4 ) ;

            statSum_dType(istats,:) = [istim, iGrp, P_dType, STATS_dType.tstat, STATS_dType.df, STATS_dType.sd];
        end
        
    end
end

fid = fopen( [ DIR.figRevision 'barStats_OKNvsBP_nGrp'  num2str(nSubGrp) '.csv' ], 'wt' );
fprintf( fid, '%s,%s,%s,%s,%s,%s\n', 'stim', 'iGrp', 'p', 'tstat', 'df', 'sd');
fprintf( fid, '%d,%d,%f,%f,%f,%f\n', statSum_dType' );
fclose(fid);

% %% CHECK STATS
% 
% statSum
% % statSum
% % [idata, icmp, istim, P(istats), tstat(istats), df(istats), sd(istats)]
% mLat1 
% % mLat -----------------------------------------------------
% % 1)idata, 2)istim, 3)isubGroup, 4)mean:50% of max, 5)SEM  |
% % ----------------------------------------------------------
% 
% %%
% set(figure(1),'PaperUnits','inches','PaperPosition',[0 0 6 8])
% 
% switch discriminability
%     case 0
%         print('-dpng', [DIR.figGrpCmp '/IM_halfMax_nGrp' num2str(nSubGrp) '.png'] )
%         save( [DIR.figGrpCmp '/IM_halfMax_nGrp' num2str(nSubGrp) '.mat'], 'statSum', 'mLat1', 'lat' )
%         
%     case 1
%         print('-dpng', [DIR.figGrpCmp '/' 'IM_disc_halfMax_nGrp' num2str(nSubGrp) '.png'] )
%         save( [DIR.figGrpCmp '/IM_disc_halfMax_nGrp' num2str(nSubGrp) '.mat'], 'statSum', 'mLat1', 'lat' )
%         
% end
% 
% 
% %% Plot bars for mean of max discriminability
% if discriminability
%     figure(5),clf
%     iloc = 0;
%     for idata = 1:2
%         switch idata
%             case 1
%                 dataType  = 'BP';
%             case 2
%                 dataType = 'OKN';
%         end
%         
%         for icmp = 1:2
%             iloc = iloc + 1;
%             
%             switch icmp
%                 case 1
%                     switch nSubGrp
%                         case 4
%                             grp = [1, 2];
%                             cmpName = {'woEle', 'wEle'};
%                         case 5
%                             grp = [2, 3]; %MAYBE [1, 2]
%                             cmpName = {'MON', 'MOF'};
%                     end
%                     
%                 case 2
%                     switch nSubGrp
%                         case 4
%                             grp = [3, 4];
%                             cmpName = {'Control', 'PD'};
%                         case 5
%                             grp = [4, 5];
%                             cmpName = {'DON', 'DOF'};
%                     end
%             end
%             
%             d1 = maxAmp2( maxAmp2(:,1)==idata & maxAmp2(:,2)==1 & maxAmp2(:,3)==grp(1),4:5);
%             d2 = maxAmp2( maxAmp2(:,1)==idata & maxAmp2(:,2)==1 & maxAmp2(:,3)==grp(2),4:5);
%             d3 = maxAmp2( maxAmp2(:,1)==idata & maxAmp2(:,2)==2 & maxAmp2(:,3)==grp(1),4:5);
%             d4 = maxAmp2( maxAmp2(:,1)==idata & maxAmp2(:,2)==2 & maxAmp2(:,3)==grp(2),4:5);
%             
%             subplot(2,2,iloc),
%             bar([d1(1), d2(1); d3(1), d4(1)])
%             k = .15;
%             hold on
%             h = errorbar([1-k, 1+k, 2-k, 2+k], [d1(1), d2(1), d3(1), d4(1)], ...
%                 [d1(2), d2(2), d3(2), d4(2)], [d1(2), d2(2), d3(2), d4(2)], '.k');
%             set(h, 'LineWidth', 3)
%             set(gca, 'XTickLabel',  {'Non-Riv','Riv'})
%             
%             legend(cmpName,'Location', 'SouthOutside')
%             ylabel('Latency(s)')
%             
%             ylim([.8 1])
%             %         xlim([0.3 2.7])
%             title(dataType, 'interpret', 'none')
%             
%         end
%     end
%     
%     set(figure(5),'PaperUnits','inches','PaperPosition',[0 0 6 8])
%     
%     switch discriminability
%         case 0
%             print('-dpng', [DIR.figGrpCmp '/IM_max_nGrp' num2str(nSubGrp) '.png'] )
%             save( [DIR.figGrpCmp '/IM_max_nGrp' num2str(nSubGrp) '.mat'], 'statSum', 'mLat1' )
%             
%         case 1
%             print('-dpng', [DIR.figGrpCmp '/' 'IM_disc_max_nGrp' num2str(nSubGrp) '.png'] )
%             save( [DIR.figGrpCmp '/IM_disc_max_nGrp' num2str(nSubGrp) '.mat'], 'maxAmp2')
%             
%     end
%     
% end


% %% Differential latency
% figure(1),
% l_BP = l(l(:,1)==1,:);
% latency_BP = [l_BP(1:nEff,5)';l_BP(nEff+1:2*nEff,5)']; % row: comparison pair, colomn: effect
%
% l_OKN = l(l(:,1)==2,:);
% latency_OKN = [l_OKN(1:nEff,5)';l_OKN(nEff+1:2*nEff,5)']; % row: comparison pair, colomn: effect
%
% latency_plot = cat(1,latency_BP,latency_OKN);
%
% effName = {'SubGrpup', 'Stimlus', 'Interaction'};
% switch nSubGrp
%     case 4
%         cmpName = {'BP-wEle/woEle', 'BP-C/PD', 'OKN-wEle/woEle', 'OKN-C/PD'};
%     case 5
%         cmpName = {'BP-MON/MOF', 'BP-DON/DOF', 'OKN-MON/MOF', 'OKN-DON/DOF'};
% end
%
% bar([1:4], latency_plot )
% % legend({'BP-wEle/woEle', 'BP-C/PD', 'OKN-wEle/woEle', 'OKN-C/PD'},'Location', 'EastOutside')
% legend(effName,'Location', 'EastOutside')
%
% hold on
% set(gca, 'XTickLabel', cmpName)
%
% ylabel('time from stim onset (sec)')
%
% set(figure(1),'PaperUnits','inches','PaperPosition',[0 0 6 4])
%
% switch discriminability
%     case 0
%         print('-dpng', [DIR.figGrpCmp '/' 'IM_diffLat_nGrp' num2str(nSubGrp) '_tp' num2str(ntp) '.png'] )
%     case 1
%         print('-dpng', [DIR.figGrpCmp '/' 'IM_disc_diffLat_nGrp' num2str(nSubGrp) '_tp' num2str(ntp) '.png'] )
% end
%
% % %% Poplulation latency
% % figure(3), clf
% %
% % % load([ DIR.allSub '/IM_latency_sepOnOff_tp' num2str(ntp) '_' type '.mat'])
% % latency_plot = [latency(1:nSubGrp,4)';latency(nSubGrp+1:2*nSubGrp,4)';...
% %     latency(2*nSubGrp+1:3*nSubGrp,4)';latency(3*nSubGrp+1:4*nSubGrp,4)'];
% %
% % switch nSubGrp
% %     case 5
% %         grpName = {'c', 'mon', 'mof', 'don', 'dof'};
% %     case 4
% %         grpName = {'wEle', 'woEle', 'c', 'PD'};
% % end
% %
% % bar([1:nSubGrp], latency_plot' )
% % legend({'BP-nonRiv', 'BP-Riv', 'OKN-nonRiv', 'OKN-Riv'},'Location', 'EastOutside')
% % hold on
% % set(gca, 'XTickLabel', grpName)
% %
% % ylabel('time from stim onset (sec)')
% % ylim([0 .45])
% %
% % set(figure(3),'PaperUnits','inches','PaperPosition',[0 0 6 4])
% %
% % switch discriminability
% %     case 0
% %         print('-dpng', [DIR.figGrpCmp '/IM_meanAmp_latency_nGrp' num2str(nSubGrp) '_tp' num2str(ntp) '.png'] )
% %     case 1
% %         print('-dpng', [DIR.figGrpCmp '/IM_meanDiscriminability_latency_nGrp' num2str(nSubGrp) '_tp' num2str(ntp) '.png'] )
% % end
% %
% %
% %
% %
% % % figure(1),
% % % set(figure(1),'PaperUnits','inches','PaperPosition',[0 0 25 15])
% % % print('-dpng', [DIR.figIMSortBPSlwPhs '/' 'IM_mSubImage_eachDataEachStimEachSubgrp_sepOnOff_tp' num2str(ntp) '_allpBP.png'] )
% %
% figure(2),
% set(figure(2),'PaperUnits','inches','PaperPosition',[0 0 25 15])
%
% if diffLat
%     switch discriminability
%         case 0
%             print('-dpng', [DIR.figGrpCmp '/' 'IM_mTrimSubPlot_diffLat_nGrp' num2str(nSubGrp) '_tp' num2str(ntp) '.png'] )
%         case 1
%             print('-dpng', [DIR.figGrpCmp '/' 'IM_mTrimSubDiscriminability_diffLat_nGrp' num2str(nSubGrp) '_tp' num2str(ntp) '.png'] )
%
%     end
% else
%
%     switch discriminability
%         case 0
%             print('-dpng', [DIR.figGrpCmp '/' 'IM_mTrimSubPlot_nGrp' num2str(nSubGrp) '_tp' num2str(ntp) '.png'] )
%         case 1
%             print('-dpng', [DIR.figGrpCmp '/' 'IM_mTrimSubDiscriminability_nGrp' num2str(nSubGrp) '_tp' num2str(ntp) '.png'] )
%
%     end
% end
%
%             if and(nSubGrp == 4, ismember(isubgrp, [2,4]) ) ||...
%                     and(nSubGrp == 5, ismember(isubgrp, [3,5]) )
%
%                 if diffLat
%                     plot([latency(1,5), latency(1,5)  ], yLim, 'Color', [0 .6 .8], 'LineWidth', 1.5) %Effect of subGroup, blue
%                     plot([latency(2,5), latency(2,5)  ], yLim, 'Color', [.7 0 0], 'LineWidth', 1.5) % Effect of stimulus, red
%                     plot([latency(3,5), latency(3,5)  ], yLim, 'Color', [0.7 0.3 0], 'LineWidth', 1.5) % Interaction between subGroupd and stim, orange
%                 else
%                     tmp1(padj_subGroup<0.05)= 1;
%                     tmp1(padj_subGroup>=0.05)= nan;
%                     tmp1(isnan(padj_subGroup))= nan;
%                     H_subGrp = plot(time, (yLim(1) + abs(yLim(1)*0) )*tmp1, '.', 'MarkerSize', 15);
%                     set(H_subGrp,'MarkerFaceColor', [0 .6 .8], 'Color', [0 .6 .8])
%
%                     tmp3(padj_stimSubGroup<0.05)= 1;
%                     tmp3(padj_stimSubGroup>=0.05)= nan;
%                     tmp3(isnan(padj_stimSubGroup))= nan;
%                     H_stimSubGroup = plot(time, (yLim(1) + abs(yLim(1)*.2) )*tmp3, '.', 'MarkerSize', 15);
%                     set(H_stimSubGroup,'MarkerFaceColor', [0.7 .3 0], 'Color', [.7 .3 0])
%
%                     tmp2(padj_stim<0.05)= 1;
%                     tmp2(padj_stim>=0.05)= nan;
%                     tmp2(isnan(padj_stim))= nan;
%                     H_stim = plot(time, (yLim(1) + abs(yLim(1)*.4) )*tmp2, '.', 'MarkerSize', 15);
%                     set(H_stim,'MarkerFaceColor', [.7 0 0], 'Color', [.7 0 0])
%                 end
%
%             end
%
%             switch discriminability
%                 case 0
%                     plot([-1 2], [0 0], '--k')
%                 case 1
%                     plot([-1 2], [.5 .5], '--k')
%             end
%             plot([0 0], yLim, '--k')
%
%             title([titlePlot '_' cmpName], 'interpret', 'none')
%
%             switch idata
%                 case 2
%                     ylabel('OKNSlwPhs(deg/sec)')
%                 case 1
%                     ylabel('Consistency')
%             end
%             xlabel('Time from stim onset (sec)')

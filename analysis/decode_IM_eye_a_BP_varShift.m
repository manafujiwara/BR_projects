function decode_IM_eye_a_BP_varShift(DIR, subID)
%% Prepare data
load( [DIR.allSub '/IM_eyeCln.mat']);
load( [DIR.allSub '/IM_BPCln.mat'])
load( [ DIR.allSub '/subInfo.mat'  ] ,'rejSub', 'accSub', 'subInfo', 'badTri')

decode_stim = 0;
% -------------------------------------------------------------------------
% subInfo: 1) subNum, 2)nallTri, 3)nBadTri(eye30)), 4)nBadTri(eye150)) ,...
%         5)rejBPeye30(rate), 6)rejBPeye150(rate)];
% -------------------------------------------------------------------------
% BPCln: 1)subID, 2)irun, 3)iblk, 4)iepc, 5)stimLabel,
% 6)subjestState, 7)BPP (-1~1), 8)label, 9)fBP, 10)fcBP, 11:910)BP(t)
% -------------------------------------------------------------------------
% eyeinfo: 1)subID, 2)irun, 3)iblk, 4)iepc, 5)stimLabel, 6)subjestState
% 7)BPP (-1~1), 8)label, 9)fBP, 10)fcBP, 11:910)eyeinfo(t)
% -------------------------------------------------------------------------

time = linspace(-1, 2, 900);
ntp = 20;
nSubGrp = 5; % 4 or 5. 4)with electrode, without electrode, control, pd. 5) c, mon, mof, don, dof
iLocPlot = 0;
iLocImage = 0;
jsub = 0;
%% REMOVE SUBJECT WHO HAVE DONE EITHER ON- OR OFF- TREATMENT
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


%% HERE FIX WRONG SUBJECT INFO
BPCln( ismember( BPCln(:,1),[162, 166] ), 6 ) = 3;
slwPhsBoxCln( ismember( slwPhsBoxCln(:,1),[162, 166] ), 6 ) = 3;

%%
loadFileNameBP = [DIR.SVM '/IM_discriminability_nGrp' num2str(nSubGrp) '_tp' num2str(ntp) '_BP_balTri.mat'];
loadFileNameOKN = [DIR.SVM '/IM_discriminability_nGrp' num2str(nSubGrp) '_tp' num2str(ntp) '_OKN30_balTri.mat'];

if exist(loadFileNameBP, 'file') && exist(loadFileNameOKN, 'file')
    load(loadFileNameBP, 'BP')
    load(loadFileNameOKN, 'OKN30')
    
else
    
    for idata = 1:2 % BP or OKN
        
        for istim = 1:2 % non-rivalry or br
            if decode_stim && istim == 2, break, end
            
            switch idata
                case 1 % data = BP
                    titleData = 'BP';
                    
                    switch istim
                        case 1 % stim = non-br 1)subID, 2)stimLabel, 3)subjectState, 4) label
                            D = BPCln( ismember( BPCln(:,5), [2,3]) , [1,5,6,8,11:end] );
                            titleStim = '_non';
                            
                        case 2 % stim = br
                            titleStim = '_br';
                            D = BPCln( ismember( BPCln(:,5), 1 ) , [1,5,6,8,11:end] );
                    end
                    
                case 2 % data = eye
                    titleData = 'OKN30';
                    
                    switch istim
                        case 1 % stim = non-br
                            titleStim = '_non';
                            D = slwPhsBoxCln( ismember( slwPhsBoxCln(:,5), [2,3]) , [1,5,6,8,11:end] );
                            
                        case 2 % stim = br
                            titleStim = '_br';
                            D = slwPhsBoxCln( ismember( slwPhsBoxCln(:,5), 1) , [1,5,6,8,11:end] );
                            yLim = [-.5 1];
                    end
            end
            
            
            for isubgrp = 1:nSubGrp % control, medOn, medOff, dbsOn, dbsOff
                
                switch nSubGrp
                    case 5
                        subgrp = isubgrp;
                    case 4
                        switch isubgrp
                            case 1
                                subgrp = [2,3];
                            case 2
                                subgrp = [4,5];
                            case 3
                                subgrp = 1;
                            case 4
                                subgrp = [2,3,4,5];
                        end
                    case 2
                        switch isubgrp
                            case 1
                                subgrp = 1;
                            case 2
                                subgrp = [3,5];
                        end
                        
                end
                
                DD = D(ismember(D(:,1), accSub) & ~ismember(D(:,1), oneStateSub) & ismember(D(:,3), subgrp), :);
                
                grpSub = unique(DD(:,1));
                nSub = length(grpSub);
                
                for isub = 1:nSub
                    
                    subNum = grpSub(isub);
                    if ~ismember (subNum, accSub), continue, end
                    
                    jsub = jsub + 1;
                    DD(isnan(DD(:,4)),:) = [];
                    
                    switch istim
                        case 1
                            % DD 1)subID, 2)stimLabel, 3)subjectState, 4) label,
                            % 5:904) smplitude
                            rdd = DD( DD(:,1)==grpSub(isub) & DD(:,2) == 2, 5:904); %r-stim trials
                            ldd = DD( DD(:,1)==grpSub(isub) & DD(:,2) == 3, 5:904); %l-stim trials
                        case 2
                            rdd = DD( DD(:,1)==grpSub(isub) & DD(:,4) == 1, 5:904); %r-labeled trials
                            ldd = DD( DD(:,1)==grpSub(isub) & DD(:,4) == -1, 5:904); %l-labeled trials
                    end
                    
                    nrTri = size(rdd,1);
                    nlTri = size(ldd,1);
                    
                    %% REVISION - Weighted SVM ( No crossvalidation ) --
                    
                    %%% Put weight and run SVM with unbalance nTrial %%%
                    %                     nTri_trn = floor(.7 * [ nlTri, nrTri ]); % 1)left, 2)right
                    %                     nTri_tst = round(.3 * [ nlTri, nrTri ]);
                    
                    %                     label_trn( 1 : nTri_trn(2) ) = {'r'};
                    %                     label_trn( nTri_trn(2)+1 : sum(nTri_trn) ) = {'l'};
                    %                     label_tst( 1 : nTri_tst(2) ) = {'r'};
                    %                     label_tst( nTri_tst(2)+1 : sum(nTri_tst) ) = {'l'};
                    %
                    %                     if nlTri >= nrTri % if n left trial was more than right
                    %                         weight_trn( 1 : nTri_trn(2) ) = max(nTri_trn)/min(nTri_trn);
                    %                         weight_trn( nTri_trn(2)+1 : sum(nTri_trn) ) = 1;
                    %                     elseif nlTri < nrTri
                    %                         weight_trn( 1 : nTri_trn(2) ) = 1; % if n right tri was more than left
                    %                         weight_trn( nTri_trn(2)+1 : sum(nTri_trn) ) = max(nTri_trn)/min(nTri_trn);
                    %                     end
                    
                    %%% Balance nTrial and run SVM with cross validation %%%
                    
                    nTri_dcd = min( [nrTri, nlTri] );
                    nTri_trn = floor(.7 * nTri_dcd);
                    nTri_tst = round(.3 * nTri_dcd);
                    
                    label_trn( 1 : nTri_trn ) = {'r'};
                    label_trn( nTri_trn+1 : 2*nTri_trn ) = {'l'};
                    label_tst( 1 : nTri_tst ) = {'r'};
                    label_tst( nTri_tst+1 : 2*nTri_tst ) = {'l'};
                    
                    for icrsVal = 1:10
                        
                        watch = clock;
                        msg = [ 'data=' titleData titleStim ', subGrp=' num2str(isubgrp) ', subNum=' num2str(subNum) ...
                            ', crossVal=' num2str(icrsVal) ', ' num2str(watch(4)) ':' num2str(watch(5)) '''' num2str(watch(6))];
                        disp(msg)
                        
                        idxrTri = randperm(nrTri);
                        idxlTri = randperm(nlTri);
                        
                        idxrTri_trn = idxrTri( 1 : nTri_trn );
                        idxrTri_tst = idxrTri( nTri_trn+1 : nTri_trn + nTri_tst );
                        
                        idxlTri_trn = idxlTri( 1 : nTri_trn );
                        idxlTri_tst = idxlTri( nTri_trn+1 : nTri_trn + nTri_tst );
                        
                        d_trn = cat( 1, rdd( idxrTri_trn, : ), ldd( idxlTri_trn, : ) );
                        d_tst = cat( 1, rdd( idxrTri_tst, : ), ldd( idxlTri_tst, : ) );
                        
                        Y = label_trn;
                        for isamp = 1:900
                            X = d_trn(:,isamp);
                            SVMModel = fitcsvm(X,Y,'KernelFunction','linear');
                            [testLabel,score] = predict( SVMModel, d_tst(:,isamp) );
                            decAcc(icrsVal,isamp) = sum( strcmp( label_tst', testLabel ) ) / length( label_tst );
                        end
                        
                    end
                    
                    clear label_trn label_tst weight_trn
%                     % 1)subID, 2)stim, 3)isubGroup, 4)ratio nlTri/nrTri, 5)dataType, 6:905)mTriAmplitude
%                     mTri_DD(jsub,:) = [ subNum, istim, isubgrp, nTri_trn(1)/nTri_trn(2), idata, mean(decAcc,1) ];
                    
                     % 1)subID, 2)stim, 3)isubGroup, 4)nTri4decoding, 5)dataType, 6:905)mTriAmplitude
                    mTri_DD(jsub,:) = [ subNum, istim, isubgrp, nTri_dcd, idata, mean(decAcc,1) ];
                    
                end
                
            end
        end
        
        switch idata
            case 1
                dataType = 'BP';
                BP = mTri_DD;
                save([ DIR.SVM '/IM_discriminability_nGrp' num2str(nSubGrp) '_tp' num2str(ntp) '_' dataType '_balTri.mat' ], 'BP' )
            case 2
                dataType = 'OKN30';
                OKN30 = mTri_DD;
                save([ DIR.SVM '/IM_discriminability_nGrp' num2str(nSubGrp) '_tp' num2str(ntp) '_' dataType '_balTri.mat' ], 'OKN30' )

        end
        jsub = 0;
        mTri_DD = [];
        
    end
        
end

% %% Compute latencies and plot descriminability
% figure(1),clf
% iLocImage = 0;
% iLocPlot  = 0;
% yLim = [ .4 1.1];
% 
% for idata = 1:2
%     
%     switch idata
%         case 1
%             mTri_DD = BP;
%             titleData = 'BP';
%         case 2
%             mTri_DD = OKN30;
%             titleData = 'OKN30';
%     end
%     
%     for istim = 1:2
%         
%         switch istim
%             case 1
%                 titleStim = '_non';
%             case 2
%                 titleStim = '_br';
%         end
%         
%         for isubgrp = 1:nSubGrp
%             iLocImage = iLocImage + 1;
%             switch nSubGrp
%                 case 5
%                     switch isubgrp
%                         case 1
%                             iLocPlot = iLocPlot+ 1;
%                             lineStyle = '-';
%                             lineColor = [.8 0 0];
%                             cmpName = 'control';
%                             
%                         case 2
%                             iLocPlot = iLocPlot+ 1;
%                             lineStyle = '-';
%                             lineColor = [0 .8 0];
%                             cmpName = 'mon-mof';
%                             
%                         case 3
%                             lineStyle = '-';
%                             lineColor = [0 .2 0];
%                             cmpName = 'mon-mof';
%                             
%                         case 4
%                             iLocPlot = iLocPlot+ 1;
%                             lineStyle = '-';
%                             lineColor = [0 0 1];
%                             cmpName = 'don-dof';
%                             
%                         case 5
%                             lineStyle = '-';
%                             lineColor = [0 0 .2];
%                             cmpName = 'don-dof';
%                     end
%                     
%                 case 4
%                     
%                     switch isubgrp
%                         case 1 % without electrode
%                             iLocPlot = iLocPlot+ 1;
%                             lineStyle = '-';
%                             lineColor = [0 .5 0];
%                             cmpName = 'w-woElectrode';
%                             
%                         case 2 % with electrode
%                             lineStyle = '-';
%                             lineColor = [0 0 .6];
%                             cmpName = 'w-woElectrode';
%                             
%                         case 3 % control
%                             iLocPlot = iLocPlot+ 1;
%                             lineStyle = '-';
%                             lineColor = [0 0 .9];
%                             cmpName = 'c-pd';
%                             
%                         case 4 % PD patients
%                             lineStyle = '-';
%                             lineColor = [.8 0 0];
%                             cmpName = 'c-pd';
%                     end
%                     
%                 case 2
%                     
%                     switch isubgrp
%                         case 1 % control
%                             iLocPlot = iLocPlot+ 1;
%                             lineStyle = '-';
%                             lineColor = [0 0 .9];
%                             cmpName = 'c-pon';
%                             
%                         case 2 % PD patients
%                             lineStyle = '-';
%                             lineColor = [.8 0 0];
%                             cmpName = 'c-pon';
%                     end
%                     
%             end
%             
%             %% Compute significance %%%
%             % 1)subID, 2)stim, 3)state, 4)nTri, 5)dataType,
%             % 6:905)mTriAmplitude(discriminability)
%             amp = mTri_DD( mTri_DD(:,2)== istim & mTri_DD(:,3)==isubgrp & mTri_DD(:,5)== idata, : );
%             
%             nSub = size(amp,1);
%             
%             mAmp_plot = nanmean(amp(:,6:905),1);
%             semAmp_plot = nanstd(amp(:,6:905),1)/sqrt(nSub);
%             %             maxAmp(iLocImage,:) = [ idata, istim, isubgrp, max(mAmp_plot) min(time(find(mAmp_plot == max(mAmp_plot)))) ];
%             %
%             %             for i = 1:size(amp,1)
%             %                 tMaxAmp(i) = time( 300 + min( find( amp( i, 306:905 ) == max( amp( i, 306:905 ) ) ) ) );
%             %                 tmp = time( 300 + min( find( amp( i, 306:905 ) >= 0.75 ) ) );
%             %                 if isempty(tmp)
%             %                     t75Amp(i) = NaN;
%             %                 else
%             %                     t75Amp(i) =  tmp;
%             %                 end
%             %             end
%             %
%             %             maxAmpSub(iLocImage).d = [repmat(idata,size(amp,1),1), repmat(istim,size(amp,1),1),...
%             %                 repmat(isubgrp,size(amp,1),1), max(amp(:,306:905),[],2), tMaxAmp', amp(:,1)];
%             %
%             %             maxAmpSub(iLocImage).pct75 = [repmat(idata,size(amp,1),1), repmat(istim,size(amp,1),1),...
%             %                 repmat(isubgrp,size(amp,1),1), repmat(0.75,size(amp,1),1), t75Amp', amp(:,1)];
%             
%             %             clear tMaxAmp t75Amp
%             
%             for i = 1:1000
%                 iv = randi([1,nSub], 1, nSub);
%                 mAmp(i,:) = nanmean(amp(iv,6:905),1);
%             end
%             
%             for itp = 1:900
%                 dist = sort(mAmp(:,itp), 1, 'descend');
%                 h1(iLocImage,itp) = dist(950);
%                 %if 95% of the means are higher than 0.5 (chance)
%                 %the discriminability is significantly higher than chance.
%             end
%             
%             increasing = round(h1(iLocImage,:)*100)/100>.5;
%             i_all = 0;
%             while 1
%                 i_all = i_all + 1;
%                 if sum( increasing( 300 + i_all : 300 + i_all + ntp-1 ) ) == ntp
%                     break
%                 elseif i_all+ 300  == length(increasing)-ntp-1
%                     i_all = nan;
%                     break
%                 end
%             end
%             
%             if isnan(i_all)
%                 latency(iLocImage,:) = [ idata, istim, isubgrp, nan ];
%             else
%                 latency(iLocImage,:) = [ idata, istim, isubgrp, time(300+i_all) ];
%             end
%             
%             
%             % Let's plot
%             figure(1), % SHADEDERRORBAR ===================================
%             
%             switch nSubGrp
%                 case 5
%                     subplot(5,3,iLocPlot)
%                 case {2,4}
%                     subplot(5,2,iLocPlot)
%             end
%             titlePlot = cat(2, titleData, titleStim);
%             
%             hh = shadedErrorBar(time, mAmp_plot, semAmp_plot, [], 1);
%             hold on
%             set(hh.mainLine, 'LineWidth', 2, 'LineStyle', lineStyle , 'Color', lineColor)
%             set(hh.edge,'Color', lineColor)
%             set(hh.patch,'FAceColor', lineColor)
%             
%             ylim(yLim)
%             
%             hhh = plot( [ latency(iLocImage, 4), latency(iLocImage, 4)], yLim );
%             set(hhh, 'LineWidth', 2 , 'Color', lineColor, 'LineStyle', lineStyle)
%             
%             plot([-1 2], [.5 .5], '--k')
%             plot([0 0], yLim, '--k')
%             
%             title([titlePlot '_' cmpName], 'interpret', 'none')
%             ylabel('Discriminability')
%             xlabel('Time from stim onset (sec)')
%         end
%     end
% end
% 
% figure(1),
% set(figure(1),'PaperUnits','inches','PaperPosition',[0 0 25 15])
% save([DIR.SVM '/IM_discriminability_nGrp' num2str(nSubGrp) '_tp' num2str(ntp) '.mat'], 'h1', 'latency', 'BP', 'OKN30')
% print('-dpng', [DIR.figRevision '/IM_discriminability_nGrp' num2str(nSubGrp) '_tp' num2str(ntp) '_rev.png'] )
% 
% 
% 
% 
% % Plot bar graph for latencies.
% figure(3), clf
% 
% latency_plot = [latency(1:nSubGrp,4)';latency(nSubGrp+1:2*nSubGrp,4)';...
%     latency(2*nSubGrp+1:3*nSubGrp,4)';latency(3*nSubGrp+1:4*nSubGrp,4)'];
% 
% switch nSubGrp
%     case 5
%         grpName = {'c', 'mon', 'mof', 'don', 'dof'};
%     case 4
%         grpName = {'wEle', 'woEle', 'c', 'PD'};
%     case 2
%         grpName = {'C', 'PON'};
% end
% 
% bar([1:nSubGrp], latency_plot' )
% legend({'BP-nonRiv', 'BP-Riv', 'OKN-nonRiv', 'OKN-Riv'},'Location', 'EastOutside')
% hold on
% set(gca, 'XTickLabel', grpName)
% 
% ylabel('time from stim onset (sec)')
% 
% switch nSubGrp
%     case 5
%         ylim([0 .5])
%     case 4
%         ylim([0 .35])
% end
% 
% set(figure(3),'PaperUnits','inches','PaperPosition',[0 0 6 4])
% print('-dpng', [DIR.figRevision '/IM_discriminability_latency_nGrp' num2str(nSubGrp) '_tp' num2str(ntp) '_rev.png'] )
% 
% 
% %% Make file to run stats in R
% 
% for icmp = 1:2
%     
%     switch nSubGrp
%         case 5
%             switch icmp
%                 case 1 %compare med-on and med-off
%                     grp1 = 2;
%                     grp2 = 3;
%                     cmpName = 'mon-mof';
%                 case 2 % compare dbs-on and dbs-off
%                     grp1 = 4;
%                     grp2 = 5;
%                     cmpName = 'don-dof';
%             end
%         case 4
%             switch icmp
%                 case 1 % compare with and without electrode
%                     grp1 = 1;
%                     grp2 = 2;
%                     cmpName = 'w-woElectrode';
%                 case 2 % compare control and PD
%                     grp1 = 3;
%                     grp2 = 4;
%                     cmpName = 'c-pd';
%             end
%     end
%     % [mTri_DD]
%     % 1)subID, 2)stim, 3)isubgroup(IN COMPARISON, NOT STATE ITSELF), 4)label, 5)nTri, 6)data
%     % type, 7:906)mTriAmplitude
%     D_cmp    = [];
%     D_stimlb = [];
%     D_grplb  = [];
%     D_type   = [];
%     clear amplitude
%     for idata = 1:2
%         switch idata
%             case 1
%                 mTri_DD = BP;
%                 type = 'BP';
%             case 2
%                 mTri_DD = OKN30;
%                 type = 'OKN30';
%         end
%         % 1)mTri_DD 1)subID, 2)stim, 3)isubgrp(NOT STATE), 4)nTri, 5)data
%         % type, 6:905)mTriAmplitude
%         
%         D_stim1_grp1 = mTri_DD( ismember( mTri_DD(:,2),  [2,3]) &...
%             mTri_DD(:,3) == grp1 & mTri_DD(:,5) == idata,...
%             [1, 6:905 ] );
%         D_stim1_grp2 = mTri_DD( ismember( mTri_DD(:,2), [2,3] ) & ...
%             mTri_DD(:,3)==grp2 & mTri_DD(:,5) == idata,...
%             [1, 6:905] );
%         D_stim2_grp1 = mTri_DD( mTri_DD(:,2) == 1 &...
%             mTri_DD(:,3)==grp1 & mTri_DD(:,5) == idata,...
%             [1, 6:905] );
%         D_stim2_grp2 = mTri_DD( mTri_DD(:,2) == 1 & ...
%             mTri_DD(:,3)==grp2  & mTri_DD(:,5) == idata,...
%             [1, 6:905] );
%         
%         D_cmp = [ D_stim1_grp1; D_stim1_grp2; D_stim2_grp1; D_stim2_grp2];
%         
%         D_type = [ repmat({type}, [ size([D_stim1_grp1; D_stim1_grp2; D_stim2_grp1; D_stim2_grp2], 1) ], 1) ];
%         
%         D_stimlb = [ repmat({'non-br'}, [ size(D_stim1_grp1, 1), 1 ] );...
%             repmat({'non-br'}, [ size(D_stim1_grp2, 1 ), 1 ] );...
%             repmat({'br'}, [ size(D_stim2_grp1, 1 ), 1 ] );...
%             repmat({'br'}, [ size(D_stim2_grp2, 1 ), 1 ] ) ] ;
%         
%         D_grplb = [ repmat({'subgrp1'}, [ size( D_stim1_grp1, 1 ), 1 ] );...
%             repmat({'subgrp2'}, [ size( D_stim1_grp2, 1 ), 1 ] );...
%             repmat({'subgrp1'}, [ size( D_stim2_grp1, 1 ), 1 ] );...
%             repmat({'subgrp2'}, [ size( D_stim2_grp2, 1 ), 1 ] ) ] ;
%         
%         nSub_type= size(D_type,1);
%         
%         amplitude2(:,1)     = num2cell(D_cmp( 1: nSub_type,1)) ;
%         amplitude2(:,2)     = D_stimlb( 1: nSub_type ) ;
%         amplitude2(:,3)     = D_grplb(1: nSub_type);
%         amplitude2(:,4:903) = num2cell(D_cmp(1: nSub_type, 2:901 ));
%         
%         
%         fileID = fopen([ '../../R/BR_rev/data/stats_disc_cmp-' cmpName '_tp' num2str(ntp) '_' type '.txt'],'w');
%         
%         formatSpec = [ '%d %s %s' strjoin( repmat( {' %d'}, [1,899] )  ) ' %d\n' ];
%         
%         [nrows,ncols] = size(amplitude2);
%         for row = 1:nrows
%             fprintf(fileID,formatSpec,amplitude2{row,:});
%         end
%         
%         fclose(fileID);
%         clear amplitude2
%     end
% end

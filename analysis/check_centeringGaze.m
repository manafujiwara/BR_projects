function check_centeringGaze(DIR,subID, Cond)
% Please note that statistic analysis later in this code is written only
% for comparison between Control vs ON-treatment PD, although it should not
% be difficult to add the other group comparisons.

%-----------------------------------------------------------------------
% eyeinfo: 1)subID, 2)irun, 3)iblk, 4)iepc, 5)stimLabel, 6)subjestState |
% 7)BPP (-1~1), 8)label, 9)fBP, 10)fcBP, 11:910)eyeinfo(t)              |
% -----------------------------------------------------------------------
whicheye = lower(Cond.whicheye);

time = linspace(-1,2,900);
switch whicheye
    case 'left'
        whichEyeColumn = 3;
    case 'right'
        whichEyeColumn = 5;
end

time2cut = Cond.time2cut;
epochWin = Cond.epochWinIM;
load( [ DIR.allSub '/subInfo.mat' ] ,'rejSub', 'accSub', 'subInfo', 'badTri')
load( [DIR.cfg '/c106_Cfg_a_BR.mat' ] )

ksub = 0;
subID_grp = unique(subInfo(:,1));
for isub = 1:length(subID_grp)
    state = unique(subInfo(subInfo(:,1)==subID_grp(isub),2));
    if length(state)==1 && state == 1
        continue
    elseif length(state)==1 && state ~= 1
        ksub = ksub + 1;
        oneStateSub(ksub,1) = subID_grp(isub);
        %         oneStateSub(ksub,2) = state;
    end
end

loadFileName =[ DIR.allSub '/IM_filtGazeMat.mat' ];

if ~exist(loadFileName,'file')
    
    jepc = 0;
    
    for isub = 1: length(subID)
        
        subName = subID{isub};
        subNum = str2num(subName(2:end));
        
        if ismember( subNum, rejSub )
            disp(['filter epoch & reshape. ' subName ' skipped as the participant is rejected from the analysis' ])
            continue
        end
        
        loadSubFileName = [ DIR.collectedData '/' subID{isub} '.mat'];
        load( loadSubFileName, 'Run' )
        
        for irun = 1 : length(Run)
            
            for iblk = 1: length(Run(irun).block)
                
                %%% Remove saccades %%%
                tMissP = [];
                
                gazeData = Run(irun).block(iblk).gazeData;
                
                if ~isempty(gazeData)
                    
                    if ismember( subNum, [127, 130, 131, 149] )
                        % These accepted 4 subjects were computed in pixels as
                        % the computation was done after modification. Now put
                        % them back in ratio of the screen
                        gazeData(:, [3 5]) = gazeData(:, [3 5]) / Cfg.width; % in pixel
                        gazeData(:, [4 6]) = gazeData(:, [4 6]) / Cfg.height; % in pixel
                    end
                    
                    blinks = Run(irun).block(iblk).blinks;
                    allsac = Run(irun).block(iblk).allsac;
                    
                    if ~isempty(blinks)
                        if ~isnan(blinks)
                            tMissP = blinks ;
                        end
                    end
                    
                    if ~isempty(Run(irun).block(iblk).allsac)
                        tMissP = [tMissP; allsac(:, [1 2]) ];
                    end
                    
                    %%% adjust the length so that they can be comparable with button press responses
                    nPoints = (Run(irun).block(iblk).stimFrames + ...
                        Run(irun).block(iblk).ISI ) * Run(irun).block(iblk).nTrials * 5;
                    
                    if length(gazeData) > nPoints
                        gazeData = gazeData(1:nPoints, :);
                    end
                    
                    rawGazeData = gazeData;
                    
                    %%% replace missing points with nan
                    tMissP = sortrows ( tMissP );
                    tMissP2cut = [ tMissP(:,1) - time2cut, tMissP(:,2) + time2cut ];
                    
                    for iMissP = 1 : size(tMissP2cut,1)
                        gazeData ( gazeData ( : , 17 ) >= tMissP2cut (iMissP,1)...
                            & gazeData ( : , 17 ) <= tMissP2cut (iMissP,2) , 3:10 ) = nan;
                    end
                    
                    if 0
                        plot_raw2CheckBlinkSacc
                    end
                    
                end
                
                filtGazeData = gazeData; %- missing points are replaced by NaN
                Run(irun).block(iblk).filtGazeData = filtGazeData;
                
                data = filtGazeData;
                if ~isempty(data)
                    
                    %% EPOCH EYE DATA BASED ON TIME ===============================
                    if ~strcmp(Run(irun).block(iblk).id,'2-1')
                        continue
                    end
                    
                    epochTimes(:,1) = Run(irun).block(iblk).trial_evTimes(:,1) + repmat(epochWin(1), Run(irun).block(iblk).nTrials, 1);
                    epochTimes(:,2) = epochTimes(:,1) + ( diff(epochWin) ) ;
                    
                    for itri = 1: Run(irun).block(iblk).nTrials
                        Run(irun).block(iblk).epoch(itri).labelDes = Run(irun).block(iblk).mixedLabelDesign(itri);
                        
                        data_epoch = data(data(:,17,1) >= epochTimes(itri, 1) & data(:,17) < epochTimes(itri, 2), : );
                        
                        if size(data_epoch,1) < diff(epochWin)*300
                            continue
                            
                        elseif size(data_epoch,1) >  diff(epochWin)*300
                            data_epoch = data_epoch( 1 : diff(epochWin)*300, : , :);
                        end
                        
                        Run(irun).block(iblk).epoch(itri).filtGazeData = data_epoch;
                        
                    end % trial
                end
                
                %% RESHAPE STRUCTURE TO MATRIX ================================
                
                switch Run(irun).block(iblk).state
                    case 'nc'
                        state= 1; %=> block level
                    case 'mon'
                        state = 2;
                    case 'mof'
                        state = 3;
                    case 'donmon'
                        state = 4;
                    case 'dofmon'
                        state = 5;
                    case 'donmof'
                        state = 4;
                    case 'dofmof'
                        state = 5;
                end
                
                if isempty(Run(irun).block(iblk).gazeData)
                    
                    baseInfo(1:20,1) = str2num(subName(2:end));
                    baseInfo(1:20,2) = irun;
                    baseInfo(1:20,3) = iblk;
                    baseInfo(1:20,4) = 1:20;
                    baseInfo(1:20,5) = nan; %label
                    baseInfo(1:20,6) = state;
                    
                    filtGaze_x(jepc + 1 : jepc + 20, : )  = [ baseInfo, nan(20,900)];
                    filtGaze_y(jepc + 1 : jepc + 20, : )  = [ baseInfo, nan(20,900)];
                    
                    jepc = jepc + 20;
                    clear baseInfo
                    continue
                end
                
                for iepc = 1: length(Run(irun).block(iblk).epoch)
                    jepc = jepc + 1;
                    
                    label= Run(irun).block(iblk).epoch(iepc).labelDes; %=> epoch level. % 1)br, 2)red/right unambiguous, 3)green/left unambiduous
                    
                    baseInfo(1) = str2num(subName(2:end));
                    baseInfo(2) = irun;
                    baseInfo(3) = iblk;
                    baseInfo(4) = iepc;
                    baseInfo(5) = label;
                    baseInfo(6) = state;
                    
                    filtGaze_x(jepc, 1:6 )  = baseInfo;
                    filtGaze_y(jepc, 1:6 )  = baseInfo;
                    
                    % ---------------------------------------------------------------
                    % 1)subID, 2)irun, 3)iblk, 4)iepc, 5)stimLabel, 6)subjestState  |
                    % ---------------------------------------------------------------
                    
                    if isempty(Run(irun).block(iblk).epoch(iepc).filtGazeData)
                        filtGaze_x(jepc,7:906) = nan(1,900);
                        filtGaze_y(jepc,7:906) = nan(1,900);
                        
                    else
                        filtGaze_x(jepc, 7:906) = Run(irun).block(iblk).epoch(iepc).filtGazeData(:,whichEyeColumn)'; %<= Cond.slowPhaseSmoothWin
                        filtGaze_y(jepc, 7:906) = Run(irun).block(iblk).epoch(iepc).filtGazeData(:,whichEyeColumn+1)'; %<= Cond.slowPhaseSmoothWin
                    end
                    
                end
                
                
            end % blk
            
            if 0
                plot_gaze2Check
            end
            
        end % run
        
        saveFileName = [ subID{isub} '_filtGaze_' num2str(time2cut*10^3) '.mat' ];
        save( [DIR.filtGaze '/' saveFileName ], 'Run' )
        
        saveFileName = [ subID{isub} '_filtGaze_' num2str(time2cut*10^3) '.mat' ];
        load( [DIR.filtGaze '/' saveFileName ], 'Run' )
        
        watch = fix(clock);
        msg = sprintf('filter epoch & reshape done with sub: %s | %s:%s''''%s',...
            subID{isub},  num2str(watch(4)), num2str(watch(5)), num2str(watch(6)) ); disp(msg);
        
        clear Run
        
    end % sub
    
    saveFileName = [ DIR.allSub '/IM_filtGazeMat.mat' ];
    save( saveFileName, 'filtGaze_x', 'filtGaze_y' )
    
end


%% COMPUTE DISTANCE FROM THE CENTRE OF THE STIMULUS =======================
load( loadFileName )
load( [DIR.allSub '/IM_BPCln.mat'], 'BPCln' ) % load to get accepted trial information

pix2deg = mean( [ Cfg.visualAngleDegPerPixelX, Cfg.visualAngleDegPerPixelY ] ) ;

stimCent = [ mean(BR.framePos([1,3],1)), mean(BR.framePos([2,4],1) ) ]; %left eye (in pixel)

filtCln_x = filtGaze_x( ismember( filtGaze_x(:,1:6), BPCln(:,1:6 ), 'rows' ), : );
filtCln_y = filtGaze_y( ismember( filtGaze_y(:,1:6), BPCln(:,1:6 ), 'rows' ), : );

if size(filtCln_x,1) == size(BPCln,1)
    disp ('Accepted trials are selected properly')
else
    disp ('Some trials are missing')
end


%%% Find trials that have already been computed in pixels
tmp_filtCln_x = filtCln_x;
tmp_filtCln_y = filtCln_y;
tmp_filtCln_x( isnan( tmp_filtCln_x ) ) = 0;
tmp_filtCln_y( isnan( tmp_filtCln_y ) ) = 0;
big_x = tmp_filtCln_x( sum( tmp_filtCln_x(:,7:906) > 100, 2 ) > 0, : );
big_y = tmp_filtCln_y( sum( tmp_filtCln_y(:,7:906) > 100, 2 ) > 0, : );
big_x(big_x==0) = nan;
big_y(big_y==0) = nan;

filtCln_x(:,7:906) = filtCln_x(:,7:906) * Cfg.windowRect(3);
filtCln_y(:,7:906) = filtCln_y(:,7:906) * Cfg.windowRect(4); %(in pixel)

%%% Check the scale of data %%%
if ~isempty(big_x) || ~isempty(big_y)
    filtCln_x( sum( tmp_filtCln_x(:,7:906) > 100, 2 ) > 0, : ) = big_x;
    filtCln_y( sum( tmp_filtCln_y(:,7:906) > 100, 2 ) > 0, : ) = big_y;
end

%%% Now all in pixels %%%
distStimCent_x = filtCln_x(:,307:906) - stimCent(1);
distStimCent_y = filtCln_y(:,307:906) - stimCent(2);

distStimCent = sqrt( distStimCent_x.^2 + distStimCent_y.^2 );

mDistStimCent = [ filtCln_x(:, 1:6), nanmean(distStimCent,2) ]; % average distance in each trial

for iset =  1: length(subInfo)
    
    subNum = subInfo(iset,1);
    state = subInfo(iset,2);
    
    if ismember( subNum, rejSub )
        mDistStimCent_set( iset, : ) = [ subNum, state, nan ];
    else
        mDistStimCent_set( iset, : ) = [ subNum, state, mean( mDistStimCent( mDistStimCent(:,1) == subNum & mDistStimCent(:,6) == state, 7 ) )];
    end
    
end

mDistStimCent_set( isnan( mDistStimCent_set( :, 3 ) ), : ) = [];
mDistStimCent_set( ismember( mDistStimCent_set(:,1), oneStateSub ), : ) = [];

for istate = 1:5
    %
    %         switch istate
    %             case 1
    %                 mDistStimCent_C = mDistStimCent_set( mDistStimCent_set(:,2)==1, 3 );
    %                 mDistStimCent_state(istate) = nanmean( mDistStimCent_C );
    %             case 2
    %                 mDistStimCent_PD = mDistStimCent_set( ismember( mDistStimCent_set(:,2), 2:5 ), 3 ) ;
    %                 mDistStimCent_state(istate) = nanmean( mDistStimCent_PD );
    %         end
    
    tmp = mDistStimCent_set( mDistStimCent_set(:,2)==istate, 3 ) *pix2deg;
    mDistStimCent_state(istate) = nanmean( tmp );
    
    switch istate
        case 1
            dist_c = tmp; % in visual angle
        case 2
            dist_mon = tmp;
        case 3
            dist_mof = tmp;
        case 4
            dist_don = tmp;
        case 5
            dist_dof = tmp;
    end
    
end


%% Plot figures % ONLY FOR nGroup=2 comparison (Control vs treatment-ON PD)

figure(1),clf
subplot(1,2,1)
hist(dist_c, .8:.2:2.8);
h = findobj(gca,'Type','patch');
h.FaceColor = [0 0 .9];
ylabel('nSubject')
xlabel('Mean distance from the center (deg)', 'fontsize', 8)
title('Control')

subplot(1,2,2)
hist([ dist_mon; dist_don ], .8:.2:2.8)
h = findobj(gca,'Type','patch');
h.FaceColor = [.8 0 0];
ylabel('nSubject')
xlabel('Mean distance from the center (deg)', 'fontsize', 8)
title('ON-treatment PD')

set( figure(1),'PaperUnits','inches','PaperPosition',[0 0 5 3])
print('-dpng', [DIR.figRevision '/IM_distribution_centGaze_nGrp2.png'] )


%% Statistic test %%%%%%%%%%%%%%%%%%%%%%%%%%%%
[H, P, CI, stats] = ttest2( dist_c, [ dist_mon; dist_don ] ); % PD On-Treatment vs controls
statData(1,:) = { mean( dist_c ), mean( [ dist_mon; dist_don ] ),...
    'Control vs ON-treatment PD',...
   P, stats.tstat, stats.df, stats.sd/sqrt(size( [ dist_mon; dist_don ] ,1)) };

[H, P, CI, stats] = ttest2( dist_mon, dist_mof ); % On- vs Off- Medication
statData(2,:) = { mean( dist_mon ), mean(dist_mof ),...
    'Med-ON vs Med-Off',...
   P, stats.tstat, stats.df, stats.sd/sqrt(size( [ dist_mon; dist_mof ] ,1)) };

[H, P, CI, stats] = ttest2( dist_don, dist_dof ); % On- vs Off- DBS
statData(3,:) = { mean( dist_don ), mean(dist_dof ),...
    'DBS-ON vs DBS-Off',...
   P, stats.tstat, stats.df, stats.sd/sqrt(size( [ dist_don; dist_dof ] ,1)) };

[H, P, CI, stats] = ttest2( mean( [dist_mon, dist_mof], 2 ), mean( [dist_don, dist_dof], 2 ) ); % With- vs Without- Electrodes

statData(4,:) = {  mean( mean( [dist_mon, dist_mof]),2 ), mean( mean( [dist_don, dist_dof], 2 ) ),...
    'WO- vs W-Electrode',...
   P, stats.tstat, stats.df, stats.sd/sqrt(size( [ dist_mon; dist_don ] ,1)) };



% Insert table for stats results %
figure(2), clf
f=figure(2);

% Create the column name in cell arrays
cnames = {'mDistance1', 'mDistance2', 'Compared groups','p','t-stat','df','sem'};


% Create the uitable
t = uitable(f,'Data', statData,...
    'ColumnName',cnames,...
    'ColumnWidth',{80, 80, 200, 60, 60, 50, 60});

subplot(1,2,1:2),
pos = get(subplot(1,2,1:2),'position');
title('stats results')
set(subplot(1,2,1:2),'yTick',[])
set(subplot(1,2,1:2),'xTick',[])

set(t,'units','normalized')
set(t,'position',pos)
set(t,'ColumnName',cnames)

set(figure(2),'PaperUnits','inches','PaperPosition',[0 0 9 1.5])
print('-dpng', [DIR.figRevision '//IM_mDistance_centGaze_stats.png'])



end


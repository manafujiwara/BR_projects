function collect_data(DIR, subID, Cond)

%% Single file preproces: Collect behavioral data and eye data for each run / subject.
%% Data is saved into 'Run' structure subjects
whicheye        = Cond.whicheye ;
whicheyeanalyze = Cond.whicheye ;

for isub = 1:length(subID)
    
    if exist([DIR.collectedData '/' subID{isub} '.mat'],'file')
        continue
    end
      
    subName = subID{isub};
    subNum = str2num(subName(2:end));
%     if ~ismember( subNum, [149] ) %127, 130, 131, 149
%         continue
%     end
    
    %% COLLECT BEHAVIOURAL DATA
    DIR.sub = [DIR.dataRaw '/' subName ];
    
    evFiles = dir(  [DIR.sub '/' subName '*_events.txt' ] );
    evFiles = evFiles(find( cellfun(@isempty, regexp( {evFiles(:).name} , '\wprac' ) ) ) );
    evFiles = evFiles(find( ~cellfun(@isempty, regexp( {evFiles(:).name} , '_r(\d)b(\d)' ) ) ) );
    
    % Load each subject once at a time,if you don't specify any return
    % variable to 'load' you'll just get 'Exp' on the workspace
    matFiles = dir(  [DIR.sub '/' subName '_*.mat' ] );
    matFiles = matFiles(find( ~cellfun(@isempty, regexp( {matFiles(:).name} , [subID{isub} '_(\d*)(\.mat)' ] , 'once' ) ) ) );
    
    for imatFiles = 1: length(matFiles)
        
        load([DIR.sub '/' matFiles(imatFiles).name])
        
        idx  = regexp( matFiles(imatFiles).name,'\d+', 'match' );
        ises = str2num(idx{2});
        
        esc  = 0;
        for irun = 1: length(Exp.Run)
            
            for iblk = 1 : length(Exp.Run(irun).block)
                
                % if the session was intermitted, go another session
                if ~exist([ DIR.sub '/' matFiles(imatFiles).name(1:end-4) '_r' num2str(irun) 'b' num2str(iblk) '.txt' ],'file')...
                        || exist([ DIR.sub '/' matFiles(imatFiles).name(1:end-4) '_r' num2str(irun) 'b' num2str(iblk) '_ESC.txt' ],'file')
                    esc = 1; %- break from the blk loop
                    continue %- it used to be 'break' instead, but changed in case the experimenter didn't escape from the session
                end
                
                % BEHAVIOR
                switch Exp.Gral.subjectCategory
                    case 'c'
                        Ses(ises).Run(irun).block(iblk).state = 'nc';
                    case 'p'
                        Ses(ises).Run(irun).block(iblk).state = Exp.Gral.subjectState;
                end
                
                Ses(ises).Run(irun).block(iblk).id                = Exp.Run(irun).block(iblk).id;
                Ses(ises).Run(irun).block(iblk).condition         = Exp.Run(irun).block(iblk).condition;
                Ses(ises).Run(irun).block(iblk).nTrials           = Exp.Run(irun).block(iblk).nTrials;
                Ses(ises).Run(irun).block(iblk).stimFrames        = Exp.Run(irun).block(iblk).stimFrames;
                Ses(ises).Run(irun).block(iblk).interStimInterval = Exp.Run(irun).block(iblk).interStimInterval;
                Ses(ises).Run(irun).block(iblk).rResp             = Exp.Run(irun).block(iblk).rResp;
                Ses(ises).Run(irun).block(iblk).lResp             = Exp.Run(irun).block(iblk).lResp;
                Ses(ises).Run(irun).block(iblk).respValidRate     = Exp.Run(irun).block(iblk).respValidRate;
                Ses(ises).Run(irun).block(iblk).eyeValidRate      = Exp.Run(irun).block(iblk).eyeValidRate;
                Ses(ises).Run(irun).block(iblk).timeTrial         = Exp.Run(irun).block(iblk).TimeTrial;
                
                
                switch  Ses(ises).Run(irun).block(iblk).condition
                    case 'br'
                        mixedRate        = [];
                        mixedLabelDesign = [];
                        borderPos        = [];
                        designFrame      = [];
                        
                    case 'mixed'
                        mixedRate        = Exp.Run(irun).block(iblk).mixedRate;
                        mixedLabelDesign = Exp.Run(irun).block(iblk).mixedLabelDesign;
                        borderPos        = [];
                        designFrame      = [];
                        
                    case 'replay'
                        mixedRate        = [];
                        mixedLabelDesign = [];
                        borderPos        = Exp.Run(irun).block(iblk).replayBorderPos;
                        designFrame      = Exp.Run(irun).block(iblk).replayDesignFrame;
                end
                
                Ses(ises).Run(irun).block(iblk).mixedRate        = mixedRate;
                Ses(ises).Run(irun).block(iblk).mixedLabelDesign = mixedLabelDesign;
                Ses(ises).Run(irun).block(iblk).borderPos        = borderPos;
                Ses(ises).Run(irun).block(iblk).designFrame      = designFrame;
                
                % Preallocate memory for later analysis
                Ses(ises).Run(irun).block(iblk).gazeData   = []; %#ok
                Ses(ises).Run(irun).block(iblk).blinks     = []; %#ok
                Ses(ises).Run(irun).block(iblk).saccades   = []; %#ok
                Ses(ises).Run(irun).block(iblk).fixations  = []; %#ok
                
            end
            
            if esc % if comes to ESC, go to the next session
                break %- break from run loop
            end
            
        end
        
    end
    
    
    %% COLLECT EVENT TIMES
    % Load eye event times
    jses = 0;
    rn_idx = 0;
    endrn = 0;
    for ievFiles = 1 : length(evFiles)
        % Search for the eye data time events
        % Loop through the text file
        
        % OPEN FILE
        evFileName = evFiles(ievFiles).name;
        if exist( [ DIR.sub '/' evFileName(1:end-4) '_ESC.txt' ])
            continue
        end
        
        idx  = regexp(evFileName,'\d+','match');
        ises = str2num(idx{2});
        jrun = str2num(idx{3});
        jblk = str2num(idx{4});
        
        newSes = ises - jses;
        
        fid_evs  = fopen([DIR.sub '/' evFileName]);
        
        while (1)
            
            jses = ises;
            
            %=== Find start of run ===
            % get line by line of text file
            tline = fgetl(fid_evs);
            
            %Skip line if it's an empty one.
            if isempty(tline),  continue, end;
            
            %- Check whether we've reached the end of the file
            if tline == -1, break, end;
            
            if ~isempty(regexp(tline, '(PRAC)','once')),
                continue, end
            
            if ~isempty(regexp(tline, '(END_EXP)','once')) , break, end
            
            % Find next corresponding start of each run (they are indexed)
            if ~isempty(regexp(tline, '(START_RUN)','match'))
                
                %                 if newSes & endrn == 0
                if newSes
                    % if the previous run ended with no event log (it
                    % happens in the last run in a session)
                    rn_idx  = 0;
                end
                
                rn_idx  = rn_idx + 1;
                C = textscan(tline, '%s %n %n');
                
                if rn_idx ~= C{3}
                    msg = ['#run mismatch between event and gaze data. sub:' subID{isub}];
                    error(msg);
                end
                
                bk_idx  = 0; % set block index
                startrn = 1;
                endrn = 0; % "endrn=1" means that the run started but has not been eneded
                
                % go to next line
                continue;
            end
            
            %=== Find start of block ===
            if ~isempty(regexp(tline, 'START_BK','once'))
                
                bk_idx = bk_idx + 1;
                
                msg  = sprintf('Subject: %s, session: %d, Run: %d, block: %d', subID{isub}, ises, jrun, jblk);
                disp(msg)
                
                C = textscan(tline, '%s %n %n'); %- C{1}= event name, C(2)= time, C(3)= num of event
                
                if bk_idx ~= C{3}
                    msg = ['#block mismatch between event and gaze data. sub:' subID(isub)];
                    error(msg);
                end
                
                if ~isempty(regexp(tline, '(REPLAY)','once'))
                    if ~isempty(Ses(ises).Run(rn_idx).block(bk_idx).condition)
                        if isempty( regexp( Ses(ises).Run(rn_idx).block(bk_idx).condition, '(replay)', 'once'))
                            msg = 'condition missmatch';
                            error(msg)
                        end
                    end
                end
                
                Ses(ises).Run(rn_idx).block_evTimes(bk_idx, 1) = C{2}; %#ok % IN SECS
                
                tr_idx = 0; % set trial index
                % go to next line
                continue;
            end
            
            %=== Find start/end of trial & Collect times of trial events ===
            if ~isempty(regexp(tline, 'START_TR','once'))
                tr_idx = tr_idx + 1;
                C      = textscan(tline, '%s %n %n');
                
                if tr_idx ~= C{3}
                    msg = ['#trial mismatch between event and gaze data. sub:' subID(isub)];
                    error(msg);
                end
                
                Ses(ises).Run(rn_idx).block(bk_idx).trial_evTimes(tr_idx, 1) = C{2}; %#ok % IN SECS
                continue;
                
            elseif ~isempty(regexp(tline, 'END_TR','once'))
                C = textscan(tline, '%s %n %n');
                
                if tr_idx ~= C{3}
                    msg = ['#trial mismatch between event and gaze data. sub:' subID(isub)];
                    error(msg);
                end
                
                Ses(ises).Run(rn_idx).block(bk_idx).trial_evTimes(tr_idx, 2) = C{2}; % IN SECS; %#ok
                continue;
                
            end
            
            if ~isempty(regexp(tline, 'END_BK','once'))
                C      = textscan(tline, '%s %n %n');
                
                if bk_idx ~= C{3}
                    msg = ['#block mismatch between event and gaze data. sub:' subID(isub)];
                    error(msg);
                end
                
                Ses(ises).Run(rn_idx).block_evTimes(bk_idx, 2) = C{2}; %#ok % IN SECS
                
                if size( Ses(ises).Run(rn_idx).block,2) == bk_idx
                    break
                else
                    continue;
                end
            end
            
            if ~isempty(regexp(tline, 'END_RUN','once'))
                endrn = 1;
            end
            
            % As the events are at the end of the .txt we need to search for them
            % here separately and add them to the structure
        end
        
    end
    
    %% COLLECT EYE TRACKING DATA
    load ([DIR.cfg '/'  subID{isub} '_Cfg_a_BR.mat'], 'Cfg', 'BR')
    
    PARAMS.doplot.calc                = 0;
    PARAMS.dodownsample               = 0;
    PARAMS.engbertkliegl.VTHRES       = 6; % Velocity threshold
    PARAMS.engbertkliegl.SAMPLING     = 300; % eall.srate;    % Sampling rate
    PARAMS.engbertkliegl.VTYPE        = 2; % Velocity types (2 = using moving average)
    PARAMS.engbertkliegl.SCREEN_RES   = [Cfg.width Cfg.height];% in pix
    PARAMS.engbertkliegl.SCREEN_SIZE  = [Cfg.xDimCm*10 Cfg.yDimCm*10]; % [Exp.Cfg.xDimCm*10 Exp.Cfg.yDimCm*10]; % in mm (LNI)
    PARAMS.engbertkliegl.VIEWING_DIST = 600; % Exp.Cfg.distanceCm*10;% in mm
    PARAMS.engbertkliegl.MININTERSACC = 0.050; %in secs
    PARAMS.engbertkliegl.MINDUR       = 3; % Minimum duration (minimum number of samples needed for a saccade, 3 usually default)
    PARAMS.extra.MINFIXDUR            = 0.006; % IN MS(<- I think this must be in sec) , only for fixation / saccade detection
    % ---------------------------------------------------------------------
    
    gzFiles = dir(  [DIR.sub '/' subID{isub} '*.txt' ] );
    gzFiles = gzFiles(find( cellfun(@isempty, regexp( {gzFiles(:).name} , '\wprac' ) ) ) );
    gzFiles = gzFiles(find( ~cellfun(@isempty, regexp( {gzFiles(:).name} , '_r(\d)b(\d).txt' ) ) ) );
    
    for igzFiles =  1:length(gzFiles)
        samp_idx = 0;
        samps    = nan(1800000,17);
        
        gzFileName = [DIR.sub '/' gzFiles(igzFiles).name(1:end)];
        
        idx  = regexp(gzFileName,'\d+','match');
        ises = str2num(idx{3});
        irun = str2num(idx{4});
        iblk = str2num(idx{5});
        
        fid_gaze = fopen(gzFileName);
        
        % == Explanation of the 17 columns ===
        % 1. time in sec (only sec part)
        % 2. time in msec (only miclosec part)
        % 3. x gaze coordinate of the left eye
        % 4. y gaze coordinate of the left eye
        % 5. x gaze coordinate of the right eye
        % 6. y gaze coordinate of the right eye
        % 7. left eye position - x coordinate
        % 8. left eye position - y coordinate
        % 9. right eye position ? x coordinate
        % 10. right camera eye position - y coordinate
        % 11. left eye validity
        % 12. right eye validity
        % 13. diameter of pupil of the left eye
        % 14. diameter of pupil of the right eye
        % 15. distance of the camera from the left eye
        % 16. distance of the camera from the right eye
        % 17. time in sec(including msec part)
        %   =>  http://hci.javiergs.com/tobii.html
        
        while (1)
            
            % get line by line of text file
            tline = fgetl(fid_gaze);
            
            if isempty(tline), continue, end
            
            if tline == -1, break, end
            
            if length(tline) >= 17 %- if it's eye gaze data
                samp_idx = samp_idx + 1;
                
                C = textscan(tline, '%n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n %n');
                samps(samp_idx, 1:17) = [C{:}];
                
                if mod(samp_idx,10000) == 0
                    msg = sprintf('collect_eyeData. sub: %s, gazeFile: %d/%d, sample: %d', subID{isub}, igzFiles, length(gzFiles),  samp_idx );
                    disp(msg);
                end
            end
        end
        fclose('all'); % clean fids
        
        % samps = samps(1:samp_idx,:);
        % CAUTION!!! have to leave NaN at the end of series of samps to detect
        % blinks and saccades, so you can NOT run the process shown above!!
        
        %% Run Kliegl
        if samp_idx ~= 0
            val = samps(1:samp_idx,11:12);
            [ival q] = find(val~=4);
            
            if length(ival)<= samp_idx/2 % if # of valid points is less than a half of all
                Ses(ises).Run(irun).block(iblk).comment = 'more than the half of gazeData samples are invalid';
                continue
            end
            
            samps(:, 17)    = samps(:,17); %- time in sec.
            samps(:, [3 5]) = samps(:, [3 5]) * Cfg.width; % in pixel
            samps(:, [4 6]) = samps(:, [4 6]) * Cfg.height; % in pixel
            eall.samples    = samps(:, [17 3 4 17 5 6]);
            
            [saccade, blink] = detect_saccades_microsaccades_mod (eall, whicheye, whicheyeanalyze, [], PARAMS);
            
            if size(Ses(ises).Run(irun).block_evTimes,1) >= iblk
                bk_times = Ses(ises).Run(irun).block_evTimes(iblk, :);
                
                Ses(ises).Run(irun).block(iblk).gazeData = samps( samps(:,17) >= bk_times(1) & samps(:,17) <= bk_times(2), :);
                
                if isnan( blink)
                    Ses(ises).Run(irun).block(iblk).blinks = blink;
                else
                    Ses(ises).Run(irun).block(iblk).blinks = blink( blink(:,1) >= bk_times(1) & ...
                        blink(:,2) <= bk_times(2), :);
                end
                
                Ses(ises).Run(irun).block(iblk).macrosac = saccade.macrosac( saccade.macrosac(:,1) >= bk_times(1) & ...
                    saccade.macrosac(:,2) <= bk_times(2), :);
                
                Ses(ises).Run(irun).block(iblk).microsac = saccade.microsac( saccade.microsac(:,1) >= bk_times(1) & ...
                    saccade.microsac(:,2) <= bk_times(2), :);
                
                Ses(ises).Run(irun).block(iblk).allsac = saccade.allsac( saccade.allsac(:,1) >= bk_times(1) & ...
                    saccade.allsac(:,2) <= bk_times(2), :);
                
                for itri = 1 : Ses(ises).Run(irun).block(iblk).nTrials
                    
                    tr_times = Ses(ises).Run(irun).block(iblk).trial_evTimes(itri, :);
                    rawData = samps(samps(:,17) >= tr_times(1) & samps(:,17) <= tr_times(2), :);
                    
                    Ses(ises).Run(irun).block(iblk).trial(itri).gazeData = rawData;
                    
                    if isnan( blink)
                        Ses(ises).Run(irun).block(iblk).trial(itri).blinks = blink;
                        
                    else
                        Ses(ises).Run(irun).block(iblk).trial(itri).blinks = blink(blink(:,1) >= tr_times(1) & ...
                            blink(:,2) <= tr_times(2), :);
                    end
                    
                    Ses(ises).Run(irun).block(iblk).trial(itri).macrosac = saccade.macrosac(saccade.macrosac(:,1) >= tr_times(1) & ...
                        saccade.macrosac(:,2) <= tr_times(2), :);
                    
                    Ses(ises).Run(irun).block(iblk).trial(itri).microsac = saccade.microsac(saccade.microsac(:,1) >= tr_times(1) & ...
                        saccade.microsac(:,2) <= tr_times(2), :);
                    
                    Ses(ises).Run(irun).block(iblk).trial(itri).allsac = saccade.allsac(saccade.allsac(:,1) >= tr_times(1) & ...
                        saccade.allsac(:,2) <= tr_times(2), :);
                    
                end
            end
            clear samp_idx samps
        end
        
    end
    
    %% REMOVE SESSION LAYER
    jrun = 0;
    for ises = 1:size(Ses,2)
        
        if isempty( Ses(ises).Run), continue, end
        
        for irun = 1: size(Ses(ises).Run,2)
            jrun = jrun + 1;
            Run(jrun).block_evTimes = Ses(ises).Run(irun).block_evTimes;
            
            for iblk = 1: size(Ses(ises).Run(irun).block,2)
                Run(jrun).block(iblk).state            = Ses(ises).Run(irun).block(iblk).state;
                Run(jrun).block(iblk).id               = Ses(ises).Run(irun).block(iblk).id;
                Run(jrun).block(iblk).condition        = Ses(ises).Run(irun).block(iblk).condition;
                Run(jrun).block(iblk).nTrials          = Ses(ises).Run(irun).block(iblk).nTrials;
                Run(jrun).block(iblk).stimFrames       = Ses(ises).Run(irun).block(iblk).stimFrames;
                Run(jrun).block(iblk).ISI              = Ses(ises).Run(irun).block(iblk).interStimInterval;
                Run(jrun).block(iblk).oriSesName      = [subID{isub} '-s' num2str(ises) '-r' num2str(irun) '-b' num2str(iblk) ];
                Run(jrun).block(iblk).mixedRate        = Ses(ises).Run(irun).block(iblk).mixedRate;
                Run(jrun).block(iblk).mixedLabelDesign = Ses(ises).Run(irun).block(iblk).mixedLabelDesign;
                Run(jrun).block(iblk).rResp            = Ses(ises).Run(irun).block(iblk).rResp;
                Run(jrun).block(iblk).lResp            = Ses(ises).Run(irun).block(iblk).lResp;
                Run(jrun).block(iblk).respValRate      = Ses(ises).Run(irun).block(iblk).respValidRate;
                Run(jrun).block(iblk).eyeValRate       = Ses(ises).Run(irun).block(iblk).eyeValidRate;
                Run(jrun).block(iblk).timeTrial        = Ses(ises).Run(irun).block(iblk).timeTrial;
                
                if ~isempty(Ses(ises).Run(irun).block(iblk).gazeData)
                    Run(jrun).block(iblk).gazeData         = Ses(ises).Run(irun).block(iblk).gazeData;
                    Run(jrun).block(iblk).blinks           = Ses(ises).Run(irun).block(iblk).blinks;
                    Run(jrun).block(iblk).allsac           = Ses(ises).Run(irun).block(iblk).allsac;
                    Run(jrun).block(iblk).macrosac         = Ses(ises).Run(irun).block(iblk).macrosac;
                    Run(jrun).block(iblk).microsac         = Ses(ises).Run(irun).block(iblk).microsac;
                    Run(jrun).block(iblk).borderPos        = Ses(ises).Run(irun).block(iblk).borderPos;
                    Run(jrun).block(iblk).designFrame      = Ses(ises).Run(irun).block(iblk).designFrame;
                    Run(jrun).block(iblk).trial_evTimes    = Ses(ises).Run(irun).block(iblk).trial_evTimes;
                    Run(jrun).block(iblk).trial            = Ses(ises).Run(irun).block(iblk).trial;
                else
                    Run(jrun).block(iblk).gazeData         = [];
                    Run(jrun).block(iblk).blinks           = [];
                    Run(jrun).block(iblk).allsac           = [];
                    Run(jrun).block(iblk).macrosac         = [];
                    Run(jrun).block(iblk).microsac         = [];
                    Run(jrun).block(iblk).borderPos        = [];
                    Run(jrun).block(iblk).designFrame      = [];
                    Run(jrun).block(iblk).trial_evTimes    = [];
                    Run(jrun).block(iblk).trial            = [];
                end
                
                if isfield(Ses(ises).Run(irun).block(iblk), 'comment')
                    if ~isempty(Ses(ises).Run(irun).block(iblk).comment)
                        Run(jrun).block(iblk).comment = Ses(ises).Run(irun).block(iblk).comment;
                    end
                end
            end
            
        end
    end
    
    
    
    %% GET ALTANATED RESPONSE DATA
    for irun = 1:length(Run)
        nBlocks = length(Run(irun).block);
        
        for iblk = 1: nBlocks
            
            if  iblk > size(Run(irun).block_evTimes,1) ;
                break
            end
            
            nTrials = Run(irun).block(iblk).nTrials;
            tEndTri = Run(irun).block(iblk).stimFrames + Run(irun).block(iblk).ISI;
            
            for itri = 1: nTrials
                
                for ieye = 1:2
                    switch Exp.BR.gratingDir(ieye)
                        case 1
                            resp(:,ieye) = Run(irun).block(iblk).rResp(itri,:);
                        case -1
                            resp(:,ieye) = Run(irun).block(iblk).lResp(itri,:);
                    end
                end
                resp(:,3) = resp(:,1) + resp(:,2); % # pressed button
                
                % in case resp is shorter than it should be (will not
                % happen if it was preallocated in the experiment
                if length(resp) < tEndTri
                    resp(end + 1, :) = nan;
                end
                
                resp( ( resp(:,1) + resp(:,3)) == 2, 4 ) = -1; % single left button press
                resp( ( resp(:,2) + resp(:,3)) == 2, 4 ) = 1; % single right button press
                resp( resp(:,3) == 0, 4 )                = 0; % no button press
                resp( resp(:,3) == 2, 4 )                = 2; % double button press
                resp( isnan( resp(:,3) ), 4)             = -2; % no record
                
                [durInterpolateNo]     = getdur(resp(:,3), 0);
                [durInterpolateDouble] = getdur(resp(:,3), 2);
                durInterpolate         = cat(1, durInterpolateNo, durInterpolateDouble);
                
                if ~isempty(durInterpolate)
                    
                    if durInterpolate(1,2) == length(resp)
                        Run(irun).block(iblk).trial(itri).resp = 'no press';
                        clear resp dur
                        continue
                    end
                    
                    durInterpolate(:,3)    = durInterpolate(:,2) - durInterpolate(:,1);
                    
                    %%% Interpolate missing button press ======================================
                    for idur = 1 : size(durInterpolate,1)
                        tStartDur = durInterpolate(idur,1);
                        tEndDur   = durInterpolate(idur,2);
                        
                        if tStartDur > 1 && tEndDur < tEndTri% if it's not the very beginning
                            
                            respBef = resp(tStartDur - 1, 4);
                            respAft = resp(tEndDur + 1, 4);
                            
                            changeButton = respBef - respAft;
                            
                            if changeButton == 0 % button presses before and after the no/double press duration are the same
                                resp(tStartDur : tEndDur, 5) = respBef;
                                
                            elseif changeButton == 2 || changeButton == 1% right -> no/double -> left
                                tHalfDur = floor((tStartDur + tEndDur)/2);
                                resp(tStartDur : tHalfDur, 5)   = 1;
                                resp(tHalfDur + 1 : tEndDur, 5) = -1;
                                
                            elseif changeButton == -2 || changeButton == -1  % left -> no/double -> right
                                tHalfDur = ceil((tStartDur + tEndDur)/2);
                                resp(tStartDur : tHalfDur, 5)   = -1;
                                resp(tHalfDur + 1 : tEndDur, 5) = 1;
                                
                            end
                            
                        elseif tStartDur > 1 && tEndDur >= tEndTri % if it's the very end
                            respBef = resp(tStartDur - 1, 4);
                            resp(tStartDur:tEndDur,5) = respBef;
                            
                            
                        elseif tStartDur == 1 && tEndDur < tEndTri
                            %%% if it's the very beginning, fill the up the no/double duration
                            %%% only with responce after the duration
                            respAft = resp(tEndDur + 1, 4);
                            resp(1:tEndDur,5) = respAft;
                            
                        elseif tStartDur == 1 && tEndDur >= tEndTri
                            respAft = resp(tEndDur, 4);
                            resp(1:tEndDur,5) = respAft;
                            
                            
                        elseif tStartDur == 0
                            error('tStart is zero, somthing is wrong')
                            
                        end
                    end
                    
                    resp(resp(:,5) == 0, 5) = resp(resp(:,5) == 0, 4); % single button press duration
                    
                else %if there is no point to interpolate
                    resp(:,5) = resp(:,4);
                end
                
                %- resp(:,1) : Left response key (1 or 0)
                %- resp(:,2) : Right response key (1 or 0)
                %- resp(:,3) : # of pressed key
                %- resp(:,4) : response type (0 = no press, 1 = right, -1 = left, 2 = double)
                %- resp(:,5) : interpolated response type (only 1 or -1)
                
                lResp_300Hz = my_interp( resp(:,1)',5);
                rResp_300Hz = my_interp( resp(:,2)',5);
                
                %% get the median for dominant perception duration ========================
                [durL] = getdur(resp(:,5), -1);
                [durR] = getdur(resp(:,5), 1);
                
                durL(:,3) = -1;
                durR(:,3) = 1;
                dur = sortrows(cat(1, durL, durR));
                dur(:,4) = dur(:,2) - dur(:,1) + 1;
                
                %- dur(:,1) : duration start time (in frame)
                %- dur(:,2) : duration end time (in frame)
                %- dur(:,3) : button response (1 = right, -1 = left)
                %- dur(:,4) : length of each duration (in frame)
                
                Run(irun).block(iblk).trial(itri).allResp       = resp;
                Run(irun).block(iblk).trial(itri).altResp       = resp(:,5);
                Run(irun).block(iblk).trial(itri).altResp_300Hz = my_interp(Run(irun).block(iblk).trial(itri).altResp, 5);
                Run(irun).block(iblk).trial(itri).durResp       = dur;
                
                if strcmp(Run(irun).block(iblk).condition,'replay') % only 60-0 condition has "replay" so block=trial in this case
                    Run(irun).block(iblk).trial(itri).borderPos_300Hz = my_interp(Run(irun).block(iblk).borderPos, 5);
                end
                
                clear resp dur
            end
        end
    end
    
    
    %% Save all data to disk
    save([DIR.collectedData '/' subName '.mat'],'Run', 'PARAMS')
    clear Exp Ses Run Cfg
end

function [dur] = getdur(data, resp)

%%%
% data: vector of the continuous response which you want to get the duration
% resp: a number or a vector containds numbers if which you want to get the
% duration.

iresp = find(data == resp);
iswitchResp = find(diff(iresp) ~= 1);

if ~isempty(iresp)
    dur(1,1) = iresp(1);
    for iiswitchResp = 1: length(iswitchResp)
        dur(iiswitchResp, 2) = iresp(iswitchResp(iiswitchResp));
        dur(iiswitchResp + 1, 1) = iresp(iswitchResp(iiswitchResp) + 1);
        
        if iiswitchResp == length(iswitchResp)
            dur(end, 2) = iresp(end); % as we can not detect the end of the last duration from diff
        end
    end
    
    if ~isempty(iresp) && isempty(iswitchResp)
        dur(1,1) = iresp(1);
        dur(1,2) = iresp(end);
    end
    
    if length(dur) == 1 || isempty(iresp)
        dur = [];
    end
    
else
    dur = [];
end

function y = my_interp(vector, factor)

y = [];
for m = 1 : length(vector)
    
    x = vector (m);
    y = cat(2, y, repmat(x, 1, factor));
    
end

y = y';

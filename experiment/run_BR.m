function run_BR
% dbstop if error
%% Modified by MF 16 Jul 2015

%% Modified by MF 11 Jun 2014
% based on run_SMF, MF added Binocular ribalry part. Also, edited action
% when escape key was pressed (line 254~);
% 1) save(sprintf('%s_ESC', fileName));
% -> save([fileName(1:end-4) '%s_ESC' fileName(end-3:end)])
% 2) break -> sca, return.

%% Modified by NT 13 Mar 11
% based on runConf_1k_SFM_new, which contains commands related to
% communication with eyelink.

% Goals
% 1. to present SFM stimuli to parkinson's disease patient.
% - preliminary studies show that their percept switches very quickly
% - test this by passive viewing (constant stimulation) with button press,
% or, by interleaved presentation, which might elicit stronger OKN at the
% onset of the stimuli.

% Pre-requisites:
% Add Psychophysics/ to the path to use 'initializeScreen.m'

% Parameters
% SFM is a structure that contains basic parameters that specify the
% spatiotemporal parameters for structure from motion
%
% T is a structure that specifies parameters that are manipulated in each
% trial
% T.direction(1) - rotation direction for top
% T.direction(2) - rotation direction for bottom
% T.disparity(1) - disparity for top
% T.disparity(2) - disparity for top

addpath('./aux_files/');

%%% For mac in oculomotor room %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('/../../../Psychophysics/'));
addpath(genpath('/Applications/Psychtoolbox'));
addpath(genpath('/Users/lisandrk/no-ppc/lib/matlab')); % to the latest version of t2t


%%% For local machine %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% addpath(genpath('/Applications/Psychtoolbox'));
% PsychJavaTrouble;

try
    % check keyboards
    Exp = setKeypads;
    checkKeypads;
    
    %% Initialize screen
    % PsychJavaTrouble; % Check there are no problems with Java
    Exp.Cfg.SkipSyncTest = 0; %This should be '0' on a properly working NVIDIA video card. '1' skips the SyncTest.
    Exp.Cfg.AuxBuffers   = 1; % '0' if no Auxiliary Buffers available, otherwise put it into '1'.
    
    % Check for OpenGL compatibility
    AssertOpenGL;
    Screen('Preference','SkipSyncTests', Exp.Cfg.SkipSyncTest);
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  PARAMETERS FOR EXPERIMENTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Exp.Cfg.WinSize    = [];  %Empty means whole screen
    Exp.Cfg.xDimCm     = 51; %Length in cm of the screen in X, Tobii
    Exp.Cfg.yDimCm     = 29; %Length in cm of the screen in Y, Tobii
    Exp.Cfg.distanceCm = 74; %Viewing distance
    
    %% PRACTICE ============================================================
    run_prac = 1;
    
    if run_prac
        prac_nBlocks      = 4;
        prac_blk_id       = {'30-0'    , '2-1',       '30-0'      '2-1'}; %- avairable only 30-0 or 2-1
        prac_cond         = {'d_replay', 'mixed',     'br',       'mixed' };
        prac_nTrials      = [1,           10,          1,          10 ];% 30 20]; % trials for each block
        stimFrames        = [1800,        120,         1800,       120];% 60 120]; % duration of stimulation -in frames-
        interStimInterval = [0,           60,          0,          60 ]; % duration of intervals between stimulation -in frame
        
        for iblk = 1 : prac_nBlocks
            Prac.Run(1).block(iblk).condition         = prac_cond{iblk};
            Prac.Run(1).block(iblk).id                = prac_blk_id{iblk};
            Prac.Run(1).block(iblk).nTrials           = prac_nTrials(iblk);
            Prac.Run(1).block(iblk).stimFrames        = stimFrames(iblk);
            Prac.Run(1).block(iblk).interStimInterval = interStimInterval(iblk);
            %                 Prac.Run(irun).block(iblk).disparity         = 0;
            Prac.Run(1).block(iblk).rResp = nan(prac_nTrials(iblk),stimFrames(iblk)+interStimInterval(iblk));
            Prac.Run(1).block(iblk).lResp = nan(prac_nTrials(iblk),stimFrames(iblk)+interStimInterval(iblk));
        end
        
    end
    
    %% MAIN ================================================================
    % subject information %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subjectCategory = {'p', 'c'};
    subjectStateDlg = {'Med ON', 'Med OFF', 'DBS ON + Med', 'DBS OFF + Med', 'DBS ON + NO Med', 'DBS OFF + NO Med'};
    subjectStateClt = {'mon', 'mof', 'donmon', 'dofmon', 'donmof', 'dofmof'};
    
    dialogName = 'Subject Information';
    prompt     = {'Subject Name:', 'Subject Number:'};
    numlines   = 1;
    answer     = inputdlg(prompt,dialogName,numlines);
    Exp.Gral.subjectName   = answer{1};
    Exp.Gral.subjectNumber = str2num(answer{2});
    
    [s, v] = listdlg('PromptString','Subject Category',...
        'SelectionMode','single',...
        'ListString', subjectCategory );
    Exp.Gral.subjectCategory = subjectCategory{s};
    
    if s == 1 % if the subject is a patient
        [s, v] = listdlg('PromptString','Subject State',...
            'SelectionMode','single',...
            'ListString', subjectStateDlg );
        Exp.Gral.subjectState = subjectStateClt{s};
    end
    
    % global constants %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Exp = initializeScreen(Exp);
    
    Exp.Cfg.WinColor   = Exp.Cfg.Color.gray; % empty for the middle gray of the screen.
    Exp.Cfg.drawEye    = 0; % Draw eye online on the screen
    Exp.Cfg.trackEye   = 1; % 1 : collect eye tracking data
    
    Exp.BR.gratingDir   = [ -1, 1 ]; % 1)left eye, 2)right eye. -1 = towards left, 1= towards right
    Exp.BR.gratingColor = [ Exp.Cfg.Color.green; Exp.Cfg.Color.red ];
    Exp.BR.replay       = 1; % run replay for each block(=1) or not(=0) % haven't had an option to run just a few blocks
    
    nRuns               = 4;
    nBlocks             = 4; % number of blocks for the whole experiment
    blk_id              = {'60-0',  '60-0',    '2-1',   '2-1' }; %ALWAYS BR & REPLAY WILL BE COMBINED (so, # of blocks will be doubled)
    blk_cond            = {'br',    'replay',  'mixed', 'mixed'};
    nTrials             = [ 1,       1,        20,       20  ];% ; % trials for each block
    stimFrames          = [ 3600,    3600,     120,      120 ];% ; % duration of stimulation -in frames-
    interStimInterval   = [ 0,       0         60,       60  ]; % duration of intervals between stimulation -in frames-
    eye_samples         = ( stimFrames + interStimInterval ) * 5; % number of eye samples to collect during stimulation + blank in each trial
    
    %------example---------------------------------------------------------
    %     nRuns               = 4;
    %     nBlocks             = 2; % number of blocks for the whole experiment
    %     blk_id              = {'60-0', '2-1'}; %ALWAYS BR & REPLAY WILL BE COMBINED (so, # of blocks will be doubled)
    %     nTrials             = [1,      20  ]; % trials for each block
    %     stimFrames          = [3600,   120 ]; % duration of stimulation -in frames-
    %     interStimInterval   = [0,      60  ]; % duration of intervals between stimulation -in frames-
    %     eye_samples         = ( stimFrames + interStimInterval ) * 5; % number of eye samples to collect during stimulation + blank in each trial
    % ---------------------------------------------------------------------
    if sum(bsxfun(@times, (stimFrames + interStimInterval), nTrials) ~= 3600)
        error('length of trial and number of trial mismatches')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %% Define Parameters
    Msg.main.startRun = 'START_RUN';
    Msg.main.endRun   = 'END_RUN';
    Msg.main.startBlk = 'START_BK';
    Msg.main.endBlk   = 'END_BK';
    Msg.main.startTri = 'START_TR';
    Msg.main.endTri   = 'END_TR';
    
    Msg.prac.startRun = 'START_RUN_PRAC';
    Msg.prac.endRun   = 'END_RUN_PRAC';
    Msg.prac.startBlk = 'START_BK_PRAC';
    Msg.prac.endBlk   = 'END_BK_PRAC';
    Msg.prac.startTri = 'START_TR_PRAC';
    Msg.prac.endTri   = 'END_TR_PRAC';
    
    
    % subjectID
    subjectID = num2str(Exp.Gral.subjectNumber); %<- needs to be 2 latters
    
    % check input parameter
    if length(subjectID) ~= 2
        subjectID = [num2str(0) subjectID];
        %         error('Subject ID must have two characters');
    end
    
    % Set experiment paths
    [Exp, fileName, ~, sessionNum] = setPaths(Exp, subjectID);
    
    % set Binocular Rivalry parameters
    Exp.BR.visibleSize        = [270 270];%round([22.2 * Exp.Cfg.pixelsPerDegree, 30 * Exp.Cfg.pixelsPerDegree]); % - should be smaller than Exp.BR.frameInsRectSize
    Exp.BR.FsCyclesPerDeg     = 0.27;%8 / (Exp.BR.visibleSize(2) * Exp.Cfg.degPerPixel); %spatial frequency : make it to have 8 cycles in a stimulus
    Exp.BR.FsCyclesPerPix     = Exp.BR.FsCyclesPerDeg * Exp.Cfg.degPerPixel; %spatial frequency. (cycle/pixel)
    Exp.BR.FsRadPerPix        = Exp.BR.FsCyclesPerPix * 2 * pi;    % frequency in radians.
    Exp.BR.speedDegPerSec     = 22.3;%/ 30 * (Exp.BR.visibleSize(2) * Exp.Cfg.degPerPixel); % degree/sec
    Exp.BR.speedCyclesPerSec  = Exp.BR.speedDegPerSec * Exp.BR.FsCyclesPerDeg;
    Exp.BR.gratingPixPerCycle = 1 / Exp.BR.FsCyclesPerPix;
    Exp.BR.gratingSizePix     = Exp.BR.gratingPixPerCycle * 10;
    
    Exp.BR.replayOverlapWidthDeg = 14; %- in degree
    
    % create grating texture
    for ieye = 1:2
        switch Exp.BR.gratingDir(ieye)
            case 1
                Exp.BR.gratingAng(ieye) = 180;
            case -1
                Exp.BR.gratingAng(ieye) = 0;
        end
    end
    
    z = meshgrid(0 : Exp.BR.gratingSizePix - 1, 1 );
    Exp.BR.grating    = Exp.Cfg.Color.inc * cos (Exp.BR.FsRadPerPix * z) + Exp.Cfg.Color.gray;
    Exp.BR.gratingTex = Screen('MakeTexture', Exp.Cfg.win, Exp.BR.grating, [], 1);
    
    Exp.BR.waitFrames       = 1;
    Exp.BR.waitDuration     = Exp.BR.waitFrames * Exp.Cfg.MonitorFlipInterval;
    Exp.BR.shiftPixPerFrame = Exp.BR.speedDegPerSec / Exp.Cfg.FrameRate * Exp.Cfg.pixelsPerDegree;
    %     Exp.BR.shiftPerFrameReplay = Exp.BR.speedCyclesPerSec * Exp.BR.gratingPixPerCycle * Exp.BR.waitDuration * 2 * pi / 60;
    Exp.BR.shiftCyclePerFrameReplay = Exp.BR.speedCyclesPerSec / 60;
    
    
    Exp.BR.rectPosFromCenter = 350;
    Exp.BR.frameWidth        = 30;
    Exp.BR.frame             = rand(Exp.BR.visibleSize + Exp.BR.frameWidth) * 255;
    Exp.BR.frame(15: Exp.BR.visibleSize(1) + Exp.BR.frameWidth - 15, 15 : Exp.BR.visibleSize(2) + Exp.BR.frameWidth - 15) = 0;
    
    %% make frames
    Exp.BR.frametex = Screen('MakeTexture', Exp.Cfg.win, Exp.BR.frame);
    
    % left frame /grating position
    Exp.BR.framePos(:,1) = CenterRectOnPointd([0 0 size(Exp.BR.frame)],...
        Exp.Cfg.centerX - Exp.BR.rectPosFromCenter, Exp.Cfg.centerY); % - for frame InsRect (L)
    
    Exp.BR.gratingDestRectPos(:,1) = CenterRectOnPointd([0 0 Exp.BR.visibleSize],...
        Exp.Cfg.centerX - Exp.BR.rectPosFromCenter, Exp.Cfg.centerY); % - for grating (L)
    
    % right frame
    Exp.BR.framePos(:,2) = CenterRectOnPointd([0 0 size(Exp.BR.frame)],...
        Exp.Cfg.centerX + Exp.BR.rectPosFromCenter, Exp.Cfg.centerY); % - for frame InsRect (R)
    
    Exp.BR.gratingDestRectPos(:,2) = CenterRectOnPointd([0 0 Exp.BR.visibleSize],...
        Exp.Cfg.centerX + Exp.BR.rectPosFromCenter, Exp.Cfg.centerY); % - for grating (R)
    
    
    % Initialize random number generators to subject dependent reconstructible
    % values
    Exp.seed = 1e5 * sessionNum + 100 * subjectID(1) + subjectID(2);
    rand('state', Exp.seed);
    randn('state', Exp.seed);
    
    
    %% Preallocate all trials.
    for irun = 1 : nRuns
        for iblk = 1 : nBlocks
            Exp.Run(irun).block(iblk).id                = blk_id{iblk};
            Exp.Run(irun).block(iblk).condition         = blk_cond{iblk};
            Exp.Run(irun).block(iblk).nTrials           = nTrials(iblk);
            Exp.Run(irun).block(iblk).stimFrames        = stimFrames(iblk);
            Exp.Run(irun).block(iblk).interStimInterval = interStimInterval(iblk);
        end
    end
    
    % Set the order of blocks
    %     Exp.randomize = 1;
    %     if Exp.randomize
    %         for irun = 1:nRuns
    %             Exp_rndblock_idxs   = randperm(nBlocks);
    %             Exp.Run(irun).block = Exp.Run(irun).block(Exp_rndblock_idxs);
    %         end
    %     end
    
    %% Present the two frames for binocular fusion
    % Draw frames for binocular fusion
    Screen('FillRect',  Exp.Cfg.win, Exp.Cfg.WinColor);
    Screen('DrawTextures', Exp.Cfg.win, Exp.BR.frametex, [], Exp.BR.framePos);
    
    Screen(Exp.Cfg.win, 'DrawText', 'Press "space"', Exp.Cfg.centerX - 450, Exp.Cfg.centerY-30, Exp.Cfg.Color.white);
    Screen(Exp.Cfg.win, 'DrawText', 'if two frames', Exp.Cfg.centerX - 450, Exp.Cfg.centerY, Exp.Cfg.Color.white);
    Screen(Exp.Cfg.win, 'DrawText', 'are merged.', Exp.Cfg.centerX - 450, Exp.Cfg.centerY+30, Exp.Cfg.Color.white);
    
    Screen(Exp.Cfg.win, 'DrawText', 'Press "space"', Exp.Cfg.centerX + 250, Exp.Cfg.centerY-30, Exp.Cfg.Color.white);
    Screen(Exp.Cfg.win, 'DrawText', 'if two frames', Exp.Cfg.centerX + 250, Exp.Cfg.centerY, Exp.Cfg.Color.white);
    Screen(Exp.Cfg.win, 'DrawText', 'are merged.', Exp.Cfg.centerX + 250, Exp.Cfg.centerY+30, Exp.Cfg.Color.white);
    
    Screen('Flip', Exp.Cfg.win,  [], Exp.Cfg.AuxBuffers);
    
    ContOnlyIfSpacePress
    
    WaitSecs(0.5);
    
    %% Initialize Tobii (Calibration)
    hostName = '169.254.7.248';
    portName = '4455';
    res      = [Exp.Cfg.width Exp.Cfg.height];
    
    if Exp.Cfg.trackEye
        % Exp.Tobii = TobiiInit(hostName, portName, Exp.Cfg.win, res, Exp);
        Exp.Tobii = TobiiInit_BR(hostName, portName, Exp.Cfg.win, res, Exp);
    end
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PRACTICE LOOP HERE                                                 %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% planed replay stimulus %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if run_prac
        irun = 1;
        fileNamePrac = [fileName(1:end-4) '_prac' fileName(end-3:end)];
        
        % ready to go?
        Screen('FillRect',  Exp.Cfg.win, Exp.Cfg.WinColor);
        Screen('DrawTextures', Exp.Cfg.win, Exp.BR.frametex, [], Exp.BR.framePos);
        Screen(Exp.Cfg.win, 'DrawText', 'Practice Exp', Exp.Cfg.centerX - 450, Exp.Cfg.centerY-60, Exp.Cfg.Color.white);
        Screen(Exp.Cfg.win, 'DrawText', 'Experimenter:', Exp.Cfg.centerX - 450, Exp.Cfg.centerY-30, Exp.Cfg.Color.white);
        Screen(Exp.Cfg.win, 'DrawText', 'press "space" ', Exp.Cfg.centerX - 450, Exp.Cfg.centerY, Exp.Cfg.Color.white);
        Screen(Exp.Cfg.win, 'DrawText', 'to start session', Exp.Cfg.centerX - 450, Exp.Cfg.centerY+30, Exp.Cfg.Color.white);
        
        Screen(Exp.Cfg.win, 'DrawText', 'Practice Exp', Exp.Cfg.centerX + 250, Exp.Cfg.centerY-60, Exp.Cfg.Color.white);
        Screen(Exp.Cfg.win, 'DrawText', 'Experimenter:', Exp.Cfg.centerX + 250, Exp.Cfg.centerY-30, Exp.Cfg.Color.white);
        Screen(Exp.Cfg.win, 'DrawText', 'press "space"', Exp.Cfg.centerX + 250, Exp.Cfg.centerY, Exp.Cfg.Color.white);
        Screen(Exp.Cfg.win, 'DrawText', 'to start session', Exp.Cfg.centerX + 250, Exp.Cfg.centerY+30, Exp.Cfg.Color.white);
        
        Screen('Flip', Exp.Cfg.win, [], Exp.Cfg.AuxBuffers);
        %         HideCursor;
        ListenChar(0);
        
        % Continue only if space is pressed
        ContOnlyIfSpacePress
        
        % Trigger the beginning of the block
        if Exp.Cfg.trackEye
            mesg = Msg.prac.startRun;
            runStart = talk2tobii('EVENT', mesg, irun);
            WaitSecs(0.2)
        end
        
        %         % Initialize Tobii (Calibration)
        %         if Exp.Cfg.trackEye == 1
        %             % Exp.Tobii = TobiiInit(hostName, portName, Exp.Cfg.win, res, Exp);
        %             message = 'Initializing Tobii...';
        %             disp(message);
        %             Exp.Tobii = TobiiCalib_BR(hostName, portName, Exp.Cfg.win, res, Exp);
        %         end
        
        Prac.Gral = Exp.Gral;
        Prac.Cfg  = Exp.Cfg;
        Prac.BR   = Exp.BR;
        Prac.addParams = Exp.addParams;
        Prac.BR.replayOverlapWidthPix = Prac.BR.replayOverlapWidthDeg * Prac.Cfg.pixelsPerDegree;
        
        
        for iblk = 1: length(Prac.Run(irun).block)
            
            switch Prac.Run(irun).block(iblk).condition
                
                case 'd_replay'
                    %% set parameters ---------------------------------------------
                    switch Prac.Run(irun).block(iblk).id
                        
                        case '30-0'
                            replayDesign{1}(:,4) = [5, 5, 3, 1, 2, 4, 1, 3, 4, 2]' * 60; % in frame
                            replayDesign{1}(:,3) = [-1, 1, -1, 1, -1, 1, -1, 1, -1, 1]'; % direction of grating
                            replayDesign{1}(1,1) = 1;
                            replayDesign{1}(1,2) = replayDesign{1}(1,4);
                            
                            Prac.Run(irun).block(iblk).replayMedianFrac{1} = 0.8 * ones(1, length(replayDesign{1})); % just use 0.8 as a fixed median fractuation for all trials and durations
                            
                            for idur = 2: size( replayDesign{1}, 1 )
                                replayDesign{1}(idur,1) = replayDesign{1}(idur-1,2) + 1;
                                replayDesign{1}(idur,2) = sum(replayDesign{1}(1:idur,4));
                            end
                            durAllTri = replayDesign{1};
                            
                            switch durAllTri(1,3)
                                case 1 % in case of the grating is going to the right, swiping will start from the left side
                                    Prac.Run(irun).block(iblk).replayBorderPos(1,1) = Prac.BR.visibleSize(1) +  Prac.BR.replayOverlapWidthPix/2 ;
                                case -1 % in case of the grating is going to the left, swiping will start from the right side
                                    Prac.Run(irun).block(iblk).replayBorderPos(1,1) = 0 - Prac.BR.replayOverlapWidthPix/2;
                            end
                            Prac.Run(irun).block(iblk).replayDesignFrame{1} = durAllTri;
                            Prac.Run(irun).block(iblk).replayMedianDuration = median (replayDesign{1}(:,4));
                    end
                    
                    %% run practice -----------------------------------------------
                    [Prac KeyFlag] = prac_run_BR_blk_replay(Prac, fileNamePrac, irun, iblk, Msg.prac);
                    if KeyFlag, break, end % break from block loop
                    clear replayDesign durAllTri
                    
                    
                    
                case 'br'
                    %%% BR stimulus %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    [Prac KeyFlag] = prac_run_BR_blk_br(Prac, fileNamePrac, irun, iblk, Msg.prac);
                    
                case 'mixed'
                    switch Prac.Run(irun).block(iblk).id
                        case '2-1'
                            
                            Prac.Run(irun).block(iblk).mixedRate = [ 0.5, 0.25, 0.25 ]; %- ratio of trial for 1)br, 2)red/right unambiguous, 3)green/left unambiduous
                            
                            [Prac KeyFlag] =  run_BR_blk_mixed(Prac, fileName, irun, iblk, Msg.prac);
                            
                            if KeyFlag, break, end % break from block loop
                            
                    end
            end
            
        end
        
        % Trigger the beginning of the block
        if Exp.Cfg.trackEye
            mesg =  Msg.prac.endRun;
            runStart = talk2tobii('EVENT', mesg, irun);
            WaitSecs(0.2)
        end
        
        fprintf('saving matlab file...');
        save(fileName, 'Prac');
        fprintf(['Matfile for Prac Run ' num2str(irun) ' has been saved' '\n']);
        
        fprintf('saving eye data ...');
        talk2tobii('SAVE_DATA', [fileNamePrac(1:end-4) '.txt'], ...
            [fileNamePrac(1:end-4) '_events.txt'], 'TRUNK');
        fprintf(['Eye data for Prac Run ' num2str(irun) ' has been saved' '\n']);
        
        
        %% End of practice: stop tracking eyes and save Prac
        
        % plot the result of Practice session
        disp('CHECKING PRACTICE PERFORMANCE>>>>>')
        check_prac_run_BR
        
        %         ListenChar(1);
        
    end
    
    
    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  MAIN LOOP HERE                                                   %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for irun = 1 : nRuns
        % Trigger the beginning of the run
        if Exp.Cfg.trackEye
            mesg = Msg.main.startRun;
            runStart = talk2tobii('EVENT', mesg, irun);
            WaitSecs(0.2)
        end
            
        % Tobii Calibration if needed
        pos = [2500 430 400 100];%- [winLeftGap bottom winWide winHigh];
        if Exp.Cfg.trackEye
            calib = 0;
            while (1)
                calib = menu2(pos,'Caliblation?','Yes','No');
                if calib ~= 0, break,end
            end
            switch calib
                case 1
                    % Exp.Tobii = TobiiInit(hostName, portName, Exp.Cfg.win, res, Exp);
                    Exp.Run(irun).TobiiCalib = TobiiCalib_BR(hostName, portName, Exp.Cfg.win, res, Exp);
                    WaitSecs(0.1)
                case 2
                    Exp.Run(irun).TobiiCalib = nan;
            end
        end
        
        for iblk = 1 : nBlocks
            
            eyeFileName   = [fileName(1:end-4) '_r' num2str(irun) 'b' num2str(iblk) '.txt' ];
            eventFileName = [fileName(1:end-4) '_r' num2str(irun) 'b' num2str(iblk) '_events.txt'];
            eyeFileNameESC = [eyeFileName(1:end-4) '_ESC.txt'];
            eventFileNameESC = [eventFileName(1:end-4) '_ESC.txt'];
            
            % Tobii draw eyes if needed
            if Exp.Cfg.trackEye
                drawEyes = 0;
                while (1)
                    drawEyes = menu2(pos,'Draw eyes?','Yes','No');
                    if drawEyes ~= 0,
                        Exp.Run(irun).block(iblk).TobiiDrawEye = drawEyes-1;
                        break
                    end
                end
                
                if drawEyes == 1
                        Exp.Run(irun).TobiiDrawEyes = 1;
                        
                        % Display eyes on the display
                        flagNotBreak = 0;
                        disp('Press ''SPACE'' to continue');
                        
                        while ~flagNotBreak
                            eyeTrack = talk2tobii('GET_SAMPLE');
                            DrawEyes_BR(Exp.Cfg.win, eyeTrack(9), eyeTrack(10),eyeTrack(11), eyeTrack(12),eyeTrack(8), eyeTrack(7),Exp);
                            
                            if( IsKey(Exp.addParams.spaceKey) )
                                flagNotBreak = 1;
                                if( flagNotBreak )
                                    break;
                                end
                            end
                        end
                        WaitSecs(0.1)
                end
            end
            
            % ready to go?
            Screen('FillRect',  Exp.Cfg.win, Exp.Cfg.WinColor);
            Screen('DrawTextures', Exp.Cfg.win, Exp.BR.frametex, [], Exp.BR.framePos);
            Screen(Exp.Cfg.win, 'DrawText', 'Main Exp', Exp.Cfg.centerX - 450, Exp.Cfg.centerY-60, Exp.Cfg.Color.white);
            Screen(Exp.Cfg.win, 'DrawText', 'Experimenter:', Exp.Cfg.centerX - 450, Exp.Cfg.centerY-30, Exp.Cfg.Color.white);
            Screen(Exp.Cfg.win, 'DrawText', 'press "space" ', Exp.Cfg.centerX - 450, Exp.Cfg.centerY, Exp.Cfg.Color.white);
            Screen(Exp.Cfg.win, 'DrawText', 'to start session', Exp.Cfg.centerX - 450, Exp.Cfg.centerY+30, Exp.Cfg.Color.white);
            
            Screen(Exp.Cfg.win, 'DrawText', 'Main Exp', Exp.Cfg.centerX + 250, Exp.Cfg.centerY-60, Exp.Cfg.Color.white);
            Screen(Exp.Cfg.win, 'DrawText', 'Experimenter:', Exp.Cfg.centerX + 250, Exp.Cfg.centerY-30, Exp.Cfg.Color.white);
            Screen(Exp.Cfg.win, 'DrawText', 'press "space"', Exp.Cfg.centerX + 250, Exp.Cfg.centerY, Exp.Cfg.Color.white);
            Screen(Exp.Cfg.win, 'DrawText', 'to start session', Exp.Cfg.centerX + 250, Exp.Cfg.centerY+30, Exp.Cfg.Color.white);
            
            Screen('Flip', Exp.Cfg.win, [], Exp.Cfg.AuxBuffers);
            %     HideCursor;
            ListenChar(0);
            
            % Continue only if space is pressed
            ContOnlyIfSpacePress
            if exist('KeyFlag', 'var') && KeyFlag==1, break, end
            
            WaitSecs(0.2);
            
            
            switch Exp.Run(irun).block(iblk).id
                case '60-0' %- 60-0 condition has BR(=normal) trials and replay condition separately in different block
                    
                    switch Exp.Run(irun).block(iblk).condition
                        case 'br'
                            [Exp KeyFlag] =  run_BR_blk_br(Exp, fileName, irun, iblk, Msg.main);
                            
                            if KeyFlag, break, end % break from block loop
                            
                        case 'replay' % mimic the previous normal block
                            if  ~isempty(find(Exp.Run(irun).block(iblk-1).rResp,1)) ...
                                    || ~isempty(find(Exp.Run(irun).block(iblk-1).lResp,1))
                                
                                [Exp KeyFlag] = run_BR_blk_replay(Exp, fileName, irun, iblk, Msg.main);
                                
                                if KeyFlag, break, end % break from block loop
                            end
                    end
                    
                case '2-1' %- 2-1 condition has BR(=ambiguous) and unambiguous trials mixed in the same block
                    
                    Exp.Run(irun).block(iblk).mixedRate = [ 0.5 0.25 0.25 ]; %- ratio of trial for 1)br, 2)red/right unambiguous, 3)green/left unambiduous
                    
                    [Exp KeyFlag] =  run_BR_blk_mixed(Exp, fileName, irun, iblk, Msg.main);
                    
                    if KeyFlag, break, end % break from block loop
                    
            end
            
            % save eye tracking files each block
            if  Exp.Cfg.trackEye
              
                talk2tobii('STOP_AUTO_SYNC');
                talk2tobii('STOP_RECORD');
                talk2tobii('STOP_TRACKING');
                talk2tobii('SAVE_DATA', eyeFileName, eventFileName , 'TRUNK');
                talk2tobii('CLEAR_DATA');
                
                save(fileName,'Exp')
                
                talk2tobii('START_TRACKING');
                WaitSecs(0.2)
                
                talk2tobii('RECORD');
                WaitSecs(0.2)
                
                % TET_API_AUTOSYNCED
                talk2tobii('START_AUTO_SYNC')
                WaitSecs(2)
                
                [eyeValidRate]  = check_eyeValidity_blk(fileName, eyeFileName, eventFileName, irun, iblk);
                [respValidRate] = check_respValidity_blk(fileName, irun, iblk);
                
                Exp.Run(irun).block(iblk).respValidRate = respValidRate;
                Exp.Run(irun).block(iblk).eyeValidRate  = eyeValidRate;
                save(fileName,'Exp')
                
                % present messagebox untill it's ready to continue
                title = 'Data Validity';
                string = {['eye gaze data: (L)' num2str(eyeValidRate(1)) '%, (R)' num2str(eyeValidRate(2)) '%'],...
                    ['button resp     : (L)' num2str(respValidRate(1)) '%, (R)' num2str(respValidRate(2)) '%'],...
                    [],...
                    'Press OK when you are ready to continue.'};
                
                h = msgbox(string,title,'modal');
                ah = get( h, 'CurrentAxes' );
                ch = get( ah, 'Children' );
                set( ch, 'FontSize', 20 );
                set(h,'Position',[2500 430 500 150])
                uiwait(h);
                
            else
                save(fileName,'Exp');
                [respValidRate] = check_respValidity_blk(fileName, irun, iblk);
                Exp.Run(irun).block(iblk).respValidRate = respValidRate;
                save(fileName,'Exp')
            end
            
        end
        
        % TAKE A PAUSE IN BETWEEN RUNS
        if Exp.Cfg.trackEye
            mesg = Msg.main.endRun;
            runEnd = talk2tobii('EVENT', mesg, irun);
            WaitSecs(0.2)
        end
        
        Screen('FillRect',  Exp.Cfg.win, Exp.Cfg.WinColor);
        Screen('DrawTextures', Exp.Cfg.win, Exp.BR.frametex, [], Exp.BR.framePos);
        
        Screen(Exp.Cfg.win, 'DrawText', ['End of Run ' num2str(irun) ' of ' num2str(length(Exp.Run)) ], Exp.Cfg.centerX - 450, Exp.Cfg.centerY-50, Exp.Cfg.Color.white);
        Screen(Exp.Cfg.win, 'DrawText',  'Take a pause.', Exp.Cfg.centerX - 450, Exp.Cfg.centerY-25, Exp.Cfg.Color.white);
        Screen(Exp.Cfg.win, 'DrawText',  'Press "space"', Exp.Cfg.centerX - 450, Exp.Cfg.centerY, Exp.Cfg.Color.white);
        Screen(Exp.Cfg.win, 'DrawText',  'to continue.', Exp.Cfg.centerX - 450, Exp.Cfg.centerY+25, Exp.Cfg.Color.white);
        
        Screen(Exp.Cfg.win, 'DrawText', ['End of Run ' num2str(irun) ' of ' num2str(length(Exp.Run)) ], Exp.Cfg.centerX + 250, Exp.Cfg.centerY-50, Exp.Cfg.Color.white);
        Screen(Exp.Cfg.win, 'DrawText', 'Take a pause.', Exp.Cfg.centerX + 250, Exp.Cfg.centerY-25, Exp.Cfg.Color.white);
        Screen(Exp.Cfg.win, 'DrawText', 'Press "space"', Exp.Cfg.centerX + 250, Exp.Cfg.centerY, Exp.Cfg.Color.white);
        Screen(Exp.Cfg.win, 'DrawText', 'to continue.', Exp.Cfg.centerX + 250, Exp.Cfg.centerY+25, Exp.Cfg.Color.white);
        
        Screen('flip',Exp.Cfg.win);
        
        % Continue only after pressing the space bar @main keyboard (by
        % experimenter) or number keys @numeric keypads (subject)
        ContOnlyIfSpacePress
        
        if KeyFlag, break, end
        
        WaitSecs(0.5)
        
        if irun ~= length(Exp.Run) % save data in case of getting an error
            fprintf('saving matlab file...');
            save(fileName, 'Exp');
            fprintf(['Matfile for Run ' num2str(irun) ' has been saved' '\n']);
            
            fprintf('saving eye data ...');
            talk2tobii('SAVE_DATA', [fileName(1:end-4) '.txt'], ...
                [fileName(1:end-4) '_events.txt'], 'TRUNK');
            fprintf(['Eye data for Run ' num2str(irun) ' has been saved' '\n']);
        end
        
    end
    
    
    %% End of experiment: close and save Exp
    %     ListenChar(1);
    
    % close all on and off Screens
    Screen('CloseAll');
    
    %Stop recording & close Tobii
    if Exp.Cfg.trackEye
        fprintf('stop recording...');
        
        TobiiClose(Exp,eyeFileName,eventFileName)
        fprintf('done.\n');
    end
    
    % save parameters
    fprintf('saving matlab file...');
    save(fileName, 'Exp');
    fprintf(['All data has been saved' '\n']);
    
    % Show Mouse pointer
    ShowCursor;
    
    
catch ME1
    rethrow(ME1)
end

function [Exp, fileName, ELName, sessionNum] = setPaths(Exp, subjectID)


if  Exp.Cfg.computer.windows == 1
    
    % Add to the path the files to communicate with the eye tracker
    if ~exist(sprintf('..\\dataRaw\\%s\\',subjectID));
        mkdir(sprintf('..\\dataRaw\\%s\\',subjectID));
    end
    
    % filename
    sessionNum = 1;
    while exist(sprintf('..\\dataRaw\\%s\\conf1C_%s%3.3d.mat',subjectID,subjectID,sessionNum));
        sessionNum = sessionNum+1;
    end
    
    fileName = sprintf('..\\dataRaw\\%s\\conf1C_%s%3.3d.mat',subjectID,subjectID,sessionNum);
    ELName   = sprintf('%s%3.3d.edf',subjectID,sessionNum);
    
elseif Exp.Cfg.computer.linux == 1 || Exp.Cfg.computer.osx == 1
    
    % Add to the path the files to communicate with the eye tracker
    dataDir = sprintf('../dataRaw/%s%s/',Exp.Gral.subjectCategory,subjectID);
    if ~exist(dataDir);
        mkdir(dataDir);
    end
    
    % filename
    sessionNum = 1;
    while exist(sprintf('../dataRaw/%s%s/%s%s_%2.2d.mat',Exp.Gral.subjectCategory,subjectID,Exp.Gral.subjectCategory,subjectID,sessionNum));
        sessionNum = sessionNum+1;
    end
    
    fileName = sprintf([dataDir '%s%s_%2.2d.mat'],Exp.Gral.subjectCategory,subjectID,sessionNum);
    ELName   = sprintf('%s%3.3d.edf',subjectID,sessionNum);
end


function ctrl=IsKey(key)

global KEYBOARD;
[keyIsDown,secs,keyCode]=PsychHID('KbCheck', KEYBOARD);
if ~isnumeric(key)
    kc = KbName(key);
else
    kc = key;
end;
ctrl=keyCode(kc);
return


function [Exp KeyFlag] = run_BR_blk_replay(Exp, fileName, irun, iblk, msg)
% 27 May 2015, MF changed a parameter as following. 
% It was: Exp(iblk).BR.replayoverlapwidthDeg = 14 
% Now it's: Exp(iblk).BR.replayoverlapwidthDeg = 3
% Realised filter overlap width is actually wider than stimulus (!)

%%% using single keyboard %%%
% rKey = KbName('0)');
% lKey = KbName('1!');

%%% using keypads %%%
rKeypad           = Exp.Cfg.rKeypad;
lKeypad           = Exp.Cfg.lKeypad;
mainKeyboard      = Exp.Cfg.mainKeyboard;
nTrials           = Exp.Run(irun).block(iblk).nTrials;
stimFrames        = Exp.Run(irun).block(iblk).stimFrames;
interStimInterval = Exp.Run(irun).block(iblk).interStimInterval;

replayMedianFrac             = [0.8, 1.2];
Exp.Gral.BlockName           = num2str(iblk);
% Exp.BR.replayoverlapwidthDeg = 14; %- in degree
Exp.BR.replayOverlapWidthPix = Exp.BR.replayOverlapWidthDeg * Exp.Cfg.pixelsPerDegree;

%%% preallocate button response
Exp.Run(irun).block(iblk).rRespReplay = nan(nTrials, stimFrames + interStimInterval);
Exp.Run(irun).block(iblk).lRespReplay = nan(nTrials, stimFrames + interStimInterval);


z = meshgrid(0:Exp.BR.visibleSize(2)-1, 1);

%% load button responce ===================================================
for itri = 1: nTrials
    
    for ieye = 1:2
        switch Exp.BR.gratingDir(ieye)
            case 1
                resp(:,ieye) = Exp.Run(irun).block(iblk-1).rResp(itri, 1 : stimFrames);
            case -1
                resp(:,ieye) = Exp.Run(irun).block(iblk-1).lResp(itri, 1 : stimFrames);
        end
    end
    resp(:,3) = resp(:,1) + resp(:,2); % # pressed button
    
    resp( ( resp(:,1) + resp(:,3)) == 2, 4 ) = -1; % single left button press
    resp( ( resp(:,2) + resp(:,3)) == 2, 4 ) = 1; % single right button press
    resp( resp(:,3) == 0, 4 ) = 0; % no button press
    resp( resp(:,3) == 2, 4 ) = 2; % double button press
    
    [durInterpolateNo]     = getdur(resp(:,3), 0);
    [durInterpolateDouble] = getdur(resp(:,3), 2);
    durInterpolate         = cat(1, durInterpolateNo, durInterpolateDouble);
    
    if ~isempty(durInterpolate)
        
        durInterpolate(:,3)    = durInterpolate(:,2) - durInterpolate(:,1);
        
        %%% Interpolate missing button press ======================================
        for idur = 1 : size(durInterpolate,1)
            tStart = durInterpolate(idur,1);
            tEnd = durInterpolate(idur,2);
            
            
            if tStart > 1 && tEnd < Exp.Run(irun).block(iblk).stimFrames% if it's not the very beginning
                
                respBef = resp(tStart - 1, 4);
                respAft = resp(tEnd + 1, 4);
                
                changeButton = respBef - respAft;
                
                if changeButton == 0 % button presses before and after the no/double press duration are the same
                    resp(tStart : tEnd, 5) = respBef;
                    
                elseif changeButton == 2 % right -> no/double -> left
                    tHalfDur = floor((tStart + tEnd)/2);
                    resp(tStart : tHalfDur, 5) = 1;
                    resp(tHalfDur + 1 : tEnd, 5) = -1;
                    
                elseif changeButton == -2 % left -> no/double -> right
                    tHalfDur = ceil((tStart + tEnd)/2);
                    resp(tStart : tHalfDur, 5) = -1;
                    resp(tHalfDur + 1 : tEnd, 5) = 1;
                    
                end
                
            elseif tStart > 1 && tEnd == Exp.Run(irun).block(iblk).stimFrames % if it's the very end
                respBef = resp(tStart - 1, 4);
                resp(tStart:tEnd,5) = respBef;
                
            elseif tStart == 1
                %%% if it's the very beginning, fill the up the no/double duration
                %%% only with responce after the duration
                respAft = resp(tEnd + 1, 4);
                resp(1:tEnd,5) = respAft;
                
            elseif tStard == 0
                error('tStart is zero, somthing is wrong')
                
            end
            
        end
        resp(resp(:,5) == 0,5) = resp(resp(:,5) == 0,4); % single button press duration
        
    else %if there is no point to interpolate
        resp(:,5) = resp(:,4);
    end

    %- resp(:,1) : Left response key (1 or 0)
    %- resp(:,2) : Right response key (1 or 0)
    %- resp(:,3) : # of pressed key
    %- resp(:,4) : response type (0 = no press, 1 = right, -1 = left, 2 = double)
    %- resp(:,5) : interpolated response type (only 1 or -1)
    
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
    
    med = median(dur(:,4)); % median of dominance duration in frames
    
    Exp.Run(irun).block(iblk).replayMedianFrac{itri} = replayMedianFrac( round( rand(1, size(dur,1) ) ) + 1);

    
    % set start position of the border depending on the very first response
    switch resp(1,5)
        case 1
            Exp.Run(irun).block(iblk).replayBorderPos(itri,1) = Exp.BR.visibleSize(1);
        case -1
            Exp.Run(irun).block(iblk).replayBorderPos(itri,1) = 0;
    end
    
    respAllTri(itri,:,:) = resp;
    durAllTri{itri}      = dur;
    
    clear resp dur
end

    Exp.Run(irun).block(iblk).replayDesignFrame   = durAllTri;
    Exp.Run(irun).block(iblk).replayAltanatedResp = respAllTri;
    %-!! note that the stimulus is going to be presented up to nFrames, not the end of the duration.
    % means if the very end of the stimulus duration exeeds nFrames, it will be
    % shortened.
    

%% Trigger the beginning of the block
if Exp.Cfg.trackEye
    mesg       = [ msg.startBlk '_REPLAY'];
    blockStart = talk2tobii('EVENT', mesg, iblk);
    WaitSecs(0.2);
end



% z = 1 : Exp.BR.visibleSize;

for itri = 1 : nTrials
    idur = 1;
    
    %% PRESENT STIMULUS ===================================================
    timeTrial = nan(1, stimFrames);
    xoffset   = 0;
    KeyFlag   = 0;
    
    for iFrame = 1 : stimFrames
        
        % Clean screen
        Screen('FillRect',  Exp.Cfg.win, Exp.Cfg.WinColor);
        
        % Draw frames for binocular fusion
        Screen('DrawTextures', Exp.Cfg.win, Exp.BR.frametex, [0 0 size(Exp.BR.frame); 0 0 size(Exp.BR.frame)]', Exp.BR.framePos );
        
        
        % make moving sinusoidal curve for stimulus------------------------
        %         xoffset = (iFrame - 1) * Exp.BR.shiftPerFrame /(2*pi); % in pixel(pi)
        xoffset = mod( (iFrame - 1) * Exp.BR.shiftCyclePerFrameReplay * 2 * pi,  2 * pi);
        Exp.BR.replayGrating(1,:) = 0.5 * cos( Exp.BR.FsRadPerPix * z - Exp.BR.gratingDir(1) * xoffset ) + 0.5;
        Exp.BR.replayGrating(2,:) = 0.5 * cos( Exp.BR.FsRadPerPix * z - Exp.BR.gratingDir(2) * xoffset ) + 0.5;
        
        
        % Create filters --------------------------------------------------
        if iFrame >= durAllTri{itri}(idur,1)
            
            transdur = Exp.Run(irun).block(iblk).replayMedianFrac{itri}(idur) * (med/2); % frames/ swipe
            borderMovPixPerFrame = (Exp.BR.visibleSize(1) + Exp.BR.replayOverlapWidthPix ) / transdur ; % pixel/ frame
            
            if iFrame <= durAllTri{itri}(idur,1) + transdur % if the frame is in the transition duration
                borderDir = durAllTri{itri}(idur,3);
                
            elseif iFrame > durAllTri{itri}(idur,1) + transdur % if the frame has past the transition duration
                borderDir = 0; % then, don't move!
                
            end
            
            if iFrame == durAllTri{itri}(idur,2)
                idur = idur + 1; % for the next frame
            end
            
        end
        
        if iFrame ~= 1
            if Exp.Run(irun).block(iblk).replayBorderPos(itri,iFrame-1) <= ( 0 - Exp.BR.replayOverlapWidthPix/2 )...
                    && borderDir == -1     % if the border is on the left end and still wants to move to left
                Exp.Run(irun).block(iblk).replayBorderPos(itri,iFrame) = 0 - Exp.BR.replayOverlapWidthPix/2;
                
            elseif Exp.Run(irun).block(iblk).replayBorderPos(itri,iFrame-1) >= (Exp.BR.visibleSize(2) + Exp.BR.replayOverlapWidthPix/2) ...
                    && borderDir == 1    % if the border is on the right end and still wants to move to right
                Exp.Run(irun).block(iblk).replayBorderPos(itri,iFrame) = Exp.BR.visibleSize(2) + Exp.BR.replayOverlapWidthPix/2;
                
            else % if the border is somewhere middle of the stimulus and move more either to left/right
                Exp.Run(irun).block(iblk).replayBorderPos(itri,iFrame) = Exp.Run(irun).block(iblk).replayBorderPos(itri, iFrame - 1) + borderDir * borderMovPixPerFrame;
            end
        end
        
        a = Exp.Run(irun).block(iblk).replayBorderPos(itri,iFrame) - Exp.BR.replayOverlapWidthPix / 2; %<- It was originally 1/2.
        b = Exp.Run(irun).block(iblk).replayBorderPos(itri,iFrame) + Exp.BR.replayOverlapWidthPix / 2;
        
        for ipix = 1 : Exp.BR.visibleSize
            
            if z(ipix) <= a
                fr(ipix) = 1;
                
            elseif z(ipix) >= b
                fr(ipix) = 0;
                
            elseif a < z(ipix) && z(ipix) < b
                fr(ipix) = 0.5 * cos( ( z(ipix) - (a+b)/2 ) * ( pi/ Exp.BR.replayOverlapWidthPix ) + pi/2 ) + 0.5;
                %- filter for red (going the right to the left)
            end
            
        end
        
        fg = 1 - fr;
        % fr : filter for from-right grating
        % fl : filter for from-left grating
        
        
        %%% Apply the filter to the stimulus ------------------------------
        stim(:,:,1) = repmat( dot(fr, Exp.BR.replayGrating(find(Exp.BR.gratingColor(1,:)~=0),:),1) * Exp.Cfg.Color.red(1), Exp.BR.visibleSize(1), 1)  ; % red
        stim(:,:,2) = repmat( dot(fg, Exp.BR.replayGrating(find(Exp.BR.gratingColor(2,:)~=0),:),1) * Exp.Cfg.Color.green(2), Exp.BR.visibleSize(1), 1) ; % green
        stim(:,:,3) = zeros(Exp.BR.visibleSize(1), Exp.BR.visibleSize(1)); % blue
        
        % Draw filtered stimulus on the screen ----------------------------
        Exp.BR.replaytex = Screen('MakeTexture', Exp.Cfg.win, stim);
        Screen('DrawTextures', Exp.Cfg.win, Exp.BR.replaytex, [], Exp.BR.gratingDestRectPos);
        
        timeTrial(iFrame) = Screen('flip',Exp.Cfg.win);
        
        % Send trigger for the beginning of each trial
        if  iFrame == 1 && Exp.Cfg.trackEye
            mesg       = [ msg.startTri '_REPLAY'];
            trialStart = talk2tobii('EVENT', mesg, itri);
            WaitSecs(0.002);
        end
        
        % COLLECT CONTINUOUS RESPONSES HERE (one each frame)
        
        [keyIsDownR, ~, keyCodeR ] = KbCheck(rKeypad); % device #1= right hand, #2= left hand
        [keyIsDownL, ~, keyCodeL ] = KbCheck(lKeypad);
        [keyIsDown, ~, keyCode ]   = KbCheck(mainKeyboard);% Experimenter
        
        if keyIsDownR || keyIsDownL
            
            rPress = sum(keyCodeR);
            lPress = sum(keyCodeL);
            
            if rPress >= 1
                rResp = 1;
            elseif rPress == 0
                rResp = 0;
            end
            
            if lPress >= 1
                lResp = 1;
            elseif lPress == 0
                lResp = 0;
            end
            
        else
            rResp = 0;
            lResp = 0;
        end
        
        Exp.Run(irun).block(iblk).rResp(itri, iFrame ) = rResp;
        Exp.Run(irun).block(iblk).lResp(itri, iFrame ) = lResp;
        
        if keyIsDown
            % ESCAPE KEY
            % check_escapeKey ()
            if strcmpi(KbName(keyCode), Exp.addParams.escapeKey)
                
                % esc was pressed
                KeyFlag = 1;
                save([fileName(1:end-4) '_ESC' fileName(end-3:end)])
                
                if Exp.Cfg.trackEye
                    TobiiClose(Exp,eyeFileNameESC, eventFileNameESC)
                end
                
                ShowCursor;
                warning('Esc pressed, experiment ended.');
                break; % break from frame loop
            end
        end
        
    end
    
    if iFrame == stimFrames
        if Exp.Cfg.trackEye
            mesg     = [ msg.endTri '_REPLAY'];
            trialEnd = talk2tobii('EVENT', mesg, itri);
            WaitSecs(0.2)
        end
    end
    
    if KeyFlag, break, end % break from frame loop
    
    % INTER STIMULUS BLANK
    for iFrame = 1 : interStimInterval
        % Clean screen
        Screen('FillRect',  Exp.Cfg.win, Exp.Cfg.WinColor);
        % Draw frames for binocular fusion
        Screen('DrawTextures', Exp.Cfg.win, Exp.BR.frametex, [0 0 size(Exp.BR.frame); 0 0 size(Exp.BR.frame)]', Exp.BR.framePos );
        
        Screen('flip',Exp.Cfg.win);
        
        % COLLECT CONTINUOUS RESPONSES HERE (one each frame)
        
        [keyIsDownR, ~, keyCodeR ] = KbCheck(rKeypad); % device #1= right hand, #2= left hand
        [keyIsDownL, ~, keyCodeL ] = KbCheck(lKeypad);
        [keyIsDown, ~, keyCode ]   = KbCheck(mainKeyboard);% Experimenter
        
        if keyIsDownR || keyIsDownL
            
            rPress = sum(keyCodeR);
            lPress = sum(keyCodeL);
            
            if rPress >= 1
                rResp = 1;
            elseif rPress == 0
                rResp = 0;
            end
            
            if lPress >= 1
                lResp = 1;
            elseif lPress == 0
                lResp = 0;
            end
     
        else
            rResp = 0;
            lResp = 0;
        end
        Exp.Run(irun).block(iblk).rResp(itri, stimFrames + iFrame ) = rResp;
        Exp.Run(irun).block(iblk).lResp(itri, stimFrames + iFrame ) = lResp;
        
        if keyIsDown
            % ESCAPE KEY
            % check_escapeKey ()
            if strcmpi(KbName(keyCode), Exp.addParams.escapeKey)
                
                % esc was pressed
                KeyFlag = 1;
                save([fileName(1:end-4) '_ESC' fileName(end-3:end)])
                
                if Exp.Cfg.trackEye
                    TobiiClose(Exp,eyeFileNameESC, eventFileNameESC)
                end
                
                ShowCursor;
                warning('Esc pressed, experiment ended.');
                break; % break from frame loop
            end
            
        end
        
    end
    
    
    Exp.Run(irun).block(iblk).replayTimeTrial(itri,:) = timeTrial;
    
    % Careful here we need to consider the no ITI trials, we
    % still need to send the triggerExp.Cfg.win
    %                 if Exp.Cfg.trackEye
    %                     if Exp.Run(rn).block(blk).interStimInterval == 0
    %                         trialstart = iViewX('message', Exp.ivx, 'END_TR');
    %                         WaitSecs(0.2)
    %                     end
    %                 end
end

% if KeyFlag, break, end % break from block loop

% Trigger the end of the block
if Exp.Cfg.trackEye
    mesg     = [ msg.endBlk '_REPLAY'];
    blockEnd = talk2tobii('EVENT', mesg, iblk);
    WaitSecs(0.2);
end

% Show screen for inter trial inteval
% Clean screen
Screen('FillRect',  Exp.Cfg.win, Exp.Cfg.WinColor);
Screen('DrawTextures', Exp.Cfg.win, Exp.BR.frametex, [], Exp.BR.framePos);
Screen(Exp.Cfg.win, 'DrawText', 'End of Block,', Exp.Cfg.centerX - 475, Exp.Cfg.centerY-30, Exp.Cfg.Color.white);
Screen(Exp.Cfg.win, 'DrawText', 'Experimenter: press', Exp.Cfg.centerX - 475, Exp.Cfg.centerY, Exp.Cfg.Color.white);
Screen(Exp.Cfg.win, 'DrawText', '"space"to continue.', Exp.Cfg.centerX - 475, Exp.Cfg.centerY+30, Exp.Cfg.Color.white);
Screen(Exp.Cfg.win, 'DrawText', 'End of Block,', Exp.Cfg.centerX + 225, Exp.Cfg.centerY-30, Exp.Cfg.Color.white);
Screen(Exp.Cfg.win, 'DrawText', 'Experimenter: press', Exp.Cfg.centerX + 225, Exp.Cfg.centerY, Exp.Cfg.Color.white);
Screen(Exp.Cfg.win, 'DrawText', '"space"to continue.', Exp.Cfg.centerX + 225, Exp.Cfg.centerY+30, Exp.Cfg.Color.white);
Screen('flip',Exp.Cfg.win);

% Continue only after pressing space
while 1
    [keyIsDown, ~, keyCode ] = KbCheck(-3);
    if keyIsDown
        if ismember(KbName(keyCode), Exp.addParams.spaceKey)
            break;
        end
    end
    WaitSecs(0.2);
end
WaitSecs(0.5);


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

function [Prac KeyFlag] = prac_run_BR_blk_replay(Prac, fileName, irun, iblk, msg)

% 27 May 2015, MF changed a parameter as following. 
% It was: Prac.BR.replayoverlapwidthDeg = 14 
% Now it's: Prac.BR.replayoverlapwidthDeg = 5
% Realised filter overlap width is actually wider than stimulus (!)

%%% using single keyboard %%%
% rKey = KbName('0)');
% lKey = KbName('1!');

%%% using keypads %%%
rKeypad      = Prac.Cfg.rKeypad;
lKeypad      = Prac.Cfg.lKeypad;
mainKeyboard = Prac.Cfg.mainKeyboard;

nTrials    = Prac.Run(irun).block(iblk).nTrials;
stimFrames = Prac.Run(irun).block(iblk).stimFrames;
med = Prac.Run(irun).block(iblk).replayMedianDuration;

z = meshgrid(0:Prac.BR.visibleSize(2)-1, 1);

%% Trigger the beginning of the block
if Prac.Cfg.trackEye
    mesg       = [ msg.startBlk '_REPLAY'];
    blockStart = talk2tobii('EVENT', mesg, iblk);
    WaitSecs(0.2);
end

for itri = 1 : nTrials
    idur = 1;
    
    durAllTri = Prac.Run(irun).block(iblk).replayDesignFrame{itri};
    %% PRESENT STIMULUS ===================================================
    timeTrial = nan(1,stimFrames);
    xoffset   = 0;
    KeyFlag   = 0;

    for iFrame = 1 : stimFrames
        
        % Clean screen
        Screen('FillRect',  Prac.Cfg.win, Prac.Cfg.WinColor);
        
        % Draw frames for binocular fusion
        Screen('DrawTextures', Prac.Cfg.win, Prac.BR.frametex, [0 0 size(Prac.BR.frame); 0 0 size(Prac.BR.frame)]', Prac.BR.framePos );
        
        
        % make moving sinusoidal curve for stimulus------------------------
        %         xoffset = (iFrame - 1) * Prac.BR.shiftPerFrame /(2*pi); % in pixel(pi)
        xoffset = mod( (iFrame - 1) * Prac.BR.shiftCyclePerFrameReplay * 2 * pi,  2 * pi);
        Prac.Run(irun).block(iblk).replayGrating(1,:) = 0.5 * cos( Prac.BR.FsRadPerPix * z - Prac.BR.gratingDir(1) * xoffset ) + 0.5;
        Prac.Run(irun).block(iblk).replayGrating(2,:) = 0.5 * cos( Prac.BR.FsRadPerPix * z - Prac.BR.gratingDir(2) * xoffset ) + 0.5;
        
        
        % Create filters --------------------------------------------------
        if iFrame >= durAllTri(idur,1)
            
            transdur = Prac.Run(irun).block(iblk).replayMedianFrac{itri}(idur) * (med/2); % frames/ swipe
            borderMovPixPerFrame = (Prac.BR.visibleSize(1) + Prac.BR.replayOverlapWidthPix ) / transdur ; % pixel/ frame
            
            if iFrame <= durAllTri(idur,1) + transdur % if the frame is in the transition duration
                borderDir = durAllTri(idur,3);
                
            elseif iFrame > durAllTri(idur,1) + transdur % if the frame has past the transition duration
                borderDir = 0; % then, don't move!
            end
            
            if iFrame == durAllTri(idur,2)
                idur = idur + 1; % for the next frame
            end
        end
        
        if iFrame ~= 1
            if Prac.Run(irun).block(iblk).replayBorderPos(itri,iFrame-1) <= ( 0 - Prac.BR.replayOverlapWidthPix/2 )...
                    && borderDir == -1     % if the border is on the left end and still wants to move to left
                Prac.Run(irun).block(iblk).replayBorderPos(itri,iFrame) = 0 - Prac.BR.replayOverlapWidthPix/2;
                
            elseif Prac.Run(irun).block(iblk).replayBorderPos(itri,iFrame-1) >= (Prac.BR.visibleSize(2) + Prac.BR.replayOverlapWidthPix/2) ...
                    && borderDir == 1    % if the border is on the right end and still wants to move to right
                Prac.Run(irun).block(iblk).replayBorderPos(itri,iFrame) = Prac.BR.visibleSize(2) + Prac.BR.replayOverlapWidthPix/2;
                
            else % if the border is somewhere middle of the stimulus and move more either to left/right
                Prac.Run(irun).block(iblk).replayBorderPos(itri,iFrame) = Prac.Run(irun).block(iblk).replayBorderPos(itri, iFrame - 1) + borderDir * borderMovPixPerFrame;
            end
        end
        
         a = Prac.Run(irun).block(iblk).replayBorderPos(itri,iFrame) - Prac.BR.replayOverlapWidthPix / 2; %<- It was originally 1/2.
        b = Prac.Run(irun).block(iblk).replayBorderPos(itri,iFrame) + Prac.BR.replayOverlapWidthPix / 2;
        
        for ipix = 1 : Prac.BR.visibleSize
            
            if z(ipix) <= a
                fr(ipix) = 1;
                
            elseif z(ipix) >= b
                fr(ipix) = 0;
                
            elseif a < z(ipix) && z(ipix) < b
                fr(ipix) = 0.5 * cos( ( z(ipix) - (a+b)/2 ) * ( pi/ Prac.BR.replayOverlapWidthPix ) + pi/2 ) + 0.5;
                %- filter for red (going the right to the left)
            end
            
        end
        
        fg = 1 - fr;
        % fr : filter for from-right grating
        % fl : filter for from-left grating
        
        %%% Apply the filter to the stimulus ------------------------------
        stim(:,:,1) = repmat( dot(fr, Prac.Run(irun).block(iblk).replayGrating(find(Prac.BR.gratingColor(1,:)~=0),:),1) * Prac.Cfg.Color.red(1), Prac.BR.visibleSize(1), 1)  ; % red
        stim(:,:,2) = repmat( dot(fg, Prac.Run(irun).block(iblk).replayGrating(find(Prac.BR.gratingColor(2,:)~=0),:),1) * Prac.Cfg.Color.green(2), Prac.BR.visibleSize(1), 1) ; % green
        stim(:,:,3) = zeros(Prac.BR.visibleSize(1), Prac.BR.visibleSize(1)); % blue
        
        % Draw filtered stimulus on the screen ----------------------------
        Prac.BR.replaytex = Screen('MakeTexture', Prac.Cfg.win, stim);
        Screen('DrawTextures', Prac.Cfg.win, Prac.BR.replaytex, [], Prac.BR.gratingDestRectPos);
        
        timeTrial(iFrame) = Screen('flip',Prac.Cfg.win);
        
        % Send trigger for the beginning of each trial
        if  iFrame == 1 && Prac.Cfg.trackEye
            mesg       = [ msg.startTri '_REPLAY'];
            trialStart = talk2tobii('EVENT', mesg, itri);
            WaitSecs(0.002);
        end
        
        % COLLECT CONTINUOUS RESPONSES HERE (one each frame)
        [keyIsDownR, ~, keyCodeR ] = KbCheck(rKeypad); % device #1= right hand, #2= left hand
        [keyIsDownL, ~, keyCodeL ] = KbCheck(lKeypad);
        [keyIsDown, ~, keyCode ]   = KbCheck(mainKeyboard);% Prac.Run(irun).block(iblk)erimenter
        
        if keyIsDownR || keyIsDownL
            
            Rpress = sum(keyCodeR);
            Lpress = sum(keyCodeL);
            
            if Rpress >= 1
                rRespPress = 1;
            elseif Rpress == 0
                rRespPress = 0;
            end
            
            if Lpress >= 1
                lRespPress = 1;
            elseif Lpress == 0
                lRespPress = 0;
            end
            
            Prac.Run(irun).block(iblk).rResp(itri, iFrame ) = rRespPress;
            Prac.Run(irun).block(iblk).lResp(itri, iFrame ) = lRespPress;
            
        end
        
%         if keyIsDown
%             % ESCAPE KEY
%             % check_escapeKey ()
%             if strcmpi(KbName(keyCode), Prac.addParams.escapeKey)
%                 
%                 % esc was pressed
%                 KeyFlag = 1;
%                 save([fileName(1:end-4) '_ESC' fileName(end-3:end)])
%                 
%                 if Prac.Cfg.trackEye
%                     TobiiClose(Prac.Run(irun).block(iblk),[fileName(1:end-4) '_ESC' fileName(end-3:end)])
%                 end
%                 
%                 ShowCursor;
%                 warning('Esc pressed, Prac.Run(irun).block(iblk)eriment ended.');
%                 break; % break from frame loop
%             end
%         end
        
    end
    
    if KeyFlag, break, end % break from frame loop
    
    % INTER STIMULUS BLANK
    for iFrame = 1 : Prac.Run(irun).block(iblk).interStimInterval
        % Clean screen
        Screen('FillRect',  Prac.Cfg.win, Prac.Cfg.WinColor);
        % Draw frames for binocular fusion
        Screen('DrawTextures', Prac.Cfg.win, Prac.BR.frametex, [0 0 size(Prac.BR.frame); 0 0 size(Prac.BR.frame)]', Prac.BR.framePos );
        
        Screen('flip',Prac.Cfg.win);
        
        % COLLECT CONTINUOUS RESPONSES HERE (one each frame)
        [keyIsDownR, ~, keyCodeR ] = KbCheck(rKeypad); % device #1= right hand, #2= left hand
        [keyIsDownL, ~, keyCodeL ] = KbCheck(lKeypad);
        [keyIsDown, ~, keyCode ]   = KbCheck(mainKeyboard);% Prac.Run(irun).block(iblk)erimenter
        
        if keyIsDownR || keyIsDownL
            
            Rpress = sum(keyCodeR);
            Lpress = sum(keyCodeL);
            
            if Rpress >= 1
                rRespPress = 1;
            elseif Rpress == 0
                rRespPress = 0;
            end
            
            if Lpress >= 1
                lRespPress = 1;
            elseif Lpress == 0
                lRespPress = 0;
            end
            
            Prac.Run(irun).block(iblk).rResp(itri, stimFrames + iFrame ) = rRespPress;
            Prac.Run(irun).block(iblk).lResp(itri, stimFrames + iFrame ) = lRespPress;
        end
        
        if keyIsDown
            % ESCAPE KEY
            % check_escapeKey ()
            if strcmpi(KbName(keyCode), Prac.addParams.escapeKey)
                
                % esc was pressed
                KeyFlag = 1;
                save([fileName(1:end-4) '_ESC' fileName(end-3:end)])
                
                if Prac.Cfg.trackEye
                    TobiiClose(Prac.Run(irun).block(iblk),[fileName(1:end-4) '_ESC' fileName(end-3:end)])
                end
                
                ShowCursor;
                warning('Esc pressed, Prac.Run(irun).block(iblk)eriment ended.');
                break; % break from frame loop
            end
            
        end
        
        if iFrame == 1
            if Prac.Cfg.trackEye
                mesg     = [ msg.endTri '_REPLAY'];
                trialEnd = talk2tobii('EVENT', mesg, itri);
                WaitSecs(0.2)
            end
        end
        
    end
    
    Prac.Run(irun).block(iblk).timeTrial(itri,:) = timeTrial;
    
    % Careful here we need to consider the no ITI trials, we
    % still need to send the triggerPrac.Cfg.win
    %                 if Prac.Cfg.trackEye
    %                     if Prac.Run(irun).block(iblk).Run(rn).block(blk).interStimInterval == 0
    %                         trialstart = iViewX('message', Prac.Run(irun).block(iblk).ivx, 'END_TR');
    %                         WaitSecs(0.2)
    %                     end
    %                 end
end

% Trigger the end of the block
if Prac.Cfg.trackEye
    mesg     = [ msg.endBlk '_REPLAY'];
    blockEnd = talk2tobii('EVENT', mesg, iblk);
    WaitSecs(0.2);
end

% Show screen for inter trial inteval
% Clean screen
Screen('FillRect',  Prac.Cfg.win, Prac.Cfg.WinColor);
Screen('DrawTextures', Prac.Cfg.win, Prac.BR.frametex, [], Prac.BR.framePos);
Screen(Prac.Cfg.win, 'DrawText', 'End of Block,', Prac.Cfg.centerX - 475, Prac.Cfg.centerY-30, Prac.Cfg.Color.white);
Screen(Prac.Cfg.win, 'DrawText', 'Experimenter: press', Prac.Cfg.centerX - 475, Prac.Cfg.centerY, Prac.Cfg.Color.white);
Screen(Prac.Cfg.win, 'DrawText', '"space"to continue.', Prac.Cfg.centerX - 475, Prac.Cfg.centerY+30, Prac.Cfg.Color.white);
Screen(Prac.Cfg.win, 'DrawText', 'End of Block,', Prac.Cfg.centerX + 225, Prac.Cfg.centerY-30, Prac.Cfg.Color.white);
Screen(Prac.Cfg.win, 'DrawText', 'Experimenter: press', Prac.Cfg.centerX + 225, Prac.Cfg.centerY, Prac.Cfg.Color.white);
Screen(Prac.Cfg.win, 'DrawText', '"space"to continue.', Prac.Cfg.centerX + 225, Prac.Cfg.centerY+30, Prac.Cfg.Color.white);
Screen('flip',Prac.Cfg.win);

% Continue only after pressing space
while 1
    [keyIsDown, ~, keyCode ] = KbCheck(-3);
    if keyIsDown
        if ismember(KbName(keyCode), Prac.addParams.spaceKey)
            break;
        end
    end
    WaitSecs(0.2);
end
WaitSecs(0.5);


function [Prac KeyFlag] = prac_run_BR_blk_br(Prac, fileName, irun, iblk, msg)
%%% using single keyboard
% rKey = KbName('0)');
% lKey = KbName('1!'); <- No need to use these if you are using numeric
% keypads

%%% using keypads
rKeypad = Prac.Cfg.rKeypad;
lKeypad = Prac.Cfg.lKeypad;
mainKeyboard = Prac.Cfg.mainKeyboard;

nTrials    = Prac.Run(irun).block(iblk).nTrials;
stimFrames = Prac.Run(irun).block(iblk).stimFrames;

%%% preallocate button response
% Prac.Run(irun).block(iblk).rResp = nan(nTrials,...
%     stimFrames + Prac.Run(irun).block(iblk).interStimInterval);
%
% Prac.Run(irun).block(iblk).lResp = nan(nTrials,...
%     stimFrames + Prac.Run(irun).block(iblk).interStimInterval);
%
% Prac.BR.grating           = Prac.Cfg.Color.gray + Prac.Cfg.Color.inc * cos( Prac.BR.Fsr * x);
% Prac.BR.gratingtex        = Screen('MakeTexture', Prac.Cfg.win, Prac.BR.grating, [], 1);

% Trigger the beginning of the block
if Prac.Cfg.trackEye
    mesg       = msg.startBlk;
    blockStart = talk2tobii('EVENT', mesg, iblk);
    WaitSecs(0.2);
end

for itri = 1 : nTrials
    % PRESENT STIMULUS
    xoffset   = 0;
    KeyFlag   = 0;
    timeTrial = nan(1, stimFrames);
    
    for iFrame = 1 : stimFrames
        % Clean screen
        Screen('FillRect',  Prac.Cfg.win, Prac.Cfg.WinColor);
        
        % Draw frames for binocular fusion
        Screen('DrawTextures', Prac.Cfg.win, Prac.BR.frametex, [0 0 size(Prac.BR.frame); 0 0 size(Prac.BR.frame)]', Prac.BR.framePos );
        
        % Draw Gratings on screen
        % Shift the grating by "shiftperframe" pixels per frame:
        xoffset = mod( (iFrame - 1) * Prac.BR.shiftPixPerFrame, Prac.BR.gratingPixPerCycle);
        %         xoffset = xoffset + Prac.BR.shiftPerFrame;
        Prac.BR.gratingSrcRectPos = [xoffset, 0, xoffset + Prac.BR.visibleSize(2), Prac.BR.visibleSize(1)];
        
        Screen('DrawTextures', Prac.Cfg.win, Prac.BR.gratingTex,...
            [ Prac.BR.gratingSrcRectPos; Prac.BR.gratingSrcRectPos]',...
            Prac.BR.gratingDestRectPos,...
            Prac.BR.gratingAng', [], [],Prac.BR.gratingColor' );
        
        timeTrial(iFrame) = Screen('flip',Prac.Cfg.win);
        
        % Send trigger for the beginning of each trial
        if iFrame == 1 && Prac.Cfg.trackEye
            trial_id = num2str(itri);
            %                         trial_id(ismember(trial_id,' ')) = [];
            
            mesg       = msg.startTri;
            trialStart = talk2tobii('EVENT', mesg, itri);
            WaitSecs(0.002);
        end
        
        % COLLECT CONTINUOUS RESPONSES HERE (one each frame)
        [keyIsDownR, ~, keyCodeR ] = KbCheck(rKeypad); % device #1= right hand, #2= left hand
        [keyIsDownL, ~, keyCodeL ] = KbCheck(lKeypad);
        [keyIsDown, ~, keyCode ]   = KbCheck(mainKeyboard);% Experimenter
        
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
        
        % ESCAPE KEY
        % check_escapeKey ()
        if keyIsDown
            if strcmpi(KbName(keyCode), Prac.addParams.escapeKey) % this can be done only by experimenter
                
                % esc was pressed
                KeyFlag = 1;
                save([fileName(1:end-4) '_ESC' fileName(end-3:end)])
                
                if Prac.Cfg.trackEye
                    TobiiClose(Exp, [fileName(1:end-4) '_ESC' fileName(end-3:end)])
                end
                
                ShowCursor;
                warning('Esc pressed, experiment ended.');
                break; % break from frame loop
            end
        end
    end
    
    if KeyFlag, break, end % break from trial loop
    
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
        [keyIsDown, ~, keyCode ] = KbCheck(mainKeyboard);% Experimenter
        
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
            if strcmpi(KbName(keyCode), Prac.addParams.escapeKey) % this can be done only by experimenter
                
                % esc was pressed
                KeyFlag = 1;
                save([fileName(1:end-4) '_ESC' fileName(end-3:end)])
                
                if Prac.Cfg.trackEye
                    TobiiClose(Exp, [fileName(1:end-4) '_ESC' fileName(end-3:end)])
                end
                
                ShowCursor;
                warning('Esc pressed, experiment ended.');
                break; % break from frame loop
            end
        end
        
        if iFrame == 1
            if Prac.Cfg.trackEye
                mesg     = msg.endTri;
                trialEnd = talk2tobii('EVENT', mesg, itri);
                WaitSecs(0.2)
            end
        end
        
    end
    
    Prac.Run(irun).block(iblk).timeTrial(itri,:) = timeTrial;
    
    % Careful here we need to consider the no ITI trials, we
    % still need to send the triggerPrac.Cfg.win
    %                 if Prac.Cfg.trackEye
    %                     if Prac.Run(rn).block(blk).interStimInterval == 0
    %                         trialstart = iViewX('message', Prac.ivx, 'END_TR');
    %                         WaitSecs(0.2)
    %                     end
    %                 end
end

% if KeyFlag, break, end % break from block loop

% Trigger the end of the block
if Prac.Cfg.trackEye
    mesg     = msg.endBlk;
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
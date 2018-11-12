function [Exp KeyFlag] = run_BR_blk_br(Exp, fileName, irun, iblk)
%%% using single keyboard
% rKey = KbName('0)');
% lKey = KbName('1!'); <- No need to use these if you are using numeric
% keypads

%%% using keypads
keypadR = Exp.Cfg.keypadR;
keypadL = Exp.Cfg.keypadL;
mainKeyboard = Exp.Cfg.mainKeyboard;

nTrials    = Exp.Run(irun).block(iblk).nTrials;
stimFrames = Exp.Run(irun).block(iblk).stimFrames;
interStimInterval = Exp.Run(irun).block(iblk).interStimInterval;

%%% preallocate button response
% Exp.Run(irun).block(iblk).rResp = nan(nTrials, stimFrames + interStimInterval);
% Exp.Run(irun).block(iblk).lResp = nan(nTrials, stimFrames + interStimInterval);
% Exp.BR.grating           = Exp.Cfg.Color.gray + Exp.Cfg.Color.inc * cos( Exp.BR.Fsr * x);
% Exp.BR.gratingtex        = Screen('MakeTexture', Exp.Cfg.win, Exp.BR.grating, [], 1);

% Trigger the beginning of the block
if Exp.Cfg.trackEye
    mesg       = 'START_BK';
    blockStart = talk2tobii('EVENT', mesg, iblk);
    WaitSecs(0.2);
end

for itri = 1 : nTrials
    % PRESENT STIMULUS
    xoffset   = 0;
    KeyFlag   = 0;
    TimeTrial = nan(1, stimFrames);
    
    for iFrame = 1 : stimFrames
        % Clean screen
        Screen('FillRect',  Exp.Cfg.win, Exp.Cfg.WinColor);
        
        % Draw frames for binocular fusion
        Screen('DrawTextures', Exp.Cfg.win, Exp.BR.frametex, [0 0 size(Exp.BR.frame); 0 0 size(Exp.BR.frame)]', Exp.BR.framePos );
        
        % Draw Gratings on screen
        % Shift the grating by "shiftperframe" pixels per frame:
        xoffset = mod( (iFrame - 1) * Exp.BR.shiftPixPerFrame, Exp.BR.gratingPixPerCycle);
        %         xoffset = xoffset + Exp.BR.shiftPerFrame;
        Exp.BR.gratingSrcRectPos = [xoffset, 0, xoffset + Exp.BR.visibleSize(2), Exp.BR.visibleSize(1)];
        
        Screen('DrawTextures', Exp.Cfg.win, Exp.BR.gratingTex,...
            [ Exp.BR.gratingSrcRectPos; Exp.BR.gratingSrcRectPos]',...
            Exp.BR.gratingDestRectPos,...
            Exp.BR.gratingAng', [], [],Exp.BR.gratingColor' );
        
        %             [ Exp.BR.gratingColor ./ repmat(max(Exp.BR.gratingColor,[],2),[1,3]) ]' );
        
        %           Screen('DrawTextures', Exp.Cfg.win, Exp.BR.gratingtex, ...
        %                         [Exp.BR.gratingSrcRect ; Exp.BR.gratingSrcRect]'  ,...
        %                         Exp.BR.gratingDestRectPos, Exp.BR.gratingAng', [], [],  Exp.BR.gratingColor );
        %
        %         xoffset = xoffset - Exp.BR.shiftperframe;
        
        %         scrRect = [xoffset, 0, xoffset + Exp.BR.visibleSize(1), Exp.BR.visibleSize(2)];
        %
        %         Screen('DrawTexture', window, gratingMaskTex, srcRect, dstRect, []);
        
        
        
        %         normalGrating(1,:) = 0.5 * cos( Exp.BR.Fsr * z - Exp.BR.gratingDir(1) * xoffset) + 0.5; % for left eye
        %         normalGrating(2,:) = 0.5 * cos( Exp.BR.Fsr * z - Exp.BR.gratingDir(2) * xoffset) + 0.5; % for right eye
        %
        %
        %         vectorL = repmat(normalGrating(1,:), [Exp.BR.visibleSize(1), 1, 3]);
        %         colorL = permute(repmat(Exp.BR.gratingColor(:,1),[1, Exp.BR.visibleSize(1), Exp.BR.visibleSize(2)]) , [2,3,1]);
        %         stimL = vectorL .* colorL;
        %
        %         vectorR = repmat(normalGrating(2,:), [Exp.BR.visibleSize(1), 1, 3]);
        %         colorR = permute(repmat(Exp.BR.gratingColor(:,2),[1, Exp.BR.visibleSize(1), Exp.BR.visibleSize(2)]) , [2,3,1]);
        %         stimR = vectorR .* colorR;
        %
        %
        %         Exp.BR.normaltexL = Screen('MakeTexture', Exp.Cfg.win, stimL);
        %         Screen('DrawTextures', Exp.Cfg.win, Exp.BR.normaltexL, [], Exp.BR.gratingDestRectPos(:,1));
        %
        %         Exp.BR.normaltexR = Screen('MakeTexture', Exp.Cfg.win, stimR);
        %         Screen('DrawTextures', Exp.Cfg.win, Exp.BR.normaltexR, [], Exp.BR.gratingDestRectPos(:,2));
        %
        
        TimeTrial(iFrame) = Screen('flip',Exp.Cfg.win);
        
        % Send trigger for the beginning of each trial
        if iFrame == 1 && Exp.Cfg.trackEye
            trial_id = num2str(itri);
            %                         trial_id(ismember(trial_id,' ')) = [];
            
            mesg       = 'START_TR';
            trialStart = talk2tobii('EVENT', mesg, itri);
            WaitSecs(0.002);
        end
        
        % COLLECT CONTINUOUS RESPONSES HERE (one each frame)
        [keyIsDownR, ~, keyCodeR ] = KbCheck(keypadR); % device #1= right hand, #2= left hand
        [keyIsDownL, ~, keyCodeL ] = KbCheck(keypadL);
        [keyIsDown, ~, keyCode ]   = KbCheck(mainKeyboard);% Experimenter
        
        if keyIsDownR || keyIsDownL
            
            Rpress = sum(keyCodeR);
            Lpress = sum(keyCodeL);
            
            if Rpress >= 1
                rResp = 1;
            elseif Rpress == 0
                rResp = 0;
            end
            
            if Lpress >= 1
                lResp = 1;
            elseif Lpress == 0
                lResp = 0;
            end
            
            Exp.Run(irun).block(iblk).rResp(itri, iFrame ) = rResp;
            Exp.Run(irun).block(iblk).lResp(itri, iFrame ) = lResp;
        end
        
        % ESCAPE KEY
        % check_escapeKey ()
        if keyIsDown
            if strcmpi(KbName(keyCode), Exp.addParams.escapeKey) % this can be done only by experimenter
                
                % esc was pressed
                KeyFlag = 1;
                save([fileName(1:end-4) '_ESC' fileName(end-3:end)])
                
                if Exp.Cfg.trackEye
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
    for iFrame = 1 : Exp.Run(irun).block(iblk).interStimInterval
        % Clean screen
        Screen('FillRect',  Exp.Cfg.win, Exp.Cfg.WinColor);
        % Draw frames for binocular fusion
        Screen('DrawTextures', Exp.Cfg.win, Exp.BR.frametex, [0 0 size(Exp.BR.frame); 0 0 size(Exp.BR.frame)]', Exp.BR.framePos );
        
        Screen('flip',Exp.Cfg.win);
        
        % COLLECT CONTINUOUS RESPONSES HERE (one each frame)
        [keyIsDownR, ~, keyCodeR ] = KbCheck(keypadR); % device #1= right hand, #2= left hand
        [keyIsDownL, ~, keyCodeL ] = KbCheck(keypadL);
        [keyIsDown, ~, keyCode ] = KbCheck(mainKeyboard);% Experimenter
        
        if keyIsDownR || keyIsDownL
            
            Rpress = sum(keyCodeR);
            Lpress = sum(keyCodeL);
            
            if Rpress >= 1
                rResp = 1;
            elseif Rpress == 0
                rResp = 0;
            end
            
            if Lpress >= 1
                lResp = 1;
            elseif Lpress == 0
                lResp = 0;
            end
            
            Exp.Run(irun).block(iblk).rResp(itri, stimFrames + iFrame ) = rResp;
            Exp.Run(irun).block(iblk).lResp(itri, stimFrames + iFrame ) = lResp;
            
        end
        
        if keyIsDown
            if strcmpi(KbName(keyCode), Exp.addParams.escapeKey) % this can be done only by experimenter
                
                % esc was pressed
                KeyFlag = 1;
                save([fileName(1:end-4) '_ESC' fileName(end-3:end)])
                
                if Exp.Cfg.trackEye
                    TobiiClose(Exp, [fileName(1:end-4) '_ESC' fileName(end-3:end)])
                end
                
                ShowCursor;
                warning('Esc pressed, experiment ended.');
                break; % break from frame loop
            end
        end
        
        if iFrame == 1
            if Exp.Cfg.trackEye
                mesg     = 'END_TR';
                trialEnd = talk2tobii('EVENT', mesg, itri);
                WaitSecs(0.2)
            end
        end
        
    end
    
    Exp.Run(irun).block(iblk).TimeTrial(itri,:) = TimeTrial;
    
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
    mesg     = 'END_BK';
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
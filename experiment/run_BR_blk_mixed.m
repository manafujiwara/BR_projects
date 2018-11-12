function [Exp KeyFlag] =  run_BR_blk_mixed(Exp, fileName, irun, iblk, msg)
                            
%%% using single keyboard
% rKey = KbName('0)');
% lKey = KbName('1!'); <- No need to use these if you are using numeric
% keypads

%%% using keypads
rKeypad      = Exp.Cfg.rKeypad;
lKeypad      = Exp.Cfg.lKeypad;
mainKeyboard = Exp.Cfg.mainKeyboard;

nTrials           = Exp.Run(irun).block(iblk).nTrials;
stimFrames        = Exp.Run(irun).block(iblk).stimFrames;
interStimInterval = Exp.Run(irun).block(iblk).interStimInterval;
nMixedTrials      = Exp.Run(irun).block(iblk).mixedRate * nTrials;
if sum(nMixedTrials ~= round(nMixedTrials))
    idecimal = find(nMixedTrials ~= round(nMixedTrials));
    i = randperm(sum(nMixedTrials ~= round(nMixedTrials)));
    nMixedTrials(idecimal(i(1))) = ceil(nMixedTrials(idecimal(i(1)))); %- ASSUMING THE NUMBER OF LABEL IS 3 FOR NOW
    nMixedTrials(idecimal(i(2))) = floor(nMixedTrials(idecimal(i(2))));
end

Exp.Run(irun).block(iblk).rResp = zeros(nTrials, stimFrames + interStimInterval);
Exp.Run(irun).block(iblk).lResp = zeros(nTrials, stimFrames + interStimInterval);

labels = 1 : length(Exp.Run(irun).block(iblk).mixedRate);

tmp = [];
for ilabels = 1 : length(labels)
    tmp = [tmp labels(ilabels)*ones(1,nMixedTrials(ilabels))];
end

mixedLabelDesign = tmp(randperm(Exp.Run(irun).block(iblk).nTrials));
Exp.Run(irun).block(iblk).mixedLabelDesign = mixedLabelDesign; % 1)br, 2)red/right unambiguous, 3)green/left unambiduous

%%% preallocate button response
% Exp.Run(irun).block(iblk).rResp = nan(nTrials, stimFrames + interStimInterval);
% Exp.Run(irun).block(iblk).lResp = nan(nTrials, stimFrames + interStimInterval);
% Exp.BR.grating    = Exp.Cfg.Color.gray + Exp.Cfg.Color.inc * cos( Exp.BR.Fsr * x);
% Exp.BR.gratingtex = Screen('MakeTexture', Exp.Cfg.win, Exp.BR.grating, [], 1);

% Trigger the beginning of the block
if Exp.Cfg.trackEye
    mesg       = msg.startBlk;
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
        gratingSrcRectPos = [xoffset, 0, xoffset + Exp.BR.visibleSize(2), Exp.BR.visibleSize(1)];
        gratingDestRectPos = Exp.BR.gratingDestRectPos;
        
        switch  mixedLabelDesign(itri)
            case 1 % br
                gratingAng         = Exp.BR.gratingAng' ;
                gratingColor       = Exp.BR.gratingColor' ;
                
            case 2 % unambiguous red/right
                gratingAng         = [Exp.BR.gratingAng(2) Exp.BR.gratingAng(2)]';
                gratingColor       = [Exp.Cfg.Color.red; Exp.Cfg.Color.red]' ;
                
            case 3 % unambiguous green/left
                gratingAng         = [Exp.BR.gratingAng(1) Exp.BR.gratingAng(1)]';
                gratingColor       = [Exp.Cfg.Color.green; Exp.Cfg.Color.green]' ;
        end
        
        Screen('DrawTextures', Exp.Cfg.win, Exp.BR.gratingTex,...
            [ gratingSrcRectPos; gratingSrcRectPos]',gratingDestRectPos,...
            gratingAng, [], [], gratingColor);
         
        TimeTrial(iFrame) = Screen('flip',Exp.Cfg.win);
        
        % Send trigger for the beginning of each trial
        if iFrame == 1 && Exp.Cfg.trackEye
            trial_id   = num2str(itri);
            mesg       = msg.startTri;
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
            
            Exp.Run(irun).block(iblk).rResp(itri, iFrame ) = rResp;
            Exp.Run(irun).block(iblk).lResp(itri, iFrame ) = lResp;
        end

        
        % ESCAPE KEY
        if keyIsDown
            if strcmpi(KbName(keyCode), Exp.addParams.escapeKey) % this can be done only by experimenter
                
                % esc was pressed
                KeyFlag = 1;
                save([fileName(1:end-4) '_ESC' fileName(end-3:end)])
                
                if Exp.Cfg.trackEye
                    TobiiClose(Exp, eyeFileNameESC, eventFileNameESC)
                end
                
                ShowCursor;
                warning('Esc pressed, experiment ended.');
                break; % break from frame loop
            end
        end
    end
    
    if KeyFlag, break, end % break from trial loop
    
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
            Exp.Run(irun).block(iblk).rResp(itri, stimFrames + iFrame ) = rResp;
            Exp.Run(irun).block(iblk).lResp(itri, stimFrames + iFrame ) = lResp;
            
        end
        
        
        if keyIsDown
            if strcmpi(KbName(keyCode), Exp.addParams.escapeKey) % this can be done only by experimenter
                
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
        
        if iFrame == 1
            if Exp.Cfg.trackEye
                mesg     = msg.endTri;
                trialEnd = talk2tobii('EVENT', mesg, itri);
                WaitSecs(0.2)
            end
        end
        
    end
    
    Exp.Run(irun).block(iblk).TimeTrial(itri,:) = TimeTrial;
end

% Trigger the end of the block
if Exp.Cfg.trackEye
    mesg     = msg.endBlk;
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
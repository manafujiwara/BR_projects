%% Continue Only If Space Key Is Pressed

while 1
    [keyIsDown, ~, keyCode ] = KbCheck(-3);%- KbCheck will check all the keypads and keyboards connected
    if keyIsDown
        if ismember(KbName(keyCode), Exp.addParams.spaceKey)
            break; 
        elseif strcmpi(KbName(keyCode), Exp.addParams.escapeKey) % this can be done only by experimenter
            
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
    
    WaitSecs(0.2);
end
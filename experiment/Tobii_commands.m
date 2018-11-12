

trial_dur = 60;
samples_total = trial_dur * 300;

% First preallocate space for the samples

sample_l = 1;

samples_x = nan(1, samples_total); % time of the trial, in samples
samples_y = nan(1, samples_total);
gazeData = nan(samples_total, 17);
samp = 1;



while time <= trial_dur
    
    
    % RECORD EYE POSITION
    if Exp.addParams.witheye
        
        % check for the presence of a new sample update
        gazeData(samp, :)= talk2tobii ('GET_SAMPLE_EXT');
        
        % Do we have valid data FOR THE LEFT EYE and is the pupil visible?
        if gazeData(samp,7) <= 1 && gazeData(samp,8) <= 1
            
            % Get current gaze position from sample
            xl = gazeData(samp, 1) * Exp.Cfg.width;
            yl = gazeData(samp, 2) * Exp.Cfg.height;
            
            xr = gazeData(samp, 3) * Exp.Cfg.width;
            yr = gazeData(samp, 4) * Exp.Cfg.height;
            
            % Here I am averaging the positions of the right and left eye,
            % but we could use only one eye for our experiment
            samples_x(sample_l) = (xl + xr) / 2;
            samples_y(sample_l) = (yl + yr) / 2;
        end
    end
    
    
       sample_l = sample_l + 1;
       samp = samp + 1;
       
       WaitSecs(0.002);
    
end


%% this is how you send Tobii a message
mesg = ['endImage_' num2str(tr)];
timeStart = talk2tobii('EVENT', mesg, 2);

%%  Code to draw online the position of the eye on the screen


if Exp.addParams.drawEye == 1
    
    pupil_size = 30; %0.05 * Exp.Cfg.width;
    color_leftEye = [0 0 255];
    LwC = samples_x(sample_l);
    LhC = samples_y(sample_l);
    rectOvalL(1) = LwC - pupil_size/2;
    rectOvalL(2) = LhC - pupil_size/2;
    rectOvalL(3) = LwC + pupil_size/2;
    rectOvalL(4) = LhC + pupil_size/2;
    Screen('FillOval',Exp.Cfg.win, color_leftEye, rectOvalL);
    %flip
    Screen('Flip', Exp.Cfg.win, [], Exp.Cfg.AuxBuffers);

end









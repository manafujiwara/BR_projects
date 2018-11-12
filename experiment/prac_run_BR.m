function prac_run_BR()

%% Written by MF, 19 May 2015
% This is a practice version for run_BR.m

addpath('./aux_files/');

%%% For mac in oculomotor room %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% addpath(genpath('/../../../Psychophysics/'));
% addpath(genpath('/Applications/Psychtoolbox'));
% addpath(genpath('/Users/lisandrk/no-ppc/lib/matlab')); % to the latest version of t2t

%%% For local machine %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath(genpath('/Applications/Psychtoolbox'));
% PsychJavaTrouble;



%% check keyboards
[KeyboardIndices productNames, allInfos] = GetKeyboardIndices();
Exp.Cfg.keyboardInfo = allInfos;
goLeft = 0;

if length(KeyboardIndices) < 3
    error('Keyboards don''t be connected propery')
    
elseif length(KeyboardIndices) > 3
    error('More than 3 keyboards are connected')
    
else
    
    for iKeyBoard = 1:length(KeyboardIndices)
        switch productNames{iKeyBoard}
            case 'Apple Internal Keyboard / Trackpad'
                Exp.Cfg.MainKeyboard = KeyboardIndices(iKeyBoard);
                
            case 'USB Compliant Keypad'
                switch goLeft
                    case 0
                        Exp.Cfg.KeypadR = KeyboardIndices(iKeyBoard);
                        goLeft = 1;
                    case 1
                        Exp.Cfg.KeypadL = KeyboardIndices(iKeyBoard);
                end
        end
    end
    
end


%% Initialize screen
% PsychJavaTrouble; % Check there are no problems with Java
Exp.Cfg.SkipSyncTest = 0; %This should be '0' on a properly working NVIDIA video card. '1' skips the SyncTest.
Exp.Cfg.AuxBuffers   = 1; % '0' if no Auxiliary Buffers available, otherwise put it into '1'.

% Check for OpenGL compatibility
AssertOpenGL;
Screen('Preference','SkipSyncTests', Exp.Cfg.SkipSyncTest);


%% PARAMETER FOR EXPERIMENTER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Exp.Prac.screen = 1; % do you want to run practice screening session?
Exp.Prac.br     = 1; % do you want to run practice BR session?
Exp.prac.replay = 1; % do you want to run practice replay sessiom? (can say yes only if BR was run)

Exp.Cfg.WinSize    = [0 0 1000 500];  %Empty means whole screen
Exp.Cfg.xDimCm     = 51; %Length in cm of the screen in X, Tobii
Exp.Cfg.yDimCm     = 29; %Length in cm of the screen in Y, Tobii
Exp.Cfg.distanceCm = 74; %Viewing distance

Exp.Gral.SubjectName    = input('Please enter subject name:\n', 's');
Exp.Gral.SubjectNumber  = input('Please enter subject number:\n');
Exp = initializeScreen(Exp);
Exp.Cfg.WinColor = Exp.Cfg.Color.gray; % empty for the middle gray of the screen.

% global constants
Exp.SFM.drawEye    = 0; % Draw eye online on the screen
Exp.SFM.trackEye   = 0; % 1 : collect eye tracking data

Exp.BR.gratingDir   = [ -1, 1 ]; % 1)left eye, 2)right eye. -1 = towards left, 1= towards right
Exp.BR.gratingColor = [ Exp.Cfg.Color.green; Exp.Cfg.Color.red];
Exp.BR.replay       = 1; % run replay for each block(=1) or not(=0) % haven't had an option to run just a few blocks

runs                = 1;
nBlocks             = 2; % number of blocks for the whole experiment
blk_id              = {'30-0', '2-1' };
nTrials             = [1, 10];% 30 20]; % trials for each block
stimFrames          = [1800, 120];% 60 120]; % duration of stimulation -in frames-
interStimInterval   = [0, 60]; % duration of intervals between stimulation -in frames-
eye_samples         = ( stimFrames + interStimInterval ) * 5; % number of eye samples to collect during stimulation + blank in each trial

%------example---------------------------------------------------------
%     runs = 9;
%     nBlocks = 5; % number of blocks for the whole experiment
%     nTrials             = [40 30  20  12  17  15  12   9     1]; % trials for each block
%     stimFrames          = [30 60 120 240  30  60 120 240  3600]; % duration of stimulation -in frames-
%     interStimInterval   = [60 60  60  60 180 180 180 180     0]; % duration of intervals between stimulation -in frames-
%     run_id              ={'0.5-1' '1-1' '2-1' '4-1' '0.5-3' '1-3' '2-3' '4-3' '60-0'};
% ---------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%






resp = 1;

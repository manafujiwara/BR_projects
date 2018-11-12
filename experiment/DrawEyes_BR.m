function DrawEyes_BR(win, left_w_camera, left_h_camera, right_w_camera, right_h_camera, left_validity, right_validity , Exp)
%DrawEyes(win, w_camera, h_camera, left_validity, right_validity)
%
% w_camera: Left eye position seen by the camera.
%           0.0 is leftmost and 1.0 is rightmost.
% h_camera: Left eye position seen by the camera.
%           0.0 is topmost and 1.0 is bottommost.
% validity: How likely is it that the left eye was found? 
%           0 - Certainly (>99%),
%           1 - Probably (80%),
%           2 - (50%),
%           3 - Likely not (20%),
%           4 - Certainly not (0%) 

LCamW = left_w_camera;
RCamW = right_w_camera;
LCamH = left_h_camera;
RCamH = right_h_camera;

%estimate screen dimension
% screenNumber=Screen('WindowScreenNumber', win);
% res=Screen('Rect', screenNumber);
% res = res(3:4);

% %use one square display proportional to the screen resolution
% width = 0.2*res(1);
% height =  width;% 0.2*res(2); %
% 
% Exp.BR.rectPosFromCenter
% 
% %draw rectangles
% rectWin = Exp.BR.gratingDestRectPos;
% % rectWin(1,1) = Exp.BR.rectPosFromCenter - ;
% % rectWin(2,1) = res(2)/2 - Exp.Cfg.centerY;
% % rectWin(3,1) = res(1)/2 + Exp.BR.rectPosFromCenter;
% % rectWin(4,1) = res(2)/2 + Exp.Cfg.centerY;
% % 
% % rectWin(1,2) = res(1)/2 - Exp.BR.rectPosFromCenter;
% % rectWin(2,2) = res(2)/2 - Exp.Cfg.centerY;
% % rectWin(3,2) = res(1)/2 + Exp.BR.rectPosFromCenter;
% % rectWin(4,2) = res(2)/2 + Exp.Cfg.centerY;
% color_frame = [0 0 0 ];
% 
% Screen( 'FillRect',win,color_frame,rectWin );


Screen('FillRect',  Exp.Cfg.win, Exp.Cfg.WinColor);
Screen('DrawTextures', Exp.Cfg.win, Exp.BR.frametex, [], Exp.BR.framePos);


%draw eyes
%decide color
switch left_validity
    case 0,
        color_leftEye = [0 255 0];
    case 1, 
        color_leftEye = [64 192 0];
    case 2,
        color_leftEye = [128 128 0];
    case 3,
        color_leftEye = [192 64 0];
    otherwise,
        color_leftEye = [255 0 0];
end


switch right_validity
    case 0,
        color_rightEye = [0 255 0];
    case 1, 
        color_rightEye = [64 192 0];
    case 2,
        color_rightEye = [128 128 0];
    case 3,
        color_rightEye = [192 64 0];
    otherwise,
        color_rightEye = [255 0 0];
end


pupil_size = 0.1*Exp.BR.visibleSize(1);

if( left_validity<4)
    LwC(1) = Exp.BR.gratingDestRectPos(1,1) + Exp.BR.visibleSize(1)*LCamW;
    LhC(1) = Exp.BR.gratingDestRectPos(2,1) + Exp.BR.visibleSize(2)*LCamH;
    LwC(2) = Exp.BR.gratingDestRectPos(1,2) + Exp.BR.visibleSize(1)*LCamW;
    LhC(2) = Exp.BR.gratingDestRectPos(2,2) + Exp.BR.visibleSize(2)*LCamH;
    rectOvalL(1,:) = LwC - pupil_size/2;
    rectOvalL(2,:) = LhC - pupil_size/2;
    rectOvalL(3,:) = LwC + pupil_size/2;
    rectOvalL(4,:) = LhC + pupil_size/2;
    Screen('FillOval',win, color_leftEye,rectOvalL);
end

if(right_validity<4)
    RwC(1) = Exp.BR.gratingDestRectPos(1,1) + Exp.BR.visibleSize(1)*RCamW;
    RhC(1) = Exp.BR.gratingDestRectPos(2,1) + Exp.BR.visibleSize(2)*RCamH;
    RwC(2) = Exp.BR.gratingDestRectPos(1,2) + Exp.BR.visibleSize(1)*RCamW;
    RhC(2) = Exp.BR.gratingDestRectPos(2,2) + Exp.BR.visibleSize(2)*RCamH;
    rectOvalR(1,:) = RwC - pupil_size/2;
    rectOvalR(2,:) = RhC - pupil_size/2;
    rectOvalR(3,:) = RwC + pupil_size/2;
    rectOvalR(4,:) = RhC + pupil_size/2;
    Screen('FillOval',win, color_rightEye,rectOvalR);
end


%flip
Screen('Flip', win );

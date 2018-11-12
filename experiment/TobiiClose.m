function TobiiClose (Exp, eyeFileName, eventFileName)

talk2tobii('STOP_AUTO_SYNC');
talk2tobii('STOP_RECORD');
talk2tobii('STOP_TRACKING');

talk2tobii('SAVE_DATA', eyeFileName,eventFileName, 'TRUNK');

talk2tobii ('CLEAR_DATA');
talk2tobii('DISCONNECT');



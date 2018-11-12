function exportData2excel(DIR,subID,Cond)

downsamp = 1; % yes->1 (60 Hz), no->0 (300 Hz)

whicheye = Cond.whicheye;
time2cut = Cond.time2cut;
epochWin = Cond.epochWinIM;

time = linspace(0,60,18000);

for isub = 8: length(subID)
    
    loadFileName = [ subID{isub} '_intp' num2str(time2cut*10^3)  '_vel_slwPhs_epIM.mat' ];
    load( [ DIR.epoched '/' loadFileName] )
    load( [ DIR.cfg '/' subID{isub} '_Cfg_a_br.mat'] )
    
    CT_mat = [];
    IM_mat = [];
    subNum = str2num( subID{isub}(2:4) );
    %% compute velocity with boxcar smoothing =============================
    for irun = 1 : length(Run)
        for iblk = 1: length(Run(irun).block)
            
            switch Run(irun).block(iblk).state
                case 'nc'
                    state= 1; %=> block level
                case 'mon'
                    state = 2;
                case 'mof'
                    state = 3;
                case 'donmon'
                    state = 4;
                case 'dofmon'
                    state = 5;
                case 'donmof'
                    state = 4;
                case 'dofmof'
                    state = 5;
            end
            
            if isempty(Run(irun).block(iblk).gazeData)
                iOKN           = nan(18000,2);
                vel_smtiOKN    = nan(18000,2);
                smtvel_smtiOKN = nan(18000,2);
            else
                iOKN           = Run(irun).block(iblk).iOKN(:,[3,5]);
                vel_smtiOKN    = Run(irun).block(iblk).iOKN_vel(:,[3,5]);
                smtvel_smtiOKN = Run(irun).block(iblk).iOKN_slwPhsDegPerSec_box(:,[3,5]);
            end
            
            
            switch Run(irun).block(iblk).condition
                case 'br'
                    cond         = ones(1,18000) * 1;
                    lResp        = my_interp( Run(irun).block(iblk).lResp, 5 );
                    rResp        = my_interp( Run(irun).block(iblk).rResp, 5 );
                    borderPos    = nan(18000,1);
                    borderPos_lb = nan(18000,1);
                    CT_mat = cat(1, CT_mat, [repmat(subNum,18000,1), repmat(state,18000,1), time', cond',lResp, rResp, iOKN, ...
                        vel_smtiOKN, smtvel_smtiOKN, borderPos, borderPos_lb ] );
                    
                    
                case 'replay'
                    cond = ones(1,18000) * 4;
                    % Run(irun).block(iblk).borderPos;
                    lResp = my_interp( Run(irun).block(iblk).lResp, 5 );
                    rResp = my_interp( Run(irun).block(iblk).rResp, 5 );
                    
                    if isempty(Run(irun).block(iblk).gazeData)
                        borderPos = nan(1,3600);
                        borderPos_lb = nan(1,3600);
                    else
                        borderPos    = Run(irun).block(iblk).borderPos;
                        borderPos_lb(borderPos<135) = -1;
                        borderPos_lb(borderPos>=135) = 1;
                    end
                    
                    CT_mat = cat(1, CT_mat, [repmat(subNum,18000,1), repmat(state,18000,1), time', cond', lResp, rResp, iOKN, ...
                        vel_smtiOKN, smtvel_smtiOKN, my_interp(borderPos,5), my_interp(borderPos_lb,5) ] );
                    
                case 'mixed'
                    tmp = []; tmp_rResp = []; tmp_lResp = [];
                    for i=1:20
                        tmp = [ tmp, Run(irun).block(iblk).mixedLabelDesign(i), Run(irun).block(iblk).mixedLabelDesign(i), 0 ];
                        tmp_rResp = [tmp_rResp, Run(irun).block(iblk).rResp(i,:)];
                        tmp_lResp = [tmp_lResp, Run(irun).block(iblk).lResp(i,:)];
                    end
                    cond = my_interp( tmp, 300 );
                    rResp = my_interp( tmp_rResp, 5 );
                    lResp = my_interp( tmp_lResp, 5 );
                    
                    IM_mat = cat(1, IM_mat, [repmat(subNum,18000,1), repmat(state,18000,1), time', cond, lResp, rResp, iOKN, ...
                        vel_smtiOKN, smtvel_smtiOKN] );                        
            end
            clear borderPos_lb
        end
    end

    
    fileName = [ DIR.excel '/'  subID{isub} '_CT.xlsx'];
    xlswrite(fileName, CT_mat)
    
    fileName = [ DIR.excel '/'  subID{isub} '_IM.xlsx'];
    xlswrite(fileName, IM_mat)
    
    save( [DIR.excel '/' subID{isub} '.mat'], 'CT_mat', 'IM_mat')
    
    if downsamp % down sample 300Hz to 60Hz
        
        x = 1:size(IM_mat,1)/5;
        y = 5*x;
        
        IM_altResp = IM_mat(:, 6) - IM_mat(:, 5);
        d_IM_mat = [ IM_mat(y,1:4), IM_altResp(y), IM_mat( y, 11 ) ];
        
        
        x = 1:size(CT_mat,1)/5;
        y = 5*x;
        
        CT_altResp = CT_mat(:, 6) - CT_mat(:, 5);
        d_CT_mat = [ CT_mat(y,1:4), CT_altResp(y), CT_mat( y, [11,14] ) ];
        
        fileName = [ DIR.excel '/'  subID{isub} '_CT_60Hz.xlsx'];
        xlswrite(fileName, d_CT_mat)
        % CT_mat: 1)subNum, 2)state, 3)time, 4)stimCondition, 5)lResp, 6)rResp, 7:8)iOKN, ...
        % 9:10)vel_smtiOKN, 11:12)smtvel_smtiOKN, 13)borderPos, 14)borderPos_lb ] );
        
        fileName = [ DIR.excel '/'  subID{isub} '_IM_60Hz.xlsx'];
        xlswrite(fileName, d_IM_mat)
        
        save( [DIR.excel '/' subID{isub} '_60Hz.mat'], 'd_CT_mat', 'd_IM_mat')
        
        
    end
    
end

function y = my_interp(vector, factor)

y = [];
for m = 1 : length(vector)
    
    x = vector (m);
    y = cat(2, y, repmat(x, 1, factor));
    
end

y = y';
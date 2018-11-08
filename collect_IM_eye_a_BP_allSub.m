function collect_IM_eye_a_BP_allSub(DIR,subID,Cond)

whicheye = lower(Cond.whicheye);
time = linspace(-1,2,900);
switch whicheye
    case 'left'
        whichEyeColumn = 3;
    case 'right'
        whichEyeColumn = 5;
end

time2cut = Cond.time2cut * 10^3;

% savedDataName = [DIR.slowPhase '/ eyeAllSub.mat'];
% % FOR RAW DATA COLLECTION. VELOCITY ETC ARE BELOW %%%%%%%%%%%%%%%%%%%%%%%%
% jepc = 0;
% for isub = 1: length(subID)
%
%     subName = subID{isub};
%     subNum = str2num(subName(2:end));
%     loadFileName = [ subName '_intp' num2str(time2cut) '_vel_slwPhs_epIM.mat' ];
%     msg = ['collect_IM_raw_eye_allSub: Loading ' subName];
%     disp(msg)
%     load( [DIR.epoched '/' loadFileName])
%
%     for irun = 1 : length(Run)
%
%         for iblk = 1 : length(Run(irun).block)
%
%             if ~strcmp(Run(irun).block(iblk).condition,'mixed')
%                 continue
%             end
%
%             switch Run(irun).block(iblk).state
%                 case 'nc'
%                     state= 1; %=> block level
%                 case 'mon'
%                     state = 2;
%                 case 'mof'
%                     state = 3;
%                 case 'donmon'
%                     state = 4;
%                 case 'dofmon'
%                     state = 5;
%                 case 'donmof'
%                     state = 4;
%                 case 'dofmof'
%                     state = 5;
%             end
%
%             if isempty(Run(irun).block(iblk).gazeData)
%
%                 baseInfo(1:20,1) = str2num(subName(2:end));
%                 baseInfo(1:20,2) = irun;
%                 baseInfo(1:20,3) = iblk;
%                 baseInfo(1:20,4) = 1:20;
%                 baseInfo(1:20,5) = nan; %label
%                 baseInfo(1:20,6) = state;
%
%                 rawGaze(jepc + 1 : jepc + 20, : )  = [ baseInfo, nan(20,900)];
%
%                 jepc = jepc + 20;
%                 clear baseInfo
%                 continue
%             end
%
%             for iepc = 1: length(Run(irun).block(iblk).epoch)
%                 jepc = jepc + 1;
%
%                 label= Run(irun).block(iblk).epoch(iepc).labelDes; %=> epoch lavel. % 1)br, 2)red/right unambiguous, 3)green/left unambiduous
%
%                 baseInfo(1) = str2num(subName(2:end));
%                 baseInfo(2) = irun;
%                 baseInfo(3) = iblk;
%                 baseInfo(4) = iepc;
%                 baseInfo(5) = label;
%                 baseInfo(6) = state;
%
%                 rawGaze(jepc, 1:6 )  = baseInfo;
%
%                 ---------------------------------------------------------------
%                 1)subID, 2)irun, 3)iblk, 4)iepc, 5)stimLabel, 6)subjestState  |
%                 ---------------------------------------------------------------
%
%                 if isempty(Run(irun).block(iblk).epoch(iepc).iOKN_slwPhsDegPerSec_box)
%                     slwPhsBox(jepc, 7:906) = nan(1,900);
%                     slwPhsMed(jepc, 7:906) = nan(1,900);
%                     iOKN(jepc,7:906)       = nan(1,900);
%                     vel(jepc,7:906)        = nan(1,900);
%
%                 else
%                     slwPhsBox(jepc, 7:906) = Run(irun).block(iblk).epoch(iepc).iOKN_slwPhsDegPerSec_box(:,whichEyeColumn,1)'; %<= Cond.slowPhaseSmoothWin
%                     slwPhsMed(jepc, 7:906) = Run(irun).block(iblk).epoch(iepc).iOKN_slwPhsDegPerSec_med(:,whichEyeColumn,1)'; %<= Cond.slowPhaseSmoothWin
%                     iOKN(jepc, 7:906)      = Run(irun).block(iblk).epoch(iepc).iOKN(:,whichEyeColumn)'; %<= Cond.slowPhaseSmoothWin
%                     vel(jepc, 7:906)       = Run(irun).block(iblk).epoch(iepc).iOKN_vel(:,whichEyeColumn)'; %<= Cond.slowPhaseSmoothWin
%
%                 end
%
%                 if isempty(Run(irun).block(iblk).epoch(iepc).lResp_300Hz)||isempty(Run(irun).block(iblk).epoch(iepc).rResp_300Hz)
%                     lResp(jepc, 7:906) = nan(1,900);
%                     rResp(jepc, 7:906) = nan(1,900);
%                     resp(jepc, 7:906)  = nan(1,900);
%
%                     lBPBox(jepc, 7:906) = nan(1,900);
%                     rBPBox(jepc, 7:906) = nan(1,900);
%                     BPBox(jepc, 7:906)  = nan(1,900);
%
%                     lBPMed(jepc, 7:906) = nan(1,900);
%                     rBPMed(jepc, 7:906) = nan(1,900);
%                     BPMed(jepc, 7:906)  = nan(1,900);
%
%                 else
%                     lResp(jepc, 7:906) = Run(irun).block(iblk).epoch(iepc).lResp_300Hz'; %<= Cond.slowPhaseSmoothWin
%                     rResp(jepc, 7:906) = Run(irun).block(iblk).epoch(iepc).rResp_300Hz'; %<= Cond.slowPhaseSmoothWin
%                     resp(jepc, 7:906)  = Run(irun).block(iblk).epoch(iepc).resp'; %<= Cond.slowPhaseSmoothWin
%
%                     lBPBox(jepc, 7:906) = Run(irun).block(iblk).epoch(iepc).lBPBox'; %<= Cond.slowPhaseSmoothWin
%                     rBPBox(jepc, 7:906) = Run(irun).block(iblk).epoch(iepc).rBPBox'; %<= Cond.slowPhaseSmoothWin
%                     BPBox(jepc, 7:906)  = Run(irun).block(iblk).epoch(iepc).BPBox'; %<= Cond.slowPhaseSmoothWin
%
%                     lBPMed(jepc, 7:906) = Run(irun).block(iblk).epoch(iepc).lBPMed'; %<= Cond.slowPhaseSmoothWin
%                     rBPMed(jepc, 7:906) = Run(irun).block(iblk).epoch(iepc).rBPMed'; %<= Cond.slowPhaseSmoothWin
%                     BPMed(jepc, 7:906)  = Run(irun).block(iblk).epoch(iepc).BPMed'; %<= Cond.slowPhaseSmoothWin
%
%                     %% CAUTION!!!!! "+1" issue is already fixed in the previous process
%                     but to avoid run the code again, I put "+1" for here.
%                 end
%
%             end
%
%         end
%     end
% end
%
% save( [DIR.allSub '/IM_eyeAllSub.mat'] ,'iOKN', 'vel', 'slwPhsBox', 'slwPhsMed','subID')
% save( [DIR.allSub '/IM_respAllSub.mat'] , 'lResp', 'rResp', 'resp', 'BPBox', 'lBPBox', 'rBPBox', 'BPMed', 'lBPMed', 'rBPMed', 'subID')


jepc = 0;

for isub = 1: length(subID)
    
    subName = subID{isub};
    subNum = str2num(subName(2:end));
    loadFileName = [ subName '_intp' num2str(time2cut) '_vel_slwPhs_epIM.mat' ];
    msg = ['collect_IM_gaze_a_BP_allSub: Loading ' subName];
    disp(msg)
    load( [DIR.epoched '/' loadFileName])
    
    subName = subID{isub};
    subNum = str2num(subName(2:end));
    
    for irun = 1 : length(Run)
        
        for iblk = 1 : length(Run(irun).block)
            
            if ~strcmp(Run(irun).block(iblk).condition,'mixed')
                continue
            end
            
            
            %             if ~isfield(Run(irun).block(iblk), 'epoch')...
            %                     || ~isfield(Run(irun).block(iblk), 'mixedLabelDesign')
            %                 continue
            %             end
            %
            %             if isempty(Run(irun).block(iblk).epoch)
            %                 continue
            %             end
            
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
                
                baseInfo(1:20,1) = str2num(subName(2:end));
                baseInfo(1:20,2) = irun;
                baseInfo(1:20,3) = iblk;
                baseInfo(1:20,4) = 1:20;
                baseInfo(1:20,5) = nan; %label
                baseInfo(1:20,6) = state;
                
                lResp(jepc + 1 : jepc + 20, : )  = [ baseInfo, nan(20,900)];
                rResp(jepc + 1 : jepc + 20, : )  = [ baseInfo, nan(20,900)];
                resp (jepc + 1 : jepc + 20, : )  = [ baseInfo, nan(20,900)];
                lBPBox(jepc + 1 : jepc + 20, : ) = [ baseInfo, nan(20,900)];
                rBPBox(jepc + 1 : jepc + 20, : ) = [ baseInfo, nan(20,900)];
                BPBox(jepc + 1 : jepc + 20, : )  = [ baseInfo, nan(20,900)];
                lBPMed(jepc + 1 : jepc + 20, : ) = [ baseInfo, nan(20,900)];
                rBPMed(jepc + 1 : jepc + 20, : ) = [ baseInfo, nan(20,900)];
                BPMed(jepc + 1 : jepc + 20, : )  = [ baseInfo, nan(20,900)];
                
                iOKN(jepc + 1 : jepc + 20, : )      = [ baseInfo, nan(20,900)]; % iOKN
                vel(jepc + 1 : jepc + 20, : )       = [ baseInfo, nan(20,900)]; % velocity
                slwPhsBox(jepc + 1 : jepc + 20, : ) = [ baseInfo, nan(20,900)]; % slow phase smoothed by boxcar
                slwPhsMed(jepc + 1 : jepc + 20, : ) = [ baseInfo, nan(20,900)]; % slow phase smoothed by median filtering
                
                jepc = jepc + 20;
                clear baseInfo
                continue
            end
            
            for iepc = 1: length(Run(irun).block(iblk).epoch)
                jepc = jepc + 1;
                
                label= Run(irun).block(iblk).epoch(iepc).labelDes; %=> epoch level. % 1)br, 2)red/right unambiguous, 3)green/left unambiduous
                
                baseInfo(1) = str2num(subName(2:end));
                baseInfo(2) = irun;
                baseInfo(3) = iblk;
                baseInfo(4) = iepc;
                baseInfo(5) = label;
                baseInfo(6) = state;
                
                lResp(jepc, 1:6 )  = baseInfo;
                rResp(jepc, 1:6 )  = baseInfo;
                resp (jepc, 1:6 )  = baseInfo;
                lBPBox(jepc, 1:6 ) = baseInfo;
                rBPBox(jepc, 1:6 ) = baseInfo;
                BPBox(jepc, 1:6 )  = baseInfo;
                lBPMed(jepc, 1:6 ) = baseInfo;
                rBPMed(jepc, 1:6 ) = baseInfo;
                BPMed(jepc, 1:6 )  = baseInfo;
                
                
                iOKN(jepc, 1:6 )      = baseInfo; % iOKN
                vel(jepc, 1:6 )       = baseInfo; % velocity
                slwPhsBox(jepc, 1:6 ) = baseInfo; % slow phase smoothed by boxcar
                slwPhsMed(jepc, 1:6 ) = baseInfo; % slow phase smoothed by median filtering
                
                
                % ---------------------------------------------------------------
                % 1)subID, 2)irun, 3)iblk, 4)iepc, 5)stimLabel, 6)subjestState  |
                % ---------------------------------------------------------------
                
                if isempty(Run(irun).block(iblk).epoch(iepc).iOKN_slwPhsDegPerSec_box)
                    slwPhsBox(jepc, 7:906) = nan(1,900);
                    slwPhsMed(jepc, 7:906) = nan(1,900);
                    iOKN(jepc,7:906)       = nan(1,900);
                    vel(jepc,7:906)        = nan(1,900);
                    
                else
                    slwPhsBox(jepc, 7:906) = Run(irun).block(iblk).epoch(iepc).iOKN_slwPhsDegPerSec_box(:,whichEyeColumn,1)'; %<= Cond.slowPhaseSmoothWin
                    slwPhsMed(jepc, 7:906) = Run(irun).block(iblk).epoch(iepc).iOKN_slwPhsDegPerSec_med(:,whichEyeColumn,1)'; %<= Cond.slowPhaseSmoothWin
                    iOKN(jepc, 7:906)      = Run(irun).block(iblk).epoch(iepc).iOKN(:,whichEyeColumn)'; %<= Cond.slowPhaseSmoothWin
                    vel(jepc, 7:906)      = Run(irun).block(iblk).epoch(iepc).iOKN_vel(:,whichEyeColumn)'; %<= Cond.slowPhaseSmoothWin
                    
                end
                
                if isempty(Run(irun).block(iblk).epoch(iepc).lResp_300Hz)||isempty(Run(irun).block(iblk).epoch(iepc).rResp_300Hz)
                    lResp(jepc, 7:906) = nan(1,900);
                    rResp(jepc, 7:906) = nan(1,900);
                    resp(jepc, 7:906)  = nan(1,900);
                    
                    lBPBox(jepc, 7:906) = nan(1,900);
                    rBPBox(jepc, 7:906) = nan(1,900);
                    BPBox(jepc, 7:906)  = nan(1,900);
                    
                    lBPMed(jepc, 7:906) = nan(1,900);
                    rBPMed(jepc, 7:906) = nan(1,900);
                    BPMed(jepc, 7:906)  = nan(1,900);
                    
                else
                    lResp(jepc, 7:906) = Run(irun).block(iblk).epoch(iepc).lResp_300Hz'; %<= Cond.slowPhaseSmoothWin
                    rResp(jepc, 7:906) = Run(irun).block(iblk).epoch(iepc).rResp_300Hz'; %<= Cond.slowPhaseSmoothWin
                    resp(jepc, 7:906)  = Run(irun).block(iblk).epoch(iepc).resp'; %<= Cond.slowPhaseSmoothWin
                    
                    lBPBox(jepc, 7:906) = Run(irun).block(iblk).epoch(iepc).lBPBox'; %<= Cond.slowPhaseSmoothWin
                    rBPBox(jepc, 7:906) = Run(irun).block(iblk).epoch(iepc).rBPBox'; %<= Cond.slowPhaseSmoothWin
                    BPBox(jepc, 7:906)  = Run(irun).block(iblk).epoch(iepc).BPBox'; %<= Cond.slowPhaseSmoothWin
                    
                    lBPMed(jepc, 7:906) = Run(irun).block(iblk).epoch(iepc).lBPMed'; %<= Cond.slowPhaseSmoothWin
                    rBPMed(jepc, 7:906) = Run(irun).block(iblk).epoch(iepc).rBPMed'; %<= Cond.slowPhaseSmoothWin
                    BPMed(jepc, 7:906)  = Run(irun).block(iblk).epoch(iepc).BPMed'; %<= Cond.slowPhaseSmoothWin
                    
                    %%% CAUTION!!!!! "+1" issue is already fixed in the previous process
                    % but to avoid run the code again, I put "+1" here.
                end
                
            end
            
        end
    end
end

save( [DIR.allSub '/IM_eyeAllSub.mat'] ,'iOKN', 'vel', 'slwPhsBox', 'slwPhsMed','subID')
save( [DIR.allSub '/IM_respAllSub.mat'] , 'lResp', 'rResp', 'resp', 'BPBox', 'lBPBox', 'rBPBox', 'BPMed', 'lBPMed', 'rBPMed', 'subID')

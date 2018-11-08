function compute_vel_a_slwPhs_a_epIM(DIR, subID, Cond)

% == 2016-4-30 MF ===
% Put computing steps in preprocess all together to save time as locading
% takes a lot of time

whicheye = Cond.whicheye;
time2cut = Cond.time2cut;
epochWin = Cond.epochWinIM;

%How to get OKN slow phase?
medFilt = 0;

s = Cond.velocitySmoothPnt;

%%%% for assesment, set default later %%%

% t = [60 3];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

i = 0;

for isub = 1%: length(subID)
    
    subName = subID{isub};
    subNum = str2num(subName(2:end));
    
%     if ~ismember( subNum, [127, 130, 131, 149] )
%         continue
%     end

    loadFileName = [ subID{isub} '_intp' num2str(time2cut*10^3) '.mat' ];
    saveFileName = [ DIR.epoched '/' loadFileName(1:end-4) '_vel_slwPhs_epIM.mat'];
    
%     if exist( saveFileName, 'file' )
%         continue
%     end
        
    load( [ DIR.interpolated '/' loadFileName] )
    load( [DIR.cfg '/' subID{isub} '_Cfg_a_br.mat'] )
    
    %% compute velocity with boxcar smoothing =============================
    for irun = 1 : length(Run)
        for iblk = 1: length(Run(irun).block)
            
            if isempty(Run(irun).block(iblk).gazeData),
                
                switch Run(irun).block(iblk).state
                    case 'nc'
                        state = 1;
                    case 'mon'
                        state = 2;
                    case 'mof'
                        state = 3;
                    case {'don','donmon', 'donmof'}
                        state = 4;
                    case {'dof','dofmon', 'dofmof'}
                        state = 5;
                end
                i = i + 1;
                skipInfo(i,:) = [ state, str2num(subID{isub}(2:4)), irun, iblk ];
                msg = sprintf('compute_velocity. Skipped sub: %s, run: %s, block: %s', subID{isub}, num2str(irun), num2str(iblk) );
                disp(msg);
                continue
            end
                        
            %%% Preallocation %%%
            Run(irun).block(iblk).iOKN_vel = Run(irun).block(iblk).iOKN; % preallocate % in pixel.
            Run(irun).block(iblk).iOKN_slwPhsDegPerSec_box = Run(irun).block(iblk).iOKN;
            Run(irun).block(iblk).iOKN_slwPhsDegPerSec_med = Run(irun).block(iblk).iOKN;

            % == Explanation of the 17 columns ===
            % 1. time in sec (only sec part)
            % 2. time in msec (only miclosec part)
            % 3. x gaze coordinate of the left eye
            % 4. y gaze coordinate of the left eye
            % 5. x gaze coordinate of the right eye
            % 6. y gaze coordinate of the right eye
            % 7. left eye position - x coordinate
            % 8. left eye position - y coordinate
            % 9. right eye position - x coordinate
            % 10. right camera eye position - y coordinate
            % 11. left eye validity
            % 12. right eye validity
            % 13. diameter of pupil of the left eye
            % 14. diameter of pupil of the right eye
            % 15. distance of the camera from the left eye
            % 16. distance of the camera from the right eye
            % 17. time in sec(including msec part)
            
            
            iOKN = Run(irun).block(iblk).iOKN(:,3:6);
            
            %% BOXCAR SMOOTHING OF iOKN ===================================
            %             x = tsmovavg(iOKN,'s', s);
            siOKN = boxcar(iOKN',30);
            
            
            %% COMPUTE VELOCITY ===========================================
            vel = diff(siOKN') * Cfg.degPerPixel * Cfg.samplingRate; % velocity of smoothed iOKN ( in deg/sec )
            Run(irun).block(iblk).iOKN_vel(2:end,3:6) = vel; %- time x eye x s x p
            
            %% SMOOTHING TO GET SLOW PHASE ================================
            % 1) Boxcar filtering
            slwPhsBox = boxcar(vel', s); % ( in deg/sec )
            
            % 2) Median filtering
            slwPhsMed = medfilt1(vel,s); % ( in deg/sec )
            
            
            Run(irun).block(iblk).iOKN_slwPhsDegPerSec_box(2:end,3:6) = slwPhsBox'; %- time x eye x sv
            Run(irun).block(iblk).iOKN_slwPhsDegPerSec_med(2:end,3:6) = slwPhsMed; %- time x eye x sv
            
            time = clock;
            msg = sprintf('compute_slowPhase. Done with sub: %s, run: %s, block: %s. %d-%d-%d. %d:%d:%d', subID{isub}, num2str(irun), num2str(iblk), round(time));
            disp(msg);
            

            %% Epoch eye data based on time ================================
            
            if ~strcmp(Run(irun).block(iblk).id,'2-1')
                continue
            end
            
            if isempty(Run(irun).block(iblk).gazeData)
                continue
            end
            
            epochTimes(:,1) = Run(irun).block(iblk).trial_evTimes(:,1) + repmat(epochWin(1), Run(irun).block(iblk).nTrials, 1);
            epochTimes(:,2) = epochTimes(:,1) + ( diff(epochWin) ) ;
            
            for iwhat = 1:4
                
                switch iwhat
                    case 1 %iOKN
                        data = Run(irun).block(iblk).iOKN;
                    case 2 %iOKN Slow Phase smoothed by box
                        data =  squeeze( Run(irun).block(iblk).iOKN_slwPhsDegPerSec_box );
                    case 3 %iOKN Slow Phase smoothed by med
                        data =  squeeze( Run(irun).block(iblk).iOKN_slwPhsDegPerSec_med );
                    case 4
                        data =  squeeze( Run(irun).block(iblk).iOKN_vel );
                end
                
                if isempty(data), continue, end
                
                for itri = 1: Run(irun).block(iblk).nTrials
                    Run(irun).block(iblk).epoch(itri).labelDes = Run(irun).block(iblk).mixedLabelDesign(itri);
                    
                    data_epoch = data(data(:,17,1) >= epochTimes(itri, 1) & data(:,17) < epochTimes(itri, 2), : );
                    
                    if size(data_epoch,1) <  diff(epochWin)*300
                        continue
                        
                    elseif size(data_epoch,1) >  diff(epochWin)*300
                        data_epoch = data_epoch( 1 : diff(epochWin)*300, : , :);
                    end
                    
                    switch iwhat
                        case 1
                            Run(irun).block(iblk).epoch(itri).iOKN = data_epoch;
                        case 2
                            Run(irun).block(iblk).epoch(itri).iOKN_slwPhsDegPerSec_box = data_epoch;
                        case 3
                            Run(irun).block(iblk).epoch(itri).iOKN_slwPhsDegPerSec_med = data_epoch;
                        case 4
                            Run(irun).block(iblk).epoch(itri).iOKN_vel = data_epoch;
                    end
                    
                end
                
            end
            
            %% Smooth button press ========================================
            %%% RAW button press up sampled%%%
            rResp_300Hz = [];
            lResp_300Hz = [];
            for itri = 1: Run(irun).block(iblk).nTrials
                rResp_300Hz = cat( 1, rResp_300Hz, Run(irun).block(iblk).trial(itri).rResp_300Hz);
                lResp_300Hz = cat( 1, lResp_300Hz, Run(irun).block(iblk).trial(itri).lResp_300Hz);
            end
            
            resp = sum(cat(1, -1*lResp_300Hz', rResp_300Hz'), 1);
            
            %%% Boxcar smoothing %%%
            lBPBox = boxcar(lResp_300Hz',s)';
            rBPBox = boxcar(rResp_300Hz',s)';
            BPBox  = boxcar(resp,s);
                        
            %%% Median filtered button press %%%
            lBPMed = medfilt1(lResp_300Hz, s);
            rBPMed = medfilt1(rResp_300Hz, s); 
            BPMed  = medfilt1(resp, s); 
            
            %% Make a list of time window for epoching behavioural report
            % (Behaviour doesn't have time with it)
            winLength = diff(epochWin) * 300;
            epochWin_300Hz(1,:) = epochWin * 300;
            
            vepoch = 1;
            while epochWin_300Hz(vepoch,2) < length(rResp_300Hz)
                vepoch = vepoch + 1;
                epochWin_300Hz(vepoch, 1) = epochWin_300Hz(vepoch - 1, 2) + 1;
                epochWin_300Hz(vepoch, 2) = epochWin_300Hz(vepoch - 1, 2) + winLength;
            end
            
            for iepc = 1: size(epochWin_300Hz,1)
                if epochWin_300Hz(iepc, 1) < 0 || epochWin_300Hz(iepc, 2) > size(rResp_300Hz,1)
                    continue
                end
                Run(irun).block(iblk).epoch(iepc).lResp_300Hz = lResp_300Hz(epochWin_300Hz(iepc, 1):epochWin_300Hz(iepc, 2));
                Run(irun).block(iblk).epoch(iepc).rResp_300Hz = rResp_300Hz(epochWin_300Hz(iepc, 1):epochWin_300Hz(iepc, 2));
                Run(irun).block(iblk).epoch(iepc).resp        = resp(epochWin_300Hz(iepc, 1):epochWin_300Hz(iepc, 2));
                Run(irun).block(iblk).epoch(iepc).lBPBox      = lBPBox(epochWin_300Hz(iepc, 1):epochWin_300Hz(iepc, 2));
                Run(irun).block(iblk).epoch(iepc).rBPBox      = rBPBox(epochWin_300Hz(iepc, 1):epochWin_300Hz(iepc, 2));
                Run(irun).block(iblk).epoch(iepc).BPBox       = BPBox(epochWin_300Hz(iepc, 1):epochWin_300Hz(iepc, 2));
                Run(irun).block(iblk).epoch(iepc).lBPMed      = lBPMed(epochWin_300Hz(iepc, 1):epochWin_300Hz(iepc, 2));
                Run(irun).block(iblk).epoch(iepc).rBPMed      = rBPMed(epochWin_300Hz(iepc, 1):epochWin_300Hz(iepc, 2));
                Run(irun).block(iblk).epoch(iepc).BPMed       = BPMed(epochWin_300Hz(iepc, 1):epochWin_300Hz(iepc, 2));
                
            end
            
            Run(irun).block(iblk).epoch_evIdxs = epochWin_300Hz(~isnan(epochWin_300Hz(:,1)),:);
            
            time = clock;
            msg = sprintf('epoch_IM. Done with sub: %s, run: %s, block: %s. %d-%d-%d. %d:%d:%d', subID{isub}, num2str(irun), num2str(iblk), round(time));
            disp(msg);
            
        end
    end
    
    % save new file
    save( saveFileName,'Run', 'Cond', '-v7.3')
    
end

function y = my_interp(vector, factor)

y = [];
for m = 1 : length(vector)
    
    x = vector (m);
    y = cat(2, y, repmat(x, 1, factor));
    
end

y = y';
function interpolate_missingPoints(DIR,subID,Cond)

% 12 June 2015, MF
% Filter data: cut saccades and blinks from the gaze data and then
% concatenate and interpolate the eye position   
% Eliminate saccades & Eliminate Blinks

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

whicheye = Cond.whicheye;
time2cut = Cond.time2cut;

for isub = 1: length(subID)
    
    subName = subID{isub};
    subNum = str2num(subName(2:end));
    
%     if ~ismember( subNum, [127, 130, 131, 149] )
%         continue
%     end
    
    loadFileName = [DIR.collectedData '/' subID{isub} '.mat'];
    load( loadFileName )

    for irun = 1 : length(Run)
        
        for iblk = 1 : length(Run(irun).block)
            tMissP = [];
            
            gazeData = Run(irun).block(iblk).gazeData;
           
            if isempty(gazeData), continue, end
            
            blinks = Run(irun).block(iblk).blinks;
            allsac = Run(irun).block(iblk).allsac;
            
            if ~isempty(blinks)
                if ~isnan(blinks)
                    tMissP = blinks ;
                end
            end

            if ~isempty(Run(irun).block(iblk).allsac)
                tMissP = [tMissP; allsac(:, [1 2]) ];
            end
             
            %% adjust the length so that they can be comparable with button press responses
            nPoints = (Run(irun).block(iblk).stimFrames + ...
                    Run(irun).block(iblk).ISI ) * Run(irun).block(iblk).nTrials * 5;
                
            if length(gazeData) > nPoints
                gazeData = gazeData(1:nPoints, :);
            end
            
            rawGazeData = gazeData;
            
            %% replace missing points with nan
            tMissP = sortrows ( tMissP ); 
            tMissP2cut = [ tMissP(:,1) - time2cut, tMissP(:,2) + time2cut ];
        
            for iMissP = 1 : size(tMissP2cut,1)
                gazeData ( gazeData ( : , 17 ) >= tMissP2cut (iMissP,1)...
                    & gazeData ( : , 17 ) <= tMissP2cut (iMissP,2) , 3:10 ) = nan;
            end
            
            %% get index for missinp points (Use index instead of time
            % because there are many ovelaps in the list tMissP2cut made by
            % time width for detecting missing points.)
            idxMissP2cut = find( isnan ( gazeData(:,3) ) );
            nSamp = length(gazeData(:,1));
            
            idxMissP2cut2(1,1) = 1;
            
            missP_n0 = 0;
            imissP = 0; 
            for iSamp = 1 : nSamp
                
                if iSamp >= idxMissP2cut(1)
                    
                    if ~ismember( iSamp, idxMissP2cut)
                        missP_n1 = 0;
                    else
                        missP_n1 = 1;
                    end
                    
                    missP_switch = missP_n0 - missP_n1;
                    %- this variable would be -1 when it gets into missing
                    %points period.
                    
                    if missP_switch == -1
                        %- if it gets into missing points period
                        imissP = imissP + 1;
                        idxMissP2cut2(imissP, 1) = iSamp;
     
                    elseif missP_switch == 1
                        %- if it gets out from missing points period
                        idxMissP2cut2(imissP, 2) = iSamp - 1;
                    end
                    
                    missP_n0 = missP_n1;
                end
            end
            
            if idxMissP2cut2(end,2) == 0
                idxMissP2cut2(end,2) = nSamp;
                % in case the block finised with missing point
            end
            
            %% move up or down data after missing points
            
%             %- case gazeData is longer than expected length, cut off the
%             %end to make it comparable.
%             if length(gazeData) > Run(irun).block(iblk).stimFrames * Run(irun).block(iblk).nTrials * 5
%                 gazeData = gazeData(1:Run(irun).block(iblk).stimFrames * Run(irun).block(iblk).nTrials * 5, :); 
%             end
            
            
            %%% preallocation
            filtGazeData = gazeData; %- missing points are replaced by NaN
            intpGazeData = gazeData; %- missing points are interpolated
            iOKN         = gazeData; %- missing points are interpolated and the data is integrated
            
            for vidxMissP2cut = 1 : size(idxMissP2cut2,1)
                cutStartIdx = idxMissP2cut2(vidxMissP2cut,1);
                cutEndIdx   = idxMissP2cut2(vidxMissP2cut,2);
                
                if cutEndIdx == length(iOKN) %- if the index gets to the end
                    gap = iOKN(cutStartIdx - 1, 3:10) - iOKN(cutEndIdx, 3:10);
                    
                elseif cutStartIdx == 1 %- if gazeData starts with missing point
                    gap = 0 - iOKN(cutEndIdx + 1, 3:10);
                    
                elseif cutEndIdx
                    gap = iOKN(cutStartIdx - 1, 3:10) - iOKN(cutEndIdx + 1, 3:10);
                end
                
                iOKN( (cutEndIdx + 1) : length(iOKN), 3:10 ) = ...
                    iOKN( (cutEndIdx + 1) : length(iOKN), 3:10 ) + repmat ( gap, length(iOKN) - cutEndIdx , 1);
                
            end

            %- missing points have been replaced by nan and gaze data after
            % missing points were moved up or down so that we can get
            % integrated OKN.
            
            
            %% interpolate linearly (replace nan)
            for icolumn = 3:10 %- eye gaze data
                x(:,icolumn) = intpGazeData(:,icolumn);
                t = linspace( 0, 2, numel(x(:,icolumn)) );
                nans = isnan(x(:,icolumn));
                x(nans,icolumn) = interp1( t(~nans), x(~nans,icolumn), t(nans) );
            end
            
            intpGazeData(:,3:10) = x(:,3:10);
            
            for icolumn = 3:10 %- eye gaze data
                x(:,icolumn) = iOKN(:,icolumn);
                t = linspace( 0, 2, numel(x(:,icolumn)) );
                nans = isnan(x(:,icolumn));
                x(nans,icolumn) = interp1( t(~nans), x(~nans,icolumn), t(nans) );
            end
            iOKN(:,3:10) =x(:,3:10);
            
            %- missing points replaced by nan have been inrerpolated
            %linearly.
            
         
            %% put the data into structure            
            for itri = 1 : Run(irun).block(iblk).nTrials
                
                % Filter the blinks and saccades from the data in each trial
                tr_times = Run(irun).block(iblk).trial_evTimes(itri, :);
                rawGazeData_tri  = rawGazeData(rawGazeData(:,17) >= tr_times(1) & rawGazeData(:,17) <= tr_times(2), : );
                filtGazeData_tri = filtGazeData(filtGazeData(:,17) >= tr_times(1) & filtGazeData(:,17) <= tr_times(2), : );
                intpGazeData_tri = intpGazeData(intpGazeData(:,17) >= tr_times(1) & intpGazeData(:,17) <= tr_times(2), : );
                iOKN_tri         = iOKN(iOKN(:,17) >= tr_times(1) & iOKN(:,17) <= tr_times(2), : );
             
                
                % size of gazeData might be bigger or smaller than expected size.
                % So, check the size and adjust into the same size in each
                % trial.
                nSamps = Run(irun).block(iblk).stimFrames * 5;
                if size(iOKN_tri,1) >= nSamps
                    rawGazeData_tri  = rawGazeData_tri(1: nSamps, :);
                    filtGazeData_tri = filtGazeData_tri(1: nSamps, :);
                    intpGazeData_tri = intpGazeData_tri(1: nSamps, :);
                    iOKN_tri         = iOKN_tri(1: nSamps, :);
                    
                    idx_missp_tri = find ( isnan ( filtGazeData ...
                        ( filtGazeData(:,17) >= tr_times(1) & filtGazeData(:,17) <= tr_times(2), 3 ) ) );
                    
                elseif  size(iOKN_tri,1) < nSamps
                    nS = size(iOKN_tri, 1);
                    idx_missp_tri = find ( isnan ( filtGazeData ...
                        ( filtGazeData(1:nS,17) >= tr_times(1) & filtGazeData(1:nS,17) <= tr_times(2), 3 ) ) );
                    
                    rawGazeData_tri(nS + 1 : nSamps, :)   = nan;
                    filtGazeData_tri(nS + 1 : nSamps, :)  = nan;
                    intpGazeData_tri(nS + 1 : nSamps, :)  = nan;
                    iOKN_tri(nS + 1 : nSamps, :)          = nan;
                end
                
                %% culcurate percentage of missing points
                filtRateRemovedP(itri,1) = sum( isnan( filtGazeData_tri( :,3 ) ) ) / length( filtGazeData_tri( :,3 ) ); %-left eye x coodinate
                filtRateRemovedP(itri,2) = sum( isnan( filtGazeData_tri( :,5 ) ) ) / length( filtGazeData_tri( :,5 ) ); %-right eye x coodinate
                
                rawRateRemovedP(itri,1) = length(find(rawGazeData_tri(1:nSamps ,11)==4)) / length( rawGazeData_tri( 1:nSamps ,3 ) ); %-left eye x coodinate
                rawRateRemovedP(itri,2) = length(find(rawGazeData_tri(1:nSamps ,12)==4)) / length( filtGazeData_tri( 1:nSamps ,5 ) ); %-right eye x coodinate
                
                
                %% Resample the response vector to match the eye data
                rResp_300hz = my_interp(Run(irun).block(iblk).rResp(itri, :), 5);
                lResp_300hz = my_interp(Run(irun).block(iblk).lResp(itri, :), 5);
                                
%                 response_times_300hz = (1: size(response_filt_300hz, 1) ) ./ 300 * 1000;
                
%                 response_filt_300hz(idx_missp_tri) = 2;
                %- as nan in response means "no bottun press", replace
                %response with "2" for gaze missing points. (nan: no bottun
                %press, 0: right, 1: left, 2: missing point )
                
                % Save data into Run
                Run(irun).block(iblk).trial(itri).rawGazeData  = rawGazeData_tri;
                Run(irun).block(iblk).trial(itri).filtGazeData = filtGazeData_tri;
                Run(irun).block(iblk).trial(itri).intpGazeData = intpGazeData_tri;
                Run(irun).block(iblk).trial(itri).iOKN         = iOKN_tri;   % this is in pixel             
                Run(irun).block(iblk).trial(itri).rResp_300Hz  = rResp_300hz;
                Run(irun).block(iblk).trial(itri).lResp_300Hz  = lResp_300hz;
                                                
                
            end
            Run(irun).block(iblk).timeMissP    = tMissP;
            Run(irun).block(iblk).idxMissP2cut = idxMissP2cut;
            Run(irun).block(iblk).rateRemovedP = filtRateRemovedP;
            Run(irun).block(iblk).iOKN         = iOKN;
            
            clear gazeData x idxMissP2cut2 rateRemovedP
        end
    end

    saveFileName = [ subID{isub} '_intp' num2str(time2cut*10^3) '.mat' ];
    save( [DIR.interpolated '/' saveFileName ], 'Run' )
    
    msg = sprintf('interpolate_missingPoints. Done with sub: %s', subID{isub}); disp(msg);
    
    clear Run
    
end

function y = my_interp(vector, factor)

y = [];
for m = 1 : length(vector)
    
    x = vector (m);
    y = cat(2, y, repmat(x, 1, factor));
    
end

y = y';


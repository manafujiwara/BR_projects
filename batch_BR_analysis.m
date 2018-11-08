function batch_BR_analysis(DIR,subID)
%% IMPORTANT ====================================================
% Run "[DIR, subID] = startAnalysis.m" before "batch_BR_analysis.m" to 
% obtain variables DIR and subID.

% <<<<<<< DATA DIMENSION >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% 1. Subject: NC(Normal Control), MON(Medication ON), MOF(Medication OFF),
%             DON(DBS ON), DOF(DBS OFF)
% 2. Data: RT(behaviour), OKN(eye movement)
% 3. Stimuli: Ambiguous, Unambiguous
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

Cond.whicheye = 'LEFT'; % or 'RIGHT'
% Right means calculate saccades and fixations only from the right eye data


%% ===== JOBS ====================================================
% if 1, the job runs, if 0, the job does not run.

% ------ Single file preproces --------------------------------------------
% Collect behavioral data and eye data for each run / subject.
% Save both Data into 'Run' structure in one file
Job.collect_Cfg_a_BR  = 0; %OK 
Job.collect_UPDRS     = 0; %OK -> Run by itself: loads excel file 'UPDRS.xlsx' and puts it in .mat file
Job.collect_data      = 0; %OK
Job.makeFiles4Sharing = 0; %OK To share collected data, make a stracture
Job.check_gazeData    = 0; %OK check RAW gaze data here. Blinks/missing gaze points are removed later in the process.
% Job.check_UPDRS = 0;    %-> Run by itself, open code and check it.

% ----- Filter eye data ---------------------------------------------------
% Remove saccades and blinks, then concatenate and interpolate in each file
Job.interpolate_missingPoints = 0; %OK
    Cond.time2cut = 0.01; %<- before and after blink/saccade/undetected points (sec)

% ----- Compute vewlocity and epoch ---------------------------------------
Job.compute_vel_a_slwPhs_a_epIM= 0; %compute velocity, slowphase and then epoch
    Cond.velocitySmoothType = 'boxcar'; % or 'gausian'
    Cond.velocitySmoothPnt  = [30]; % Need to be common divisor of 60 and 300
    Cond.convFiltWidth      = 0.1;  % in sec
    Cond.slowPhaseSmoothWin = [30 150]; % nPoints for smoothing
    Cond.epochWinIM         = [-1 2]; % from button stimulus onset (sec)

Job.exportData2excel = 0;
    
%% === INTERMITTENT CONDITION ====================================
Job.collect_IM_eye_a_BP_allSub = 0; %OK
Job.label_IM                   = 0; %OK label BP and eye => get "BPinfo", "eye30Info", "eye150info"
Job.rej_upsmp_smth_IM          = 0; %OK reject subjects, updample and smooth trials & plot LEFT/RIGHT LABEL SEPARATELY
Job.check_vel_a_slwPhs_a_epIM  = 0; %OK 
Job.check_IM_BPtime            = 0; %OK 

% === FOR REVISION ========================================================
Job.check_centeringGaze             = 0; %OK Histogram is presented only for Control vs ON-treatment PD, but stats is done for all the group comparisons
Job.decode_IM_eye_a_BP_varShift     = 0; %OK
Job.check_cmp_subGrp                = 0; %OK
Job.check_reachedDecAcc_absCri      = 0; %OK


%% === CONTINUOUS CONDITION ======================================
Job.compute_perceptTendency = 0;

% ----- Continuous condition ----------------------------------------------
Job.exportData2excel                       = 0; %OK
Job.reject_smooth_CN                       = 0; %OK
Job.decode_SVM_CN_eye_a_BP_allSub          = 0; %OK
Job.check_decode_SVM_CN_eye_a_BP_allSub    = 0; %OK







%% ===== RUN =====================================================
if Job.collect_Cfg_a_BR
    collect_Cfg_a_BR(DIR, subID)
end

if Job.collect_UPDRS
    [statSum_UPDRS] = collect_UPDRS  % statSum_UPDRS: [group, p, tstat, df, sd]
end


if Job.collect_data
    collect_data(DIR,subID,Cond)
end

if Job.check_gazeData
    check_gazeData(DIR,subID)
end

if Job.interpolate_missingPoints
    interpolate_missingPoints(DIR,subID,Cond)
end

if Job.compute_vel_a_slwPhs_a_epIM
    compute_vel_a_slwPhs_a_epIM(DIR,subID,Cond)
end

if Job.exportData2excel
    exportData2excel(DIR,subID,Cond)
end

if Job.compute_perceptTendency
    compute_perceptTendency(DIR,subID,Cond)
end

if Job.check_IM_eye_alSub
    check_IM_eye_alSub(DIR,subID,Cond)
end

if Job.collect_IM_eye_a_BP_allSub
    collect_IM_eye_a_BP_allSub(DIR,subID,Cond)
end

if Job.label_IM
    label_IM(DIR,Cond)
end

if Job.rej_upsmp_smth_IM
    rej_upsmp_smth_IM(DIR,subID)
end

if Job.check_vel_a_slwPhs_a_epIM 
    check_vel_a_slwPhs_a_epIM(DIR)
end

if Job.check_IM_BPtime
    check_IM_BPtime(DIR, subID);
end

if Job.check_centeringGaze
    check_centeringGaze(DIR, subID, Cond);
end

if Job.decode_IM_eye_a_BP_varShift
    decode_IM_eye_a_BP_varShift(DIR,subID);
end

if Job.check_cmp_subGrp 
    check_cmp_subGrp(DIR)
end

if Job.check_reachedDecAcc_absCri
    check_reachedDecAcc_absCri(DIR);
end


if Job.makeFiles4Sharing
    makeFiles4Sharing(DIR)
end

if Job.reject_smooth_CN
    reject_smooth_CN(DIR, subID)
end

if Job.decode_SVM_CN_eye_a_BP_allSub
    decode_SVM_CN_eye_a_BP_allSub(DIR, subID)
end

if Job.check_decode_SVM_CN_eye_a_BP_allSub
    check_decode_SVM_CN_eye_a_BP_allSub(DIR, subID)
end

if Job.decode_SVM_CN_eye_a_BP_allSub_cmpEyeBP
    decode_SVM_CN_eye_a_BP_allSub_cmpEyeBP(DIR,subID)
end

if Job.compute_mean_a_latencies
    compute_mean_a_latencies(DIR, subID)
end
   
if Job.decode_IM_eye_a_BP
    decode_IM_eye_a_BP(DIR, subID)
end

if Job.check_cmp_inR
    check_cmp_inR(DIR, subID)
end

if Job.check_corr_UPDRS_a_lat
    check_corr_UPDRS_a_lat(DIR)
end

end

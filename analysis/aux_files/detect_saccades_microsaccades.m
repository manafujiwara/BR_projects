function [saccade blink] = detect_saccades_microsaccades (eall, whicheye, whicheyeanalyze, subjid, PARAMS, DIR)




% Microsaccade detection using Engbert algorithm
% INPUT: eall from read_eyetracker_data
% SACCADE OUTPUT (allsac, microsac, macrosac)
%  sac(1:num,1)   onset of saccade (msec)
%  sac(1:num,2)   end of saccade (msec)
%  sac(1:num,3)   duration of saccade (msec)
%  sac(1:num,4)   x ini sacc (pixs - absolute (0,0)=(left,top))
%  sac(1:num,5)   y ini sacc (pixs - absolute (0,0)=(left,top))
%  sac(1:num,6)   x end sacc (pixs - absolute (0,0)=(left,top))
%  sac(1:num,7)   y end sacc (pixs - absolute (0,0)=(left,top))
%  sac(1:num,8)   amplitude (deg)
%  sac(1:num,9)   peak velocity of saccade (vpeak) (deg/sec)
%  sac(1:num,10)  saccade angle (deg)
%  sac(1:num,11)  alternative saccade amplitude (in deg)

% VERSION 30-12-09
% - options to decide which eye analize, and which eye is going to be the
% output (in the case of binocular data)
% - now I remove eall the operations refered to the EyeLink events, because
% it is working too slow

% VERSION 08-01-10
% - whicheye: which eye is going to be the output ('right' or 'left')
% - whicheyeanalize: which eye analize in the case of binocular data
% ('right', 'left' or 'both')
% - subjid: subject name to save the figures
% - doplot: 0: No figures (default), 1: Macro and Microsaccades figures

% VERSION 13/04/2010
% Sampling rate read from eall structure (since binocular data is acquired
% at 1000Hz

% VERSION 20/04/2010
% Previous versions were detecting RIGHT-LEFT but the output of saccades
% and fixations was always left (d(:,2:3)).


try
    
    
    doplot  = PARAMS.doplot.calc;
    dodownsample    = PARAMS.dodownsample;
    
    
    %% Saccade detection parameters
    
    if isfield(PARAMS,'engbertkliegl')
        VTHRES          = PARAMS.engbertkliegl.VTHRES;          % Velocity threshold
        SAMPLING        = PARAMS.engbertkliegl.SAMPLING;        % eall.srate;    % Sampling rate
        VTYPE           = PARAMS.engbertkliegl.VTYPE;           % Velocity types (2 = using moving average)
        SCREEN_RES      = PARAMS.engbertkliegl.SCREEN_RES;      % in pix
        SCREEN_SIZE     = PARAMS.engbertkliegl.SCREEN_SIZE;     % in mm (LNI)
        VIEWING_DIST    = PARAMS.engbertkliegl.VIEWING_DIST;    % in mm
        MININTERSACC    = PARAMS.engbertkliegl.MININTERSACC;
        MINDUR          = PARAMS.engbertkliegl.MINDUR;          % Minimum duration (number of samples)
        MINFIXDUR       = PARAMS.extra.MINFIXDUR;
    else
        VTHRES      = 6;                  % Velocity threshold
        SAMPLING    = 300; %eall.srate;    % Sampling rate
        VTYPE       = 2;                  % Velocity types (2 = using moving average)
        SCREEN_RES  = [1920,1080];         % in pix
        SCREEN_SIZE = [255*2,145*2];      % in mm (LNI)
        VIEWING_DIST= 600;                % in mm
        % SCREEN_SIZE=[376,292];      % in mm (neurolab)
        % VIEWING_DIST=600; %in mm
        MININTERSACC = 0.050; % in secs
        MINDUR      = 3;   % Minimum duration (number of samples)
        % MINDUR      = 3*SAMPLING/250;   % Minimum duration (samples) % Engbert uses 3s=6ms/12ms for SR=500/250Hzend
        MINFIXDUR   = 0.006; % IN MS, only for fixation / saccade detection
    end
    
    % Derived from parameters
    mmPP        = mean(SCREEN_SIZE ./ SCREEN_RES); %size of 1pix in mm
    DPP         = atand(mmPP / VIEWING_DIST) * 60;   %=2.13 Angular minutes per pixel
    
    % Time, xl, yl, time, xr, yr
    d           = eall.samples;  % THE SAMPLES, times are in secs
    
    
    
    
    %% Find blinks
    % ahora tengo que hacerlos distintos porque no tengo nans. Empiezo probando
    % con una saturacion negativa
    
    % mind = min(d) con el minimo no funciona bien... voy a ponerle los limites
    % de la pantalla
    
    ind = (d(:,2) < 0 | d(:,2) > SCREEN_RES(1));
    d(ind,2) = nan;
    ind = (d(:,3) < 0 | d(:,3) > SCREEN_RES(2));
    d(ind,3) = nan;
    if size(d, 2) > 4
        ind = (d(:,5) < 0 | d(:,5) > SCREEN_RES(1));
        d(ind,5) = nan;
        ind = (d(:,6) < 0 | d(:,6) > SCREEN_RES(2));
        d(ind,6) = nan;
    end
    
    % figure;
    %     subplot(2,1,1)
    %         hold on;
    %             plot(d(:,1),1024*isnan(d(:,2)),'r');
    %             plot(d(:,1),d(:,2),'b');
    %         hold off
    %         ylim([0 PARAMS.engbertkliegl.SCREEN_RES(1)])
    %     subplot(2,2,3)
    %         hold on;
    %             plot(d(isnan(d(:,2)),1),ones(sum(isnan(d(:,2))),1),'r.');
    %             plot(d(isnan(d(:,3)),1),0.5*ones(sum(isnan(d(:,3))),1),'g.');
    %         hold off
    %         ylim([0 1.5])
    %     subplot(2,2,4)
    %         hold on;
    %             plot(d(isnan(d(:,5)),1),ones(sum(isnan(d(:,5))),1),'r.');
    %             plot(d(isnan(d(:,6)),1),0.5*ones(sum(isnan(d(:,6))),1),'g.');
    %         hold off
    %         ylim([0 1.5])
    
    
    if (strcmpi(whicheyeanalyze,'BOTH'))
        index_nan   = find(isnan(d(:,2)) |  isnan(d(:,3)) | isnan(d(:,5)) | isnan(d(:,6)))';  % nanX = nanY = nan %find blinks
        %     index_nan   = find(any(isnan(d(:,[2 3 5 6]))'));
    elseif (strcmpi(whicheyeanalyze,'RIGHT'))
        %     index_nan   = find(isnan(d(:,2)));                  % nanX = nanY = nan %find blinks
        index_nan   = find(any(isnan(d(:,[5 6]))'));
        
    elseif (strcmpi(whicheyeanalyze,'LEFT'))
        index_nan   = find(any(isnan(d(:,[2 3]))'));
    end
    
    if isempty(index_nan)
        %     blink_end   = [];
        %     blink_start = []; %middle ones
        %     blink_start = []; %the extremes
        blink       = [];
        index_data  = (1:size(d,1))'; %find good data lines
    else
        blink_end   = d(index_nan(diff(index_nan) > 1) );
        blink_start = [d(index_nan(find(diff(index_nan) > 1) + 1))']'; %middle ones
        blink_start = cat(1, d(index_nan(1)), blink_start(1 : end - 1)')'; %the extremes
        blink       = [blink_start' blink_end'];
        index_data  = setdiff(1 : size(d, 1), index_nan)'; %find good data lines
    end
    % CAUTION: 'd' is continuous across samples and time, 'index_data' are the
    % samples without NaNs in both eyes (binocular data), or in one eye
    % (monocular data).
    
    
    %% Running the Kliegl Algorithm
    % Binocular data, analyzed as binocular, and reported as Right eye as
    % default
    if strcmpi(whicheyeanalyze,'BOTH')
        
        ix=5; iy=6; %will calculate binocular but report right eye fixations/saccades by default
        xl = DPP*d(index_data,2:3);
        xl(:,1) = xl(:,1) - mean(xl(:,1));
        xl(:,2) = xl(:,2) - mean(xl(:,2));
        xr = DPP*d(index_data,5:6);
        xr(:,1) = xr(:,1) - mean(xr(:,1));
        xr(:,2) = xr(:,2) - mean(xr(:,2));
        
        % Compute 2D velocity vectors
        vl = vecvel(xl,SAMPLING,VTYPE);
        vr = vecvel(xr,SAMPLING,VTYPE);
        
        % Detection of microsaccades
        sacl = microsacc(xl,vl,VTHRES,MINDUR);
        sacr = microsacc(xr,vr,VTHRES,MINDUR);
        
        % Testing for binocular saccades via temporal overlap
        [sac monol monor] = binsacc(sacl,sacr);
        
        
        
        if strcmpi(whicheye,'LEFT')
            ix=2;iy=3; %change monocular data to be reported
            temp(:,1:15)    = sac(:,16:30);
            temp(:,16:30)   = sac(:,1:15);
            sac             = temp;
        end
        isbinocular.eall = (abs(diff(sac(:,[1 16])')')>10); % I never used this and it usually cause troubles
        
        % Binocular data, analyzed as monocular (Left), and reported as Left eye
    elseif strcmpi(whicheyeanalyze,'LEFT') % binocular
        
        ix = 2; iy = 3;
        xl = DPP * d(index_data, 2:3);
        xl(:,1) = xl(:,1) - mean(xl(:, 1));
        xl(:,2) = xl(:,2) - mean(xl(:, 2));
        
        % Compute 2D velocity vectors
        vl = vecvel(xl, SAMPLING, VTYPE);
        
        % Detection of microsaccades
        sacl = microsacc(xl, vl, VTHRES, MINDUR);
        sacr = [];
        sac = sacl;
        isbinocular.eall = ones(size(sac, 1) - 1, 1); % I never used this and it usually cause troubles, here I have to create a vector
        
        % Binocular data, analyzed as monocular (Right), and reported as Right eye
    elseif strcmpi(whicheyeanalyze,'RIGHT') % binocular
        
        ix = 5; iy = 6; % THIS IS ONLY FOR THE DATA IN LEICESTER!
        xr = DPP * d(index_data, 5:6);
        xr(:,1) = xr(:,1) - mean(xr(:,1));
        xr(:,2) = xr(:,2) - mean(xr(:,2));
        
        % Compute 2D velocity vectors
        vr = vecvel(xr, SAMPLING, VTYPE);
        
        % Detection of microsaccades
        sacr = microsacc(xr, vr, VTHRES, MINDUR);
        sacl = [];
        
        % If a chunk of trial presents no saccades lets just return the empty
        % structure
        if isempty(sacr)
            
            saccade.description_11col = 'onset end dur x_ini y_ini x_end y_end amp(deg) peak_V(deg/s) angle(deg) otheramp(deg)';
            saccade.allsac = [];
            saccade.macrosac = [];
            saccade.microsac = [];
            saccade.macrofix = [];
            saccade.macrofixprops = [];
            saccade.fix     = [];
            saccade.fixprops= [];
            
            blink = [];
            return
        else
            
            sac = sacr;
            isbinocular.eall = ones(size(sac,1) - 1, 1); % I never used this and it usually cause troubles, here I have to create a vector
        end
    else
        disp('Problem identifying tracker eye or it is monocular')
    end
    
    % CAUTION: 'd' is continuous across samples and time, 'index_data' are the
    % samples without NaNs in both eyes (binocular data), or in one eye
    % (monocular data).
    % AND NOW... 'sac' runs in samples of clean data (i.e. d(index_data,:)) so,
    % if I want to recover the original data later, I have to write something
    % like: 'd(index_data(sac(i)),:)'
    
    
    %% Remove artificial saccades right before each blink
    
    ssac = sac(:,1); %saccade onset
    esac = sac(:,2); %saccade end
    
    % ind_artsac = [];
    ind_artsac = zeros(1, size(sac,1) * size(blink,1)); % preallocate memory
    count = 1;
    for ii = 1 : size(sac, 1)
        for jj = 1 : size(blink, 1)
            
            if ~isempty(blink) && size(blink,2) == 2
                if  (blink(jj, 1) > d(index_data(ssac(ii), 1)) && ...
                        blink(jj, 2) < d(index_data(esac(ii), 1)))
                    %             ind_artsac = [ind_artsac ii];
                    ind_artsac(count) = ii;
                    count = count + 1;
                end
            end
        end
        if mod(ii, 100) == 0
            disp(ii)
        end
    end
    
    ind_artsac = ind_artsac(ind_artsac > 0);
    
    sac(ind_artsac,:) = [];
    
    
    %% Saccade main sequence (peak velocity vs magnitude)
    if doplot
        [R,PVAL] = corr(sac(:,13),sac(:,3));
        figure;
        if ~ischar(eall.resacc); subplot(2,1,1); end
        %loglog(sac(:,13),sac(:,3),'r.','linewidth',1);
        plot(sac(:,13),sac(:,3),'r.','linewidth',1);
        h_legend=legend('offline detection');
        set(h_legend,'FontSize',16);
        legend('boxoff');
        xlabel('Magnitude [deg]');
        ylabel('Peak Velocity [deg/s]');
        %set(gca,'YLim',[5 5000],'XLim',[0.02 100]);
        set(gca,'YLim',[5 200],'XLim',[0.02 1]);
        grid on
        %         hold on
        if ~ischar(eall.resacc); subplot(2,1,2);
            %loglog(eall.lesacc(:,8),eall.lesacc(:,9),'g.','linewidth',1);
            plot(eall.resacc(:,8),eall.resacc(:,9),'g.','linewidth',1);
            h_legend=legend('online detection');
            %h_legend=legend('clean saccades');
            set(h_legend,'FontSize',16);
            legend('boxoff');
            xlabel('Magnitude [deg]');
            ylabel('Peak Velocity [deg/s]');
            grid on
            %         hold off
            %set(gca,'YLim',[5 5000],'XLim',[0.02 100]);
            set(gca,'YLim',[5 200],'XLim',[0.02 1]);
        end
        outfile = sprintf('%s%s_%s', DIR.figuresOUT, subjid,'mainsequence_compare_test_linear');
        print(gcf,'-dpng',outfile);
    end
    
    %% Separate microsaccades and eall saccades and write an array with the
    %% first nine columns compatible with eall.lefix
    
    if isempty(sac)
        
        saccade.description_11col = 'onset end dur x_ini y_ini x_end y_end amp(deg) peak_V(deg/s) angle(deg) otheramp(deg)';
        saccade.allsac = [];
        saccade.macrosac = [];
        saccade.microsac = [];
        saccade.macrofix = [];
        saccade.macrofixprops = [];
        saccade.fix     = [];
        saccade.fixprops= [];
        
        blink = [];
        
        disp('here')
        return
    end
    
    
    allsac = zeros(size(sac, 1), 9);
    for s = 1 : size(sac, 1)
        idx = sac(s, 1) : sac(s, 2);
        allsac(s,1) = d(index_data(idx(1)), 1);    %onset in secs
        allsac(s,2) = d(index_data(idx(end)), 1);  %end in secs
        allsac(s,4) = d(index_data(idx(1)), ix);    %x ini sacc
        allsac(s,5) = d(index_data(idx(1)), iy);    %y ini sacc
        allsac(s,6) = d(index_data(idx(end)), ix);  %x end sacc
        allsac(s,7) = d(index_data(idx(end)), iy);  %y end sacc
    end
    allsac(:,3) = sac(:, 15)  / SAMPLING;        %duration in secs
    allsac(:,8) = sac(:, 13);                      %amplitude (in deg)
    allsac(:,9) = sac(:, 3);                       %peak velocity (in deg/sec)
    
    %the following are new fields
    allsac(:, 10) = sac(:, 14);                     %saccade angle
    allsac(:, 11) = sac(:, 12) / (60 * SAMPLING);       %alternative saccade amplitude (in deg)
    
    %foo=-allsac(1:end-1,2)+allsac(2:end,1);
    %min(foo)
    
    %% Correlation between saccade amplitude (Kliegl, based on displacements in
    %% x, y) and alternative saccade amplitude (estimated from average velocity)
    if doplot
        figure;
        plot(allsac(:, 8), allsac(:, 11), 'b.')
        %h_legend=legend('online detection');
        %set(h_legend,'FontSize',12);
        %set(h_legend,'Location','NorthEast');
        %xlim([d(index_d_start) d(index_d_end)]);
        xlabel('Amplitude [deg]');
        ylabel('Velocity-based amplitude [deg]');
        legend('boxoff');
        hold off;
        outfile = sprintf('%s%s_%s', DIR.figuresOUT, subjid, 'correl_ampltitude');
        print(gcf, '-dpng', outfile);
    end
    
    %% Detect (separate) microsaccades and macrosaccades
    % We separate micro and macro saccades as in Otero-Millan et al.,etc with a
    % fixed threshold in the amplitude of 1 degree,
    % We also add a Minimum Intersaccadic Interval, this was very variable
    % across bibliography (see Matias presentation)
    % - 50ms as in Dimigen, Kliegl, 2009;
    % - or 20ms as in Otero-Millan et al.,etc
    
    f1 = ( allsac(1:end,8) < 1 );
    f2 = [1; ((allsac(2:end, 1) - allsac(1:end-1, 2)) > MININTERSACC)];
    
    % figure;
    % yh = allsac(2:end,1)-allsac(1:end-1,2);
    % [yh xh] = hist(yh(f1(1:end-1)),0:0.01:2);
    % hold on
    % plot(xh,yh,'k-')
    % plot([MININTERSACC MININTERSACC],[0 max(yh)],'r-')
    % hold off
    % xlim([0 2])
    
    index_microsac = find(f1 & f2);       % identify microsaccades
    if (allsac(end,8) < 1);           % Add last sacc, if it was a microsacc
        index_microsac = [index_microsac; size(allsac, 1)];
    end
    microsac = allsac(index_microsac, :);  % select microsaccades
    if strcmpi(whicheyeanalyze, 'BOTH');
        isbinocular.microsac = isbinocular.eall(index_microsac);
    end
    
    %index_macrosac=find(~f1 & f2);      % identify macrosaccades (this has
    %the problem of leaving out some macrosaccades which are in between short
    %microsaccades
    index_macrosac = find(~f1);      % identify macrosaccades
    if (allsac(end,8) >= 1);          % Add last sacc, if it was a macrosacc
        index_macrosac = [index_macrosac; size(allsac, 1)];
    end
    macrosac = allsac(index_macrosac, :);  % select macrosaccades
    if strcmpi(whicheyeanalyze, 'BOTH');
        isbinocular.macrosac = isbinocular.eall(index_macrosac);
    end
    
    % CAUTION: 'index_microsac' and 'index_macrosac' are very useful later, do
    % not rewrite it!!
    
    %% Store output in a new structure "saccade"
    % Kliegl Saccades
    saccade.description_11col   = 'onset end dur x_ini y_ini x_end y_end amp(deg) peak_V(deg/s) angle(deg) otheramp(deg)';
    saccade.allsac              = allsac;
    saccade.macrosac            = macrosac;
    saccade.microsac            = microsac;
    if strcmpi(whicheyeanalyze, 'BOTH')
        saccade.isbinocular     = isbinocular;
        saccade.monol           = monol;
        saccade.monor           = monor;
    end
    
    % EyeLink Saccades for comparisons
    % saccade.description_9col    = 'onset end dur x_ini y_ini x_end y_end amp(deg) peak_V(deg/s)';
    % saccade.lesacc              = eall.lesacc;
    % saccade.resacc              = eall.resacc;
    
    
    %% Fixations from Kliegl saccades
    % Macrofixations
    for ii = 1 : size(sac, 1) - 1
        if ismember(ii, index_macrosac);
            if (ii < max(index_macrosac))     % I have to check that it is not
                % the last saccade because I use
                % the nth and (n+1)th saccade to
                % the calculation of the fixation
                icurr   = ii;
                inext   = index_macrosac(find(index_macrosac==ii)+1);
                sini    = sac(icurr,2);     % real samples (?)
                send    = sac(inext,1);     % real samples (?)
                idx     = sini+1:send-1;    % limits (in samples)
                % Main properties as in the EyeLink output...
                macrofix(ii,1)  = allsac(icurr,2) + 1 / SAMPLING;                 % START FIXATION (in secs) (agrego un sample)
                macrofix(ii,2)  = allsac(inext,1) - 1 / SAMPLING;                 % END FIXATION (in secs) (quito un sample)
                macrofix(ii,3)  = macrofix(ii,2) - macrofix(ii,1) + 1 / SAMPLING; % DURATION (in secs)
                macrofix(ii,4)  = nanmean(d(index_data(idx), ix));                % x (in pixs)
                macrofix(ii,5)  = nanmean(d(index_data(idx), iy));                % y (in pixs)
                macrofix(ii,6)  = nanmean(d(index_data(idx), 4));                 % pupil
                
                % Some extra properties
                macrofixprops(ii,1) = nanstd(d(index_data(idx),ix)); % stdx
                macrofixprops(ii,2) = nanstd(d(index_data(idx),iy)); % stdy
                macrofixprops(ii,3) = nanmin(d(index_data(idx),ix)); % minx
                macrofixprops(ii,4) = nanmin(d(index_data(idx),iy)); % miny
                macrofixprops(ii,5) = nanmax(d(index_data(idx),ix)); % maxx
                macrofixprops(ii,6) = nanmax(d(index_data(idx),iy)); % maxy
                
                %             % Some extra extra properties quantiles
                %                 tempx = quantile(d(index_data(idx),2),[0.025 0.25 0.5 0.75 0.975]);
                %                 tempy = quantile(d(index_data(idx),3),[0.025 0.25 0.5 0.75 0.975]);
                %                 macrofixprops2(ii,:) = [tempx tempy];
                
            else
                
                % I didn't solve this, if the last fixation of the experiment turns
                % important I'll do it.
                % Last fixation added: Matias, 13/4/2010
                icurr   = ii;
                sini    = sac(icurr,2);     % real samples (?)
                % Main properties as in the EyeLink output...
                macrofix(ii,1)  = allsac(icurr,2) + 1 / SAMPLING;                 % START FIXATION (in secs)
                macrofix(ii,2)  = d(index_data(end), 1);                           % END FIXATION (in secs)
                macrofix(ii,3)  = macrofix(ii,2) - macrofix(ii,1) + 1 / SAMPLING;   % DURATION (in secs)
                macrofix(ii,4)  = nanmean(d(index_data(sini+1:end), ix));          % x (in pixs)
                macrofix(ii,5)  = nanmean(d(index_data(sini+1:end), iy));          % y (in pixs)
                macrofix(ii,6)  = nanmean(d(index_data(sini+1:end), 4));           % pupil
                % Some extra properties
                macrofixprops(ii,1) = nanstd(d(index_data(sini+1:end), ix)); % stdx
                macrofixprops(ii,2) = nanstd(d(index_data(sini+1:end), iy)); % stdy
                macrofixprops(ii,3) = nanmin(d(index_data(sini+1:end), ix)); % minx
                macrofixprops(ii,4) = nanmin(d(index_data(sini+1:end), iy)); % miny
                macrofixprops(ii,5) = nanmax(d(index_data(sini+1:end), ix)); % maxx
                macrofixprops(ii,6) = nanmax(d(index_data(sini+1:end), iy)); % maxy
            end
        end
    end
    
    % Filtering macrofixations,
    % i.e. exclude fixations with blinks or fixations shorter than 6ms
    ishortmacrofix          = zeros(size(allsac, 1) - 1, 1); % 1=short fixation, 0=good fixation
    isblinkbetmacrosacc     = zeros(size(allsac, 1) -1, 1); % 1=bad fixation due to blink between saccades, 0=good fixation
    for ii = 1 : size(allsac, 1) - 1
        if ismember(ii, index_macrosac);
            if ii < max(index_macrosac) && size(blink,2) == 2
                tini = allsac(ii, 2);
                tend = allsac(index_macrosac(find(index_macrosac == ii) + 1), 1);
                for jj=1:size(blink, 1)
                    if (blink(jj, 1) > tini && blink(jj, 2) < tend)
                        isblinkbetmacrosacc(ii) = 1;
                    end
                end
                if tend < tini + MINFIXDUR
                    ishortmacrofix(ii) = 1;
                end
            else
            end
        end
    end
    
    ismacrosac              = zeros(size(allsac,1) - 1, 1);
    
    %ismacrosac = zeros(size(macrofix,1),1);
    %   ismacrosac(index_macrosac(index_macrosac<size(allsac,1)))=1; %3/3/10
    %   this line was commented...why? the next line can give an error if there
    %   are two macrosaccades in the very last one...so I comment it
    %ismacrosac_foo=zeros(size(allsac,1)-1,1);
    %ismacrosac_foo(index_macrosac(1:end-1,1))=1;
    %ismacrosac(index_macrosac(1:end-1,1))=1; %%%identical except for the
    %last element!!!
    ismacrosac(index_macrosac(index_macrosac<size(allsac, 1))) = 1;
    %ismacrosac(index_macrosac(index_macrosac<size(macrofix,1)))=1;
    
    
    
    if isempty(ismacrosac) || isempty(ishortmacrofix) || isempty(isblinkbetmacrosacc) || ~exist('macrofix', 'var')
        
        macrofix = [];
        macrofixprops = [];
    else
        
        macrofix                = macrofix(ismacrosac & ~ishortmacrofix & ~isblinkbetmacrosacc,:);
        macrofixprops           = macrofixprops(ismacrosac & ~ishortmacrofix & ~isblinkbetmacrosacc,:);
        % macrofixprops2          = macrofixprops2(ismacrosac & ~ishortmacrofix & ~isblinkbetmacrosacc,:);
    end
    
    % and saving...
    saccade.macrofix        = macrofix;
    saccade.macrofixprops   = macrofixprops;
    % saccade.macrofixprops2  = macrofixprops2;
    
    %CAUTION: consecutive saccades lead to spurious short fixations.
    %with eyelink, the shortest fixations are 6ms long.
    %I introduced the same filter here.
    
    %% recalculate fixations from Kliegl saccades
    
    if size(allsac) > 1
        
        for ii = 1 : size(allsac) - 1
            idx = sac(ii,2) + 1 : sac(ii + 1, 1) - 1;                % limits (in samples)
            % Main properties as in the EyeLink output...
            fix(ii,1)   = allsac(ii, 2) + 1 / SAMPLING;           % START FIXATION (in secs)
            fix(ii,2)   = allsac(ii + 1, 1) - 1 / SAMPLING;         % END FIXATION (in secs)
            fix(ii,3)   = fix(ii, 2) - fix(ii, 1) + 1 / SAMPLING;    % DURATION
            fix(ii,4)   = nanmean(d(index_data(idx), ix));
            fix(ii,5)   = nanmean(d(index_data(idx), iy));
            fix(ii,6)   = nanmean(d(index_data(idx), 4));
            
            % Some extra properties
            fixprops(ii,1)  = nanstd(d(index_data(idx), ix)); % stdx
            fixprops(ii,2)  = nanstd(d(index_data(idx), iy)); % stdy
            fixprops(ii,3)  = nanmin(d(index_data(idx), ix)); % minx
            fixprops(ii,4)  = nanmin(d(index_data(idx), iy)); % miny
            fixprops(ii,5)  = nanmax(d(index_data(idx), ix)); % maxx
            fixprops(ii,6)  = nanmax(d(index_data(idx), iy)); % maxy
        end
        
        % Filtering fixations,
        % i.e. exclude fixations with blinks or fixations shorter than 6ms
        
        ishortfix       = ( fix(:,3) < MINFIXDUR );
        
        isblinkbetsacc  = zeros(size(allsac,1) - 1 , 1); %0=bad fixation due to blink between saccades, 1=good fixation
        
        if size(blink, 2) ==2
            for ii = 1 : size(allsac, 1) - 1
                for jj = 1 : size(blink, 1)
                    if (blink(jj, 1) > allsac(ii, 2) && blink(jj, 2) < allsac(ii + 1, 1))
                        isblinkbetsacc(ii) = 1;
                    end
                end
            end
        end
        
        fix             = fix(~ishortfix & ~isblinkbetsacc, :);
        fixprops        = fixprops(~ishortfix & ~isblinkbetsacc, :);
        
        % and saving...
        saccade.fix     = fix;
        saccade.fixprops= fixprops;
        
    else
        saccade.fix     = [];
        saccade.fixprops= [];
    end
    
    %CAUTION: consecutive saccades lead to spurious short fixations.
    %with eyelink, the shortest fixations are 6ms long.
    %I introduced the same filter here.
    
    
catch ME1
    rethrow(ME1)
end




% %% Plot raw data and saccades for a short time window
% if doplot
%     nsac_to_plot    = 15; %number of saccades to plot
%     initial_saccade = 1;
%     %window_saccade_start=index_data(sac(initial_saccade,1),1); %sample of the first saccade to plot
%     %window_saccade_end=index_data(sac(nsac_to_plot+initial_saccade,2),1); %sample of the last saccade to plot
%     %tminsacc=min(find(allsac(initial_saccade,1)>d(:,1)));
%     %tmaxsacc=max(find(eall.lesacc(:,1)<d(window_saccade_end,1)));
%
%     window_saccade_start= allsac(initial_saccade,1);                % initial time of first saccade of interest (in ms)
%     window_saccade_end  = allsac(initial_saccade+nsac_to_plot,1);   % final time of last saccade of interest (in ms)
%
%     index_d_start       = min(find(d(:,1)>window_saccade_start)); %lower limit (in samples)
%     index_d_end         = max(find(d(:,1)<window_saccade_end)); %upper limit (in samples)
%
%     figure;
%     if ~ischar(eall.resacc); subplot(3,1,1) %raw data
%     else                    subplot(2,1,1)
%     end
%     index_window=index_d_start:index_d_end;
%     plot(d(index_window,1),d(index_window,ix),'b-');
%     h_legend=legend('raw data');
%     set(h_legend,'FontSize',12);
%     legend('boxoff');
%     xlim([d(index_d_start) d(index_d_end)]);
%
%     if ~ischar(eall.resacc); subplot(3,1,2) %kliegl saccades
%     else                    subplot(2,1,2)
%     end
%     hold on;
%     for s=initial_saccade:nsac_to_plot+initial_saccade
%         s_window_saccade_start=allsac(s,1); %initial time of first saccade of interest (in ms)
%         s_window_saccade_end=allsac(s,2); %final time of last saccade of interest (in ms)
%         s_index_d_start=min(find(d(:,1)>s_window_saccade_start)); %lower limit (in samples)
%         s_index_d_end=max(find(d(:,1)<s_window_saccade_end)); %upper limit (in samples)
%         s_index_window=s_index_d_start:s_index_d_end;
%         plot(d(s_index_window,1),d(s_index_window,ix),'r-','linewidth',1.5);
%     end
%     h_legend=legend('offline detection');
%     set(h_legend,'FontSize',12);
%     legend('boxoff');
%     xlim([d(index_d_start) d(index_d_end)]);
%     ylim([0 PARAMS.engbertkliegl.SCREEN_RES(1)]);
%     hold off;
%
%     if ~ischar(eall.resacc); subplot(3,1,3) %clean eyelink saccades
%         hold on;
%         sampleminsacc=min(find(eall.resacc(:,1)>window_saccade_start));
%         samplemaxsacc=max(find(eall.resacc(:,1)<window_saccade_end));
%         for s=sampleminsacc:samplemaxsacc;
%             index_inisacc_sample=find(eall.samples(:,1)>=eall.resacc(s,1),1,'first');
%             index_endsacc_sample=find(eall.samples(:,1)>=eall.resacc(s,2),1,'first');
%             idx = index_inisacc_sample:index_endsacc_sample;
%             plot(eall.samples(idx,1),eall.samples(idx,2),'g-','linewidth',1.5);
%         end
%         h_legend=legend('online detection');
%         set(h_legend,'FontSize',12);
%         set(h_legend,'Location','NorthEast');
%         ylim([0 PARAMS.engbertkliegl.SCREEN_RES(1)]);
%         xlim([d(index_d_start) d(index_d_end)]);
%         xlabel('Time [ms]');
%         ylabel('X-position [pix]');
%         legend('boxoff');
%         hold off;
%         outfile = sprintf('%s%s_%s', DIR.figuresOUT, subjid, 'x_timewindow');
%         print(gcf,'-dpng',outfile);
%     end
%     figure;
%     if ~ischar(eall.resacc); subplot(3,1,1) % raw data
%     else                    subplot(2,1,1)
%     end
%     index_window=index_d_start:index_d_end;
%     plot(d(index_window,1),d(index_window,iy),'b-');
%     h_legend=legend('raw data');
%     set(h_legend,'FontSize',12);
%     legend('boxoff');
%     xlim([d(index_d_start) d(index_d_end)]);
%
%     if ~ischar(eall.resacc); subplot(3,1,2) % kliegl saccades
%     else                    subplot(2,1,2)
%     end
%     hold on;
%     for s=initial_saccade:nsac_to_plot+initial_saccade
%         s_window_saccade_start=allsac(s,1); %initial time of first saccade of interest (in ms)
%         s_window_saccade_end=allsac(s,2); %final time of last saccade of interest (in ms)
%         s_index_d_start=min(find(d(:,1)>s_window_saccade_start)); %lower limit (in samples)
%         s_index_d_end=max(find(d(:,1)<s_window_saccade_end)); %upper limit (in samples)
%         s_index_window=s_index_d_start:s_index_d_end;
%         plot(d(s_index_window,1),d(s_index_window,iy),'r-','linewidth',1.5);
%     end
%     h_legend=legend('offline detection');
%     set(h_legend,'FontSize',12);
%     legend('boxoff');
%     ylim([0 PARAMS.engbertkliegl.SCREEN_RES(2)]);
%     xlim([d(index_d_start) d(index_d_end)]);
%     hold off;
%
%     if ~ischar(eall.resacc); subplot(3,1,3)  % clean eyelink saccades
%         hold on;
%         sampleminsacc=min(find(eall.resacc(:,1)>window_saccade_start));
%         samplemaxsacc=max(find(eall.resacc(:,1)<window_saccade_end));
%         for s=sampleminsacc:samplemaxsacc;
%             index_inisacc_sample=find(eall.samples(:,1)>=eall.resacc(s,1),1,'first');
%             index_endsacc_sample=find(eall.samples(:,1)>=eall.resacc(s,2),1,'first');
%             idx = index_inisacc_sample:index_endsacc_sample;
%             plot(eall.samples(idx,1),eall.samples(idx,3),'g-','linewidth',1.5);
%         end
%         h_legend=legend('online detection');
%         set(h_legend,'FontSize',12);
%         set(h_legend,'Location','NorthEast');
%         ylim([0 PARAMS.engbertkliegl.SCREEN_RES(2)]);
%         xlim([d(index_d_start) d(index_d_end)]);
%         xlabel('Time [ms]');
%         ylabel('Y-position [pix]');
%         legend('boxoff');
%         hold off;
%     end
%     outfile = sprintf('%s%s_%s', DIR.figuresOUT, subjid, 'y_timewindow');
%     print(gcf,'-dpng',outfile);
% end
%
%
% %% saccade main sequence (peak velocity vs magnitude)
% if doplot
%     figure;
%     if ~ischar(eall.resacc); subplot(2,1,1); end
%     loglog(macrosac(:,8),macrosac(:,9),'r.','linewidth',1);
%     %loglog(eall.lesac(:,8),eall.lesac(:,9),'r.','linewidth',1);
%     h_legend=legend('offline detection');
%     set(h_legend,'Location','NorthWest');
%     set(h_legend,'FontSize',16);
%     legend('boxoff');
%     xlabel('Magnitude [deg]');
%     ylabel('Peak Velocity [deg/s]');
%     set(gca,'YLim',[5 5000],'XLim',[0.02 100]);
%     grid on
%     hold on
%     if ~ischar(eall.resacc); subplot(2,1,2);
%         loglog(eall.resacc(:,8),eall.resacc(:,9),'g.','linewidth',1);
%         h_legend=legend('online detection');
%         set(h_legend,'Location','NorthWest');
%         %h_legend=legend('clean saccades');
%         set(h_legend,'FontSize',16);
%         legend('boxoff');
%         xlabel('Magnitude [deg]');
%         ylabel('Peak Velocity [deg/s]');
%         grid on
%         hold off
%         set(gca,'YLim',[5 5000],'XLim',[0.02 100]);
%     end
%     outfile = sprintf('%s%s_%s', DIR.figuresOUT, subjid, 'mainsequence_compare');
%     print(gcf,'-dpng',outfile);
% end
%
% %% Microsaccade main sequence (peak velocity vs magnitude)
% if doplot
%     figure;
%     subplot(2,2,1)
%     plot(microsac(:,8),microsac(:,9),'g.');
%     %h_legend=legend('eall saccades');
%     %set(h_legend,'FontSize',16);
%     %legend('boxoff');
%     xlabel('Magnitude [deg]');
%     ylabel('Peak Velocity [deg/s]');
%     set(gca,'YLim',[0 100],'XLim',[0 1]);
%     hold on
%     subplot(2,2,2)
%     hist(microsac(:,9),50,'g');
%     h = findobj(gca,'Type','patch');
%     set(h,'FaceColor','g','EdgeColor','w')
%     xlabel('Peak Velocity [deg/s]');
%     ylabel('Number of Microsaccades');
%     hold off
%     set(gca,'XLim',[0 100]);
%     hold on
%     subplot(2,2,3)
%     hist(microsac(:,8),50);
%     h = findobj(gca,'Type','patch');
%     set(h,'FaceColor','g','EdgeColor','w')
%     xlabel('Magnitude [deg]');
%     ylabel('Number of Microsaccades');
%     hold off
%     set(gca,'XLim',[0 1]);
%     hold on
%     subplot(2,2,4)
%     hist(1000*microsac(:,3),15);
%     h = findobj(gca,'Type','patch');
%     set(h,'FaceColor','g','EdgeColor','w')
%     xlabel('Duration [ms]');
%     ylabel('Number of Microsaccades');
%     hold off
%     set(gca,'XLim',[0 40]);
%     outfile = sprintf('%s%s_%s', DIR.figuresOUT, subjid,'mainsequence_microsac');
%     print(gcf,'-dpng',outfile);
% end
%
% %%
% if doplot
%     T{1} = [allsac(:,10)    allsac(:,8)];
%     T{2} = [macrosac(:,10)  macrosac(:,8)];
%     T{3} = [microsac(:,10)  microsac(:,8)];
%
%     polarplot(T,'rose');
%     outfile = sprintf('%s%s_%s', DIR.figuresOUT, subjid, 'angular_distribution');
%     print(gcf,'-dpng',outfile);
%
%     polarplot(T,'polar');
%     outfile = sprintf('%s%s_%s', DIR.figuresOUT, subjid, 'spatial_distribution');
%     print(gcf,'-dpng',outfile);
% end



%% plots from an old version
% figure
% %plot(xl(:,1),xl(:,2),'b-');
% set(gca,'FontSize',16);
% axis square;
% hold on;
% for s=1:size(sac,1)
%     idx = sac(s,1):sac(s,2);
%     plot(xl(idx,1),xl(idx,2),'r.','linewidth',2);
% end
% hold off;
% xlim([-30 30]);
% ylim([-30 30]);
% set(gca,'xtick',[-30 -15 0 15 30]);
% set(gca,'ytick',[-30 -15 0 15 30]);
% xlabel('x [min arc]');
% ylabel('y [min arc]');
%%

%
% eyelink_intsacc=SAMPLING*(eall.lesac(:,2)-eall.lesac(:,1))/1000;
% macrosac_intsacc=macrosac(:,2)-macrosac(:,1);
%
%
% figure;
% plot(xl(:,1),xl(:,2),'b-');
% set(gca,'FontSize',16);
% axis square;
% hold on;
% for s=1:size(sac,1)
%     idx = sac(s,1):sac(s,2);
%     plot(xl(idx,1),xl(idx,2),'r.-','linewidth',2);
% end
% hold off;
%xlim([-30 30]);
%ylim([-30 30]);
%set(gca,'xtick',[-30 -15 0 15 30]);
%set(gca,'ytick',[-30 -15 0 15 30]);
%xlabel('x [min arc]');
%ylabel('y [min arc]');
%






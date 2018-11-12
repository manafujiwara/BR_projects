function label_IM(DIR,Cond)

%%% here we make a variable "BPinfo", and fix index of epc %%%%%%%%%%%%%%%%

time2see = [0, 1.5]; % in sec from stimulus onset
tp2see = ( time2see + 1 ) * 300;

load( [DIR.allSub '/IM_eyeAllSub.mat'] ,'iOKN', 'vel', 'slwPhsBox', 'slwPhsMed','subID')
load( [DIR.allSub '/IM_respAllSub.mat'] , 'lResp', 'rResp', 'resp', 'BPBox', 'lBPBox', 'rBPBox', 'BPMed', 'lBPMed', 'rBPMed', 'subID')

%% UNTIL CALCULATION OF BOXCAR WILL BE DONE
% tmprBoxcar = boxcar(rResp(:,7:906),30);
% tmplBoxcar = boxcar(lResp(:,7:906),30);
% respBoxcar = boxcar(resp(:,7:906),30);
% 
% rBPBox = [rResp(:,1:6),tmprBoxcar];
% lBPBox = [lResp(:,1:6),tmplBoxcar];
% BPBox  = [resp(:,1:6),respBoxcar];

%%
% --------------------------------------------------
% rResp/lResp: 1)subID, 2)irun, 3)iblk, 4)iepc,    |
% 5)stimLabel, 6)subjestState, 7:906)responseData  |
% --------------------------------------------------

%%% correct number of epoch %%%
% rResp(:,4) = rResp(:, 4)+1;
% lResp(:,4) = lResp(:, 4)+1;

%% CAUTION !!!!! %%%%%%%%%%%%%%%%%%%%%%%
% This needs to be removed once you run previous processes. The issue is
% already fixed in the codes, but here for now, to avoid wasting time, I
% just add this lines to fix the bug.
% rResp(:, 7:906) = [nan(1,900); rResp(1:end-1, 7:906)];
% lResp(:, 7:906) = [nan(1,900); lResp(1:end-1, 7:906)];

%% --- Button Press -------------------------------------------------------
% pull actuall report part only
BPt = resp( :, tp2see(1) + 1 + 6: tp2see(2) + 6 ); % button press (t)
BPt ( abs( BPt(:,1) ) == 1, : ) = nan ; % remove trials the subject kept pressing the button over the blank

%--- label ---------------------------------------------------------------
pBP = nansum(BPt,2)/ diff(tp2see); % Button Prss Proportion
label( pBP > 0 ) = 1;
label( pBP < 0 ) = -1;
label( isnan (pBP) | pBP == 0 ) = nan;

%-- RT (First Button Press & First CONSISTENT Button Press) ---------------
for itri = 1: size(BPt, 1)
    if isnan(label(itri))
        fBP(itri) = nan;
        fcBP(itri) = nan;
    else
        fBP(itri) = find ( BPt (itri, :), 1, 'first' );
        fcBP(itri) = find( BPt (itri, :) == label(itri), 1, 'first' );
    end
end

BPInfo   = [ BPBox(:, 1:6), pBP, label', fBP', fcBP', BPBox(:,7:906) ];
respInfo = [ BPBox(:, 1:6), pBP, label', fBP', fcBP', resp(:,7:906) ];

% -------------------------------------------------------------------------
% BPinfo: 1)subID, 2)irun, 3)iblk, 4)iepc, 5)stimLabel,                   |
% 6)subjestState, 7)BPP (-1~1), 8)label, 9)fBP, 10)fcBP, 11:910)BP(t)     |
% -------------------------------------------------------------------------

% quick check %%%%
rBP = BPInfo(BPInfo(:,8)==1,11:910);
lBP = BPInfo(BPInfo(:,8)==-1,11:910);
figure(2), clf
plot(mean(rBP,1), 'r')
hold on
plot(mean(lBP,1), 'b')
%%%%%%%%%%%%%%%%%%

save([DIR.allSub '/IM_BPinfo.mat'], 'BPInfo', 'respInfo')



%% --- eye data -----------------------------------------------------------
% ---------------------------------------------------------------
% 1)subID, 2)irun, 3)iblk, 4)iepc, 5)stimLabel, 6)subjestState  |
% 7)BPP (-1~1), 8)label, 9)fBP, 10)fcBP, 11:910)eyeinfo(t)      |
% ---------------------------------------------------------------

slwPhsBoxInfo = [ BPInfo(:,1:10),  slwPhsBox(:,7:906) ];
slwPhsMedInfo = [ BPInfo(:,1:10),  slwPhsMed(:,7:906) ];
iOKNInfo      = [ BPInfo(:,1:10),  iOKN(:,7:906) ];
velInfo       = [ BPInfo(:,1:10),  vel(:,7:906) ];

save([DIR.allSub '/' 'IM_eyeinfo.mat'], 'slwPhsBoxInfo', 'slwPhsMedInfo', 'iOKNInfo', 'velInfo')

end
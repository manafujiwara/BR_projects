
nState        = 2;

load([DIR.cfg 'UPDRS.mat'])

rejSubID = [ 105, 110, 118, 120, 121, 136, 137, 151];
tmp = score;
subID_1     = unique( tmp( tmp(:,2) == 1, 1 ) );
subID_2     = unique( tmp( tmp(:,2) == 2, 1 ) );
subID_3     = unique( tmp( tmp(:,2) == 3, 1 ) );
subID_4     = unique( tmp( tmp(:,2) == 4, 1 ) );
subID_5     = unique( tmp( tmp(:,2) == 5, 1 ) );
subID_com23 = subID_2( ismember(subID_2, subID_3) );
subID_com45 = subID_4( ismember(subID_4, subID_5) );

% subID_all = [subID_1; subID_com23; subID_com45 ];
% subID = subID_all( ~ismember(subID_all, rejSubID) );

% Set subject by hand for now
subID_valNC = [ 106, 108, 109, 111, 122, 126, 128, 141, 143, 146, 147];
subID_valPD = [ 102, 114, 115, 138, 139, 107, 116, 130, 144, 145, 150];
subID = [subID_valNC, subID_valPD];

for i = 1:size(UPDRS,2)
    switch UPDRS(i).state
        case 'nc'
            state= 1; %=> block level
        case 'mon'
            state = 2;
        case 'mof'
            state = 3;
        case 'don'
            state = 4;
        case 'dof'
            state = 5;
    end
    
    score(i,:) = [UPDRS(i).subID, state, UPDRS(i).totalScore];
end

for isub = 1:length(subID)
    
    if ismember(subID(isub), subID_valNC)
        state = 1;
    elseif ismember(subID(isub), subID_valPD)
        state = 2;
    end
    mscoreSub(isub, : ) = [ subID(isub), state, mean( score( score(:,1) == subID(isub), 3 ) ), nanmean( score( score(:,1) == subID(isub), 3 ) ) ];
    % No nanmean-ing to find subject who doesn't have score for both
end



switch nState
    case 2 % NC vs PD
        group(1).state = 1;
        group(2).state = [2, 3, 4, 5];
        legendName = {'NC', 'PD'};
        colourPlot = [1 .4 0; 0 .8 .8];
        
    case 3 % NC vs Med vs DBS
        group(1).state = 1;
        group(2).state = [2, 3];
        group(3).state = [4, 5];
        legendName = {'NC', 'Med', 'DBS'};
        colourPlot = [1 .4 0; 0 .8 0; 0 0 .8];

        
    case 5 % NC vs Med-on vs Med-off vs DBS-on vs DBS-off
        group(1).state = 1;
        group(2).state = 2;
        group(3).state = 3;
        group(4).state = 4;
        group(5).state = 5;
        legendName = {'NC', 'MON', 'MOF', 'DON', 'DOF'};
        colourPlot = [1 .4 0; 0 1 0; 0 .4 0; 0 0 1; 0 0 .4];

end

istate = 2;
tmp = mscoreSub( mscoreSub(:,2) == istate, 3 );
mscoreSubState(1,:) = [ istate, nanmean(tmp) , ...
    nanstd(tmp) / sqrt(sum(~isnan(tmp)))];



tmp = nan(nState,20);
for istate = 1:nState
    tmpScore = score(  ismember(score(:,1), subID ) & ismember(score(:,2), group(istate).state), 3 );
    tmp(istate, 1:length(tmpScore)) = tmpScore';
    mScore(istate) = nanmean(tmp(istate,:));
    nSub(istate) = length(tmp(istate,~isnan(tmp(istate,:))));
    semmScore(istate) = nanstd(tmp(istate,:))/sqrt(nSub(istate));
end




%% ttest

switch nState
    case 2
        [H12 P12] = ttest2( tmp( 1, ~isnan( tmp(1,:) ) ), tmp( 2, ~isnan( tmp(2,:) ) ) );
        
    case 5
        [H23 P23] = ttest2( tmp( 2, ~isnan( tmp(2,:) ) ), tmp( 3, ~isnan( tmp(3,:) ) ) );
        [H45 P45] = ttest2( tmp( 4, ~isnan( tmp(4,:) ) ), tmp( 5, ~isnan( tmp(5,:) ) ) );
end



figure(1),clf, set(gca,'FontSize', 20)
for istate = 1:nState
    bar(istate, mScore(istate))
    hold on
    errorbar(istate,mScore(istate),semmScore(istate),semmScore(istate), '.k', 'LineWidth', 2)
end

h = get(gca,'Children');

for istate = 1:nState
    iLegend(istate) = length(h)-(2*(istate-1));
    legendNamePlot{istate} = [legendName{istate} ':' num2str(nSub(istate))];
    set(h(length(h)-(2*(istate-1))), 'FaceColor', colourPlot(istate,:))
    set( h( length(h)-(2*(istate-1)+1) ), 'Color', [0 0 0])
end

legend(h(iLegend),legendNamePlot, 'Location', 'EastOutside')
set(gca, 'xTick', [1:nState], 'xTickLabel', legendName)

xlabel('state', 'FontSize', 20)
ylabel('UPDRS score', 'FontSize', 20)
title('meanUPDRS_eachState', 'interpret', 'none')




print('-dpng', [DIR.figIMResp 'alSub_UPDRS_manEachState_nState' num2syr(nState)])

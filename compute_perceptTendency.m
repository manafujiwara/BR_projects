function compute_perceptTendency(DIR,subID,Cond)
%%
load( [DIR.allSub '/IM_BPCln.mat'])
load( [ DIR.allSub '/subInfo.mat'  ] ,'rejSub', 'accSub', 'subInfo', 'badTri')


% ----------------------------------------------------------------------
% BPCln: 1)subID, 2)irun, 3)iblk, 4)iepc, 5)stimLabel,                 |
% 6)subjestState, 7)pBP (-1~1), 8)label, 9)fBP, 10)fcBP, 11:910)BP(t)  |
% ----------------------------------------------------------------------
    
d = BPCln( :, 1:8 ); %stimLabel 1=BR, 2=L, 3=R

n = 1:size(d,1)-1;

dd = [ diff(d(:,3))~=0, d(n,6), d(n,5), d(n+1,5), d(n,7:8), d(n+1,7:8), ];

ddd = dd( dd(:,1)~=1 & dd(:,4)==1, 2:end );

figure(1),clf
for istatus = 1:8
    
    switch istatus
        case 1
            gName = 'allSub';
            dPairs = ddd;
        case 2
            gName = 'Control';
            dPairs = ddd( ddd(:,1) == 1, : );
        case 3
            gName = 'Med-On';
            dPairs = ddd( ddd(:,1) == 2, : );
        case 4
            gName = 'Med-Off';
            dPairs = ddd( ddd(:,1) == 3, : );
        case 5
            gName = 'DBS-On';
            dPairs = ddd( ddd(:,1) == 4, : );
        case 6
            gName = 'DBS-Off';
            dPairs = ddd( ddd(:,1) == 5, : );
        case 7
            gName = 'PD-On';
            dPairs = ddd( ismember( ddd(:,1), [2,4]), : );
        case 8
            gName = 'PD-Off';
            dPairs = ddd( ismember( ddd(:,1), [3,5]), : );
    end
    

    %dPair~~; 1)subject state, 2)pBP(1st tri), 3)label(1st), 4)pBP(2nd), 5)label(2nd)
    dPair11 = dPairs( dPairs(:,2)==1, [1,4:end] );
    dPair1r1 = dPair11( dPair11(:,3)==1, : );
    dPair1l1 = dPair11( dPair11(:,3)==-1, : );
    dPair21 = dPairs( dPairs(:,2)==2, [1,4:end] );
    dPair31 = dPairs( dPairs(:,2)==3, [1,4:end] );
    
    mLabel11 = mean(dPair11(:,5));
    mLabel1r1 = mean(dPair1r1(:,5));
    mLabel1l1 = mean(dPair1l1(:,5));
    mLabel21 = mean(dPair21(:,5));
    mLabel31 = mean(dPair31(:,5));
    
    sdLabel11 = std(dPair11(:,5));
    sdLabel1r1 = std(dPair1r1(:,5));
    sdLabel1l1 = std(dPair1l1(:,5));
    sdLabel21 = std(dPair21(:,5));
    sdLabel31 = std(dPair31(:,5));
    
    args = [ mLabel21, mLabel31, mLabel11, mLabel1l1, mLabel1r1 ];
    sds  =  [ sdLabel21, sdLabel31, sdLabel11, sdLabel1l1, sdLabel1r1 ];
    
    H21 = ttest( dPair21(:,5) );
    H31 = ttest( dPair31(:,5) );
    H11 = ttest( dPair11(:,5) );
    H1l1 = ttest( dPair1l1(:,5) );
    H1r1 = ttest( dPair1r1(:,5) );

    Hs = [ H21, H31, H11, H1l1, H1r1 ];
    Hs(Hs==0) = NaN;
    
    figure(1), 
    subplot(2,4,istatus)
    bar(args)
    hold on
    h = errorbar( args, sds );
    h.Color = [0 0 0];
    h.LineStyle = 'none';
    h.LineWidth = 2;
    %     title('Percept in BR tri grouped by prev tri')
    plot([1:5], Hs .* args + 1.2 * sds, 'k*')

    title(gName)
    
    ylabel('Trend (L-R)')
    xlabel('Previous tri label')
    xticklabels( { 'Left', 'Right', 'BR(all)', 'BR(Left)', 'BR(Right)' } )
    xtickangle(30)
    set(gca, 'Fontsize', 14)
    
    
end

print( [ DIR.figBRPerceptDrag 'IM_BRPerceptDrag'] , '-dpng')


stimLabel = dPairs(:,2);
respLabel = dPairs(:,5);
subStatus = dPairs(:,1);

[P,T,STATS,TERMS] = anovan( dPairs(:,7), { stimLabel, respLabel } ,...
       'varnames', { 'stimLabel', 'respLabel' })
[P,T,STATS,TERMS] = anovan( dPairs(:,7), { stimLabel, respLabel, subStatus },...
    'varnames', { 'stimLabel', 'respLabel', 'subStatus' })








%% collectUPDRS
dbstop if error

[num,txt,raw] = xlsread('../../dataRaw/UPDRS.xlsx');

x = [ num(:,1:2), num(:,4:6) ]; 
age = x(:,1:3);
% x(x(:,2)==1 | x(:,5)==1 ,:) = [];% Remove control AND patients with only one session
x(x(:,3)==1 ,:) = []; %Remove controls
x(:,6) = nanmean(x(:,4:5),2);

% x: 1)subID, 2)age 3)subDrp, 4)UPDRS-on, 5)UPDRS-off 6)average UPDRS
save('../../config/UPDRS.mat','x','age')

for icmp = 1:2
    switch icmp
        case 1
            grp = 2; %Med
        case 2
            grp = 3; %DBS
    end
    y = x ( x(:,3) == grp,:);
    [H,P_UPDRS,CI,STATS_UPDRS] = ttest( y(:,4), y(:,5) );
    statSum_UPDRS(icmp,:) = [grp, P_UPDRS, STATS_UPDRS.tstat, STATS_UPDRS.df, STATS_UPDRS.sd];
end


% 
% nstate = size(raw,2) - 2;
% 
% for istate = 1:nstate
%     
%     subID = raw{ 1, 2 + istate};
%     
%     med = raw{ 2, 2 + istate };
%     dbs = raw{ 3, 2 + istate };
%     
%     if strcmp(med,'-') && strcmp(dbs,'-')
%         state = 'nc';
%         
%     elseif strcmp(med,'-') && dbs == 0
%         state = 'dof';
%         
%     elseif strcmp(med,'-')&& dbs == 1
%         state = 'don';
%         
%     elseif med == 0 && strcmp(dbs,'-')
%         state = 'mof';
%         
%     elseif med == 1 && strcmp(dbs,'-')
%         state = 'mon';
%         
%     else
%         error('incorrect state discription')
%         
%     end
%     
%     eachScore =  [raw{ 4:36, istate + 2 }]';
%     
%     UPDRS(istate).subID = subID;
%     UPDRS(istate).state = state;
%     UPDRS(istate).eachScore = eachScore;
%     UPDRS(istate).totalScore = raw{ 37, istate + 2 };
%     
% end
% 
% 

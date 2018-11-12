function check_vel_a_slwPhs_a_epIM(DIR)
load( [DIR.allSub '/IM_eyeAllSub.mat'] ,'iOKN', 'vel', 'slwPhsBox', 'slwPhsMed','subID')
load( [DIR.allSub '/IM_respAllSub.mat'] , 'lResp', 'rResp', 'resp', 'BPBox', 'lBPBox', 'rBPBox', 'BPMed', 'lBPMed', 'rBPMed', 'subID')

load([DIR.allSub '/IM_BPinfo.mat'], 'BPInfo')
load([DIR.allSub '/IM_eyeinfo.mat'], 'slwPhsBoxInfo', 'slwPhsMedInfo', 'iOKNInfo', 'velInfo')
% -------------------------------------------------------------------------
% BPinfo: 1)subID, 2)irun, 3)iblk, 4)iepc, 5)stimLabel,                   |
% 6)subjestState, 7)BPP (-1~1), 8)label, 9)fBP, 10)fcBP, 11:910)BP(t)     |
% -------------------------------------------------------------------------
% eyeInfo: 1)subID, 2)irun, 3)iblk, 4)iepc, 5)stimLabel, 6)subjestState   |
% 7)BPP (-1~1), 8)label, 9)fBP, 10)fcBP, 11:910)eyeinfo(t)                |
% -------------------------------------------------------------------------

%%  Choose c109 as representative subject
subID = 159;
time = linspace(-1,2,900);
% ---------------------------------------------------------------
% 1)subID, 2)irun, 3)iblk, 4)iepc, 5)stimLabel, 6)subjestState  |
% ---------------------------------------------------------------

iplot = 0;
figure(1),clf
figure(2),clf

for idata = 1:6
    
    switch idata %1-3:BP, 4-6:gaze
        case 1
            d = resp;
            tit = 'BP-NoSmooth';
            c = {'k'};
        case 2
            d = BPBox;
            tit = 'BP-Boxcar';
            c = {'r'};
        case 3
            d = BPMed;
            tit = 'BP-Median';
            c = {'b'};
        case 4
            iplot = 0;
            figure(1),clf
            figure(2),clf
            d = velInfo;
            tit = 'Vel-NoSmooth';
            c = {'k'};
        case 5
            d = slwPhsBoxInfo;
            tit = 'Vel-Boxcar';
            c = {'r'};
        case 6
            d = slwPhsMedInfo;
            tit = 'Vel-Median';
            c = {'b'};
    end
    
    if idata<=3
        i = 7:906;
        ylab = 'Resp Direction';
        yLim = [-1 1];
    elseif idata >3
        i = 11:910;
        ylab = 'Speed (deg/sec)';
        yLim = [-10 15];
    end
    
    
    for istim = 1:2
        iplot = iplot + 1;
        
        switch istim
            case 1 % label based on physical stimulus
                col = 5;
                lr = 2; % label r
                ll = 3; % label l
                sr = 2; % stim r
                sl = 3; % stim l
                dd = d;
            case 2 % label based on response 
                col = 8;
                lr = 1;
                ll = -1;
                sr = 1;
                sl = -1;
                
                if ismember(idata,[1,2,3])
                    dd = BPInfo;
                else
                    dd = d;
                end
        end
        
        
        rAmp = d( dd(:,1) == subID & dd(:,5) == sr & dd(:,col) == lr, i );
        lAmp = d( dd(:,1) == subID & dd(:,5) == sl & dd(:,col) == ll, i );
        
        
        cons = cat(1, rAmp, lAmp*-1 );
        
        md(iplot,:) = nanmean(cons,1);
        stdd(iplot, :) = nanstd(cons,1);
        
        figure(1)
        subplot(3,2,iplot)
        plot([-1 2], [0 0], '--k')
        hold on
        plot([0 0], yLim, '--k')

        h = plot(time, cons);
        ylim(yLim)
        title(tit)
        xlabel('Time from Stim Onset(s)')
        ylabel(ylab)
        set(gca,'FontSize', 23)
        
        figure(2)
        if idata == 1
            % Plot randomly selected trial as an exemplar
            trinum(istim) = randperm(size(cons,1),1);
        end
        
        subplot(3,2,istim)
        if ismember(iplot,[1,2])
            plot([-1 2], [0 0], '--k')
            hold on
            plot([0 0], yLim, '--k')
        end
        
        h = plot(time,cons( trinum(istim),:), ['' c{:} '']);
        set(h,'LineWidth', 1.5)
        ylim(yLim)
        xlim([-1 2])
        
        xlabel('Time from Stim Onset(s)')
        ylabel(ylab)
        set(gca,'FontSize', 23)
    end
    
    
    if idata == 3
        set(figure(1),'PaperUnits','inches','PaperPosition',[0 0 25 15])
        print('-dpng', [DIR.figGrpCmp '/' num2str(subID) 'cmpFiltering_BP.png'] )
        
        set(figure(2),'PaperUnits','inches','PaperPosition',[0 0 25 15])
        print('-dpng', [DIR.figGrpCmp '/' num2str(subID) 'cmpFiltering_BP_3Tri.png'] )
        
        yLim2 = [-.5 1.5];
        figure(3),clf
        c = {'k', 'k', 'r', 'r', 'b','b'};
        for i = 1:6
            switch i
                case {1,2}
                    tit = 'BP-NoSmooth';
                case {3,4}
                    tit = 'BP-Boxcar';
                case {5,6}
                    tit = 'BP-Median';
            end
            
            subplot(4,2,i)
            h = shadedErrorBar(time, md(i,:), stdd(i,:),['' c{i} ''] , 1);
            hold on
            plot([-1 2], [0 0], '--k')
            plot([0 0], yLim2, '--k')
            ylim( yLim2)
            xlim([-1 2])
            ylabel('BP Consistency')
            title(tit)
            set(gca,'FontSize', 23)
        end
        
        c = {'k', 'k', 'r', 'r', 'b','b'};
        for i = 1:6
            switch i
                case {1,3,5}
                    subplot(4,2,7)
                case {2,4,6}
                    subplot(4,2,8)
            end
            
            h = shadedErrorBar(time, md(i,:), stdd(i,:),['' c{i} ''] , 1);
            hold on
            plot([-1 2], [0 0], '--k')
            plot([0 0], yLim2, '--k')
            ylim( yLim2 )
            xlim([-1 2])
            xlabel('Time from Stim Onset(s)')
            ylabel('BP Consistency')
            set(gca,'FontSize', 23)
        end
        
        set(figure(3),'PaperUnits','inches','PaperPosition',[0 0 25 20])
        print('-dpng', [DIR.figGrpCmp '/' num2str(subID) 'cmpFiltering_mBP.png'] )
        
        
    elseif idata == 6
        set(figure(1),'PaperUnits','inches','PaperPosition',[0 0 25 15])
        print('-dpng', [DIR.figGrpCmp '/' num2str(subID) 'cmpFiltering_vel.png'] )
        
        set(figure(2),'PaperUnits','inches','PaperPosition',[0 0 25 15])
        print('-dpng', [DIR.figGrpCmp '/' num2str(subID) 'cmpFiltering_vel_3Tri.png'] )
        
        
        figure(3),clf
        yLim2 = [-3 10];
        c = {'k', 'k', 'r', 'r', 'b','b'};
        for i = 1:6
            switch i
                case {1,2}
                    tit = 'Vel-NoSmooth';
                case {3,4}
                    tit = 'Vel-Boxcar';
                case {5,6}
                    tit = 'Vel-Median';
            end
            
            subplot(4,2,i)
            h = shadedErrorBar(time, md(i,:), stdd(i,:),['' c{i} ''] , 1);
            hold on
            plot([-1 2], [0 0], '--k')
            plot([0 0], yLim2, '--k')
            ylim( yLim2)
            xlim([-1 2])
            ylabel('Speed (deg/sec)')
            title(tit)
            set(gca,'FontSize', 23)
            
        end
        
        c = {'k', 'k', 'r', 'r', 'b','b'};
        for i = 1:6
            switch i
                case {1,3,5}
                    subplot(4,2,7)
                case {2,4,6}
                    subplot(4,2,8)
            end
            
            h = shadedErrorBar(time, md(i,:), stdd(i,:),['' c{i} ''] , 1);
            hold on
            plot([-1 2], [0 0], '--k')
            plot([0 0], yLim2, '--k')
            ylim(yLim2)
            xlim([-1 2])
            xlabel('Time from Stim Onset(s)')
            ylabel('Speed (deg/sec)')
            set(gca,'FontSize', 23)
        end
        
        set(figure(3),'PaperUnits','inches','PaperPosition',[0 0 25 20])
        print('-dpng', [DIR.figGrpCmp '/' num2str(subID) 'cmpFiltering_mVel.png'] )
        
        
    end
    
end

end

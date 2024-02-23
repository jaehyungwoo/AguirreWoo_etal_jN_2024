function plot_metrics_timeseries_subjects(all_output, taskType, group_label,show_plots)

    plot_Behavioral_met = show_plots(1);
    plot_Cond_ent = show_plots(2);
    plot_Mutual_info = show_plots(3);
    plot_supp_probs = show_plots(4);
        
if plot_Behavioral_met    
    %% 1. behavioral metrics
    parameters = ["pwin","pbetter","pstay", "winstay","loseswitch","matching_measure", "RI_BW","RI_B","RI_W"];
    parameterLabels = ["prob(Win)","prob(Better)","prob(Stay)", "Win-Stay","Lose-Switch","dev.matching", "RI_{BW}","RI_B","RI_W"];

    figure; clf
    set(gcf,'Units','normalized','Position',[0,0,1,1],'color','w');
    sgtitle(taskType+" task");

    for i = 1:numel(parameters)
        subplot(3,3,i);
%     for g = 1%:numel(lesion_groups)
%         group_label = lesion_groups{g};
        thisDat = all_output.(taskType).(group_label);
        subLabels = unique(thisDat.subjectLabel);
        DataToPlot = thisDat.(parameters(i));
        
        % set x-axis: each phase
        xx = [0:8];  % 0=D, 1=R1(100/0), 2=R2(100/0), 3=R1(90/10), 4=R2(90/10), ..., 7=R1(70/30), 8=R2(70/30)
        % plot data for each subject
        for s = 1:numel(subLabels)
            thisSub_idx = (thisDat.subjectLabel==subLabels(s));
            yy = nan(size(xx)); SD = nan(size(xx)); t = 0;
            for x = [xx]
                t = t + 1;
                % set phase index
                if x==0
                    thisPhase_idx = logical(thisDat.D);     % Discrimination phase
                elseif mod(x,2)==1
                    thisPhase_idx = logical(thisDat.R1);    % R1
                elseif mod(x,2)==0  
                    thisPhase_idx = logical(thisDat.R2);    % R2
                end
                % set schedule index
                if x<=2
                    thisPhase_idx = thisPhase_idx&thisDat.prob1000;
                elseif x<=4
                    thisPhase_idx = thisPhase_idx&thisDat.prob9010;
                elseif x<=6
                    thisPhase_idx = thisPhase_idx&thisDat.prob8020;
                elseif x<=8
                    thisPhase_idx = thisPhase_idx&thisDat.prob7030;
                end
                thisIdx = thisSub_idx&thisPhase_idx;
                yy(t) = mean(DataToPlot(thisIdx));
                SD(t) = std(DataToPlot(thisIdx));
            end
            E = errorbar(xx,yy,SD,'LineWidth',1); hold on;
%             E.Color = lesion_colors{g};
        end
            
        ax = gca;
        ax.XLim = [-.1 9];
        ax.XTick = xx;
        xrow1 = {'D','R1','R2','R1','R2','R1','R2','R1','R2'};
        xrow2 = {'','100/0','','90/10','','80/20','','70/30',''};
        labelArray = [xrow1; xrow2];
        labelArray = strjust(pad(labelArray),'center');
        xlbls = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
        ax.XTickLabel = xlbls;
        ax.YLabel.String = parameterLabels(i);
        
        if i<=5; ax.YLim = [0 1]; end
        ax.Box = 'off';
        if i==1; title([group_label(1:3),' - ',group_label(5:end)]); end
        set(ax,'FontName','Helvetica','FontSize',15,'FontWeight','normal','LineWidth',3)
        
%     end
    end
end

if plot_Cond_ent
%% 2. Conditional Entropy
    parameters = ["ERDS", "ERDS_win", "ERDS_lose"; "EODS", "EODS_better", "EODS_worse"; "ERODS", "ERODS_winbetter", "ERODS_winworse"; "BLANK", "ERODS_losebetter", "ERODS_loseworse"];
    parameterLabels = ["ERDS", "ERDS_+", "ERDS_-"; "EODS", "EODS_B", "EODS_W"; 'ERODS', "ERODS_{B+}", "ERODS_{W+}"; 'BLANK', "ERODS_{B-}", "ERODS_{W-}"];
    
    if strcmp(taskType,"SO")
        parameterLabels(1,1) = "ERDS_{Stim}";
    elseif strcmp(taskType,"AO")
        parameterLabels(1,1) = "ERDS_{Loc}";
    end
    
    figure; clf
    set(gcf,'Units','normalized','Position',[0,0,1,1],'color','w');
    sgtitle(taskType+" task");

for jj = 1:size(parameters,1)
    for kk = 1:size(parameters,2)
        if (jj==4 && kk==1); continue; end  
        i = kk + (jj-1)*3;

        subplot(4,3,i);
        thisDat = all_output.(taskType).(group_label);
        subLabels = unique(thisDat.subjectLabel);
        DataToPlot = thisDat.(parameters(jj,kk));
        
        % set x-axis: each phase
        xx = [0:8];  % 0=D, 1=R1(100/0), 2=R2(100/0), 3=R1(90/10), 4=R2(90/10), ..., 7=R1(70/30), 8=R2(70/30)
        % plot data for each subject
        for s = 1:numel(subLabels)
            thisSub_idx = (thisDat.subjectLabel==subLabels(s));
            yy = nan(size(xx)); SD = nan(size(xx)); t = 0;
            for x = [xx]
                t = t + 1;
                % set phase index
                if x==0
                    thisPhase_idx = logical(thisDat.D);     % Discrimination phase
                elseif mod(x,2)==1
                    thisPhase_idx = logical(thisDat.R1);    % R1
                elseif mod(x,2)==0  
                    thisPhase_idx = logical(thisDat.R2);    % R2
                end
                % set schedule index
                if x<=2
                    thisPhase_idx = thisPhase_idx&thisDat.prob1000;
                elseif x<=4
                    thisPhase_idx = thisPhase_idx&thisDat.prob9010;
                elseif x<=6
                    thisPhase_idx = thisPhase_idx&thisDat.prob8020;
                elseif x<=8
                    thisPhase_idx = thisPhase_idx&thisDat.prob7030;
                end
                thisIdx = thisSub_idx&thisPhase_idx;
                yy(t) = nanmean(DataToPlot(thisIdx));
                SD(t) = nanstd(DataToPlot(thisIdx));
            end
            E = errorbar(xx,yy,SD,'LineWidth',1); hold on;
%             E.Color = lesion_colors{g};
        end
            
        ax = gca;
        ax.XLim = [-.1 9];
        ax.XTick = xx;
        xrow1 = {'D','R1','R2','R1','R2','R1','R2','R1','R2'};
        xrow2 = {'','100/0','','90/10','','80/20','','70/30',''};
        labelArray = [xrow1; xrow2];
        labelArray = strjust(pad(labelArray),'center');
        xlbls = strtrim(sprintf('%s\\newline%s\n', labelArray{:}));
        ax.XTickLabel = xlbls;
        ax.YLabel.String = parameterLabels(jj,kk);
        
        ax.YLim = [0 1];
        ax.Box = 'off';
        if i==1; title([group_label(1:3),' - ',group_label(5:end)]); end
        set(ax,'FontName','Helvetica','FontSize',15,'FontWeight','normal','LineWidth',3)
        
    end
end
    fff = gcf;
    mysub = fff.Children(5);
    mysub.Position = mysub.Position - [0 0.12 0 0];
    
end
end
function plot_metrics_timeseries_group_mean(all_output, taskType, lesion_region, show_plots)
    
    plot_Behavioral_met = show_plots(1);
    plot_Cond_ent = show_plots(2);
    plot_Mutual_info = show_plots(3);
    plot_supp_probs = show_plots(4);
    
lesion_groups = ["ACC_eGFP", "ACC_hm4Di", ...
                 "BLA_eGFP", "BLA_hm4Di", ...
                 "OFC_eGFP", "OFC_hm4Di"]; 
lesion_colors = {[0 0.4470 0.7410],[0.0157 0.2157 0.3490],...
                [0.6350 0.0780 0.1840],[0.3686 0.0980 0.1451],...
                [0.9290 0.6940 0.1250],[0.4706 0.3412 0.0431]};
            
% select lesion_groups to plot             
f = [find(contains(lesion_groups,lesion_region))];
lesion_groups = lesion_groups(f);            
lesion_colors = lesion_colors(f);    % These will be R1-CNO and R2-CNO colors     

% if numel(lesion_groups)~=2; error("Check # of groups"); end
% if numel(lesion_colors)~=2; error("Check # of groups (colors)"); end
% disp([lesion_colors{1}, lesion_colors{2}]);

markersSet = {'o','s'};

if plot_Behavioral_met    
    tic
    %% 1. behavioral metrics
    parameters = ["pwin","pbetter","pstay", "winstay","loseswitch","matching_measure", "RI_BW","RI_B","RI_W"];
    parameterLabels = ["prob(Win)","prob(Better)","prob(Stay)", "Win-Stay","Lose-Switch","dev.matching", "RI_{BW}","RI_B","RI_W"];

    figure; clf
    set(gcf,'Units','normalized','Position',[0,0,1,1],'color','w');
    sgtitle(taskType+" task");

for i = 1:numel(parameters)
    subplot(3,3,i);
    for g = 1:numel(lesion_groups)
        group_label = lesion_groups(g);
        if contains(group_label,"eGFP") 
            groupColor = 'g';               % green for fluorescent
        end
        
        thisDat = all_output.(taskType).(group_label);
        subLabels = unique(thisDat.subjectLabel);   
        N(g) = numel(subLabels);                % track subject N for title
        DataToPlot = thisDat.(parameters(i));
        
        % divide hm4Di group into two: CNO first or VEH first
        CNO_first_idx = false(numel(subLabels),1);
        VEH_first_idx = false(numel(subLabels),1);
        for s = 1:numel(subLabels)
            thisSub_idx = (thisDat.subjectLabel==subLabels(s));
            if sum(thisDat.R1(thisSub_idx)==thisDat.CNO(thisSub_idx))==sum(thisSub_idx)
                CNO_first_idx(s) = true;
            elseif sum(thisDat.R1(thisSub_idx)==thisDat.VEH(thisSub_idx))==sum(thisSub_idx)
                VEH_first_idx(s) = true;
            end
        end
        CNO12_groups = {subLabels(CNO_first_idx), subLabels(VEH_first_idx)};
        
        % plot data for each (two) group
        for S = 1:2
            thisGroup = CNO12_groups{S};
            NN{g}(S) = numel(thisGroup);    % track subject N for legend
            if ~contains(group_label,"eGFP"); groupColor = lesion_colors{S}; end
            
            % set x-axis: 0=D, 1=R1(100/0), 2=R2(100/0), ..., 7=R1(70/30), 8=R2(70/30)
            xx = [0:8];
            yy = {};
            clear SEM y E
            for s = 1:numel(thisGroup)
                thisSub_idx = (thisDat.subjectLabel==thisGroup(s));
                yy{s} = nan(size(xx)); t = 0; 
                for x = [xx]
                    t = t + 1;      if s==1; Y{t} = []; end
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
                    thisIdx = thisSub_idx&thisPhase_idx;    %disp(s+": "+sum(thisIdx));
                    yy{s}(t) = nanmean(DataToPlot(thisIdx));
%                     Y{t} = [Y{t}; DataToPlot(thisIdx)]; % bin for Method 2
                end
            end
            y = nanmean(cell2mat(yy'),1);   % METHOD 1: Mean of subject means
            SEM = std(cell2mat(yy'),0,1,'omitnan')./sqrt(sum(~isnan(cell2mat(yy')),1));  % divide by # of non-NaN valued subjects
            
            if isempty(thisGroup); y = nan(size(xx)); SEM = y; end          % if no data, fill with NaN (for legend to work)
            
            if length(y)~=length(xx);   error("check y length: "+length(y)); end
            if length(SEM)~=length(xx);   error("check y length: "+length(SEM)); end
            
            % 1. Shaded Error bar
%             E = shadedErrorBar(xx,y,SEM,'lineProps',{'LineWidth', 1.5, 'Color',groupColor}); hold on;
%             E.mainLine.MarkerSize = 8;
%             if S==1; E.mainLine.MarkerIndices = [2:2:9]; E.mainLine.Marker = 'o'; end      % CNO first
%             if S==2; E.mainLine.MarkerIndices = [3:2:9]; E.mainLine.Marker = 's'; end      % VEH first
%             if ~contains(group_label,"eGFP") 
%                 E.mainLine.MarkerEdgeColor = 'm';
%                 E.mainLine.MarkerSize = 10; 
%             end
            % Normal error bar
            E0 = errorbar(xx,y,SEM,'Color',groupColor,'LineStyle','none','Color',groupColor,'HandleVisibility','off'); hold on;
            P0 = plot(xx,y,'Color',groupColor,'LineWidth',2,'Marker',markersSet{S},'MarkerSize',8,'MarkerEdgeColor',groupColor,'MarkerFaceColor',groupColor,'MarkerIndices',[1+S:2:9]);
            if ~contains(group_label,"eGFP") 
                P0.MarkerEdgeColor = 'm'; %P0.MarkerFaceColor = P0.MarkerEdgeColor;
            end
            
    %         % METHOD 2: entire subjects average (not mean of means)
    %         for t = 1:length(Y)
    %            y2(t) = nanmean(Y{t}); SD(t) = nanstd(Y{t})/sqrt(length(Y{t}));
    %         end
    %         E2 = errorbar(xx,y2,SD); hold on;
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
    end
    
    ax.Box = 'off';
    set(ax,'FontName','Helvetica','FontSize',15,'FontWeight','normal','LineWidth',3)
    if i<=5; ax.YLim = [0 1]; plot(ax.XLim, [0.5, 0.5],':k'); end
    if i==1
        title([lesion_region,'   (N = ',num2str(N(1)),', ',num2str(N(2)),')']);
        L = legend({['R1-CNO (eGFP; n=',num2str(NN{1}(1)),')'],['R2-CNO (eGFP; n=',num2str(NN{1}(2)),')'], ...
            ['R1-CNO (hm4Di; n=',num2str(NN{2}(1)),')'],['R2-CNO (hm4Di; n=',num2str(NN{2}(2)),')']}); 
        L.Box = 'off';
    end
    
end
    toc
end
    
if plot_Cond_ent
    tic
%% 2. Conditional Entropy
    
    parameters = ["ERDS", "ERDS_win", "ERDS_lose"; "EODS", "EODS_better", "EODS_worse"; "ERODS", "ERODS_winbetter", "ERODS_winworse"; "BLANK", "ERODS_losebetter", "ERODS_loseworse"];
    parameterLabels = ["ERDS", "ERDS_+", "ERDS_-"; "EODS", "EODS_B", "EODS_W"; 'ERODS', "ERODS_{B+}", "ERODS_{W+}"; 'BLANK', "ERODS_{B-}", "ERODS_{W-}"];
 
    if strcmp(taskType,"SO")
        parameterLabels(1,1) = "ERDS_{Stim(-), Loc(:)}";
        
%         % additional plot: ERDS_loc together w/ ERDS_stim for SO task
%         parameters(4,1) = "ERDS_loc";
        
    elseif strcmp(taskType,"AO")
        parameterLabels(1,1) = "ERDS_{Loc}";
    end    
    
%     parameters = ["ERDS_loc", "ERDS_loc_win", "ERDS_loc_lose"; "EODS", "EODS_better", "EODS_worse"; "ERODS", "ERODS_winbetter", "ERODS_winworse"; "BLANK", "ERODS_losebetter", "ERODS_loseworse"];
%     parameterLabels = ["ERDS_{Loc}", "ERDS+_{Loc}", "ERDS_-_{Loc}"; "EODS", "EODS_B", "EODS_W"; 'ERODS', "ERODS_{B+}", "ERODS_{W+}"; 'BLANK', "ERODS_{B-}", "ERODS_{W-}"];
    
    
    
    figure; clf
    set(gcf,'Units','normalized','Position',[0,0,1,1],'color','w');
    sgtitle(taskType+" task");

for jj = 1:size(parameters,1)
    for kk = 1:size(parameters,2)
        i = kk + (jj-1)*3;
        
        if (jj==4 && kk==1)
%             if strcmp(taskType,"AO")
                continue; 
%             else
%                 % additional plot: ERDS_loc together w/ ERDS_stim for SO task
%                 i = 1;  % plot with ERDS_stim
%             end
        end  
        
        subplot(4,3,i);  hold on;
        for g = 1:numel(lesion_groups)
            group_label = lesion_groups(g);
            % green for fluorescent
            if contains(group_label,"eGFP"); groupColor = 'g'; end

            thisDat = all_output.(taskType).(group_label);
            subLabels = unique(thisDat.subjectLabel);
            N(g) = numel(subLabels);                % track subject N for title
            DataToPlot = thisDat.(parameters(jj,kk));

            % divide hm4Di group into two: CNO first or VEH first
            CNO_first_idx = false(numel(subLabels),1);
            VEH_first_idx = false(numel(subLabels),1);
            for s = 1:numel(subLabels)
                thisSub_idx = (thisDat.subjectLabel==subLabels(s));
                if sum(thisDat.R1(thisSub_idx)==thisDat.CNO(thisSub_idx))==sum(thisSub_idx)
                    CNO_first_idx(s) = true;
                elseif sum(thisDat.R1(thisSub_idx)==thisDat.VEH(thisSub_idx))==sum(thisSub_idx)
                    VEH_first_idx(s) = true;
                end
            end
            CNO12_groups = {subLabels(CNO_first_idx), subLabels(VEH_first_idx)};

            % plot data for each (two) group
            for S = 1:2
                thisGroup = CNO12_groups{S};
                NN{g}(S) = numel(thisGroup);    % track subject N for legend
                if ~contains(group_label,"eGFP"); groupColor = lesion_colors{S}; end

                % set x-axis: 0=D, 1=R1(100/0), 2=R2(100/0), ..., 7=R1(70/30), 8=R2(70/30)
                xx = [0:8];
                yy = {};
                clear SEM y E
                for s = 1:numel(thisGroup)
                    thisSub_idx = (thisDat.subjectLabel==thisGroup(s));
                    yy{s} = nan(size(xx)); t = 0; 
                    for x = [xx]
                        t = t + 1;      if s==1; Y{t} = []; end
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
                        thisIdx = thisSub_idx&thisPhase_idx;    %disp(s+": "+sum(thisIdx));
                        yy{s}(t) = nanmean(DataToPlot(thisIdx));
    %                     Y{t} = [Y{t}; DataToPlot(thisIdx)]; % bin for Method 2
                    end
                end
                y = nanmean(cell2mat(yy'),1);   % METHOD 1: Mean of subject means
                SEM = std(cell2mat(yy'),0,1,'omitnan')./sqrt(sum(~isnan(cell2mat(yy')),1));  % divide by # of non-NaN valued subjects
                if isempty(thisGroup); y = nan(size(xx)); SEM = y; end          % if no data, fill with NaN (for legend to work)
                % Shaded error bar
%                 E = shadedErrorBar(xx,y,SEM,'lineProps',{'LineWidth', 1, 'Color',groupColor}); hold on;
%                 E.mainLine.MarkerSize = 8; 
%                 if S==1; E.mainLine.MarkerIndices = [2:2:9]; E.mainLine.Marker = 'o'; end      % CNO first
%                 if S==2; E.mainLine.MarkerIndices = [3:2:9]; E.mainLine.Marker = 's'; end      % VEH first
%                 if ~contains(group_label,"eGFP") 
%                     E.mainLine.MarkerEdgeColor = 'm';
%                     E.mainLine.MarkerSize = 10; 
%                 end
                % Normal error bar
                E0 = errorbar(xx,y,SEM,'Color',groupColor,'LineStyle','none','Color',groupColor,'HandleVisibility','off'); hold on;
                P0 = plot(xx,y,'Color',groupColor,'LineWidth',2,'Marker',markersSet{S},'MarkerSize',8,'MarkerEdgeColor',groupColor,'MarkerFaceColor',groupColor,'MarkerIndices',[1+S:2:9]);
                if ~contains(group_label,"eGFP") 
                    P0.MarkerEdgeColor = 'm'; %P0.MarkerFaceColor = P0.MarkerEdgeColor;
                end
                if jj==4&&kk==1; P0.LineStyle = ':'; end

        %         % METHOD 2: entire subjects average (not mean of means)
        %         for t = 1:length(Y)
        %            y2(t) = nanmean(Y{t}); SD(t) = nanstd(Y{t})/sqrt(length(Y{t}));
        %         end
        %         E2 = errorbar(xx,y2,SD); hold on;
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
            if ~(jj==4&&kk==1); ax.YLabel.String = parameterLabels(jj,kk); end
        end
        ax.Box = 'off';
        set(ax,'FontName','Helvetica','FontSize',15,'FontWeight','normal','LineWidth',3)
        if i<=7; ax.YLim = [0 1]; end
        if i==1
            title([lesion_region,'   (N = ',num2str(N(1)),', ',num2str(N(2)),')']);
            L = legend({['R1-CNO (eGFP; n=',num2str(NN{1}(1)),')'],['R2-CNO (eGFP; n=',num2str(NN{1}(2)),')'], ...
                ['R1-CNO (hm4Di; n=',num2str(NN{2}(1)),')'],['R2-CNO (hm4Di; n=',num2str(NN{2}(2)),')']}); 
            L.Box = 'off';
        end
    end
end

toc
end

end
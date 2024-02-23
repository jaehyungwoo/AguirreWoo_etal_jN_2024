function plot_metrics_across_days_group_mean(all_output, taskType, group_label, phaseType, show_plots, shadedError, lesion_colors)

phaseTypes = {["D"],["R1","prob1000"],["R2","prob1000"],["R1","prob9010"],["R2","prob9010"],...
    ["R1","prob8020"],["R2","prob8020"],["R1","prob7030"],["R2","prob7030"]};

    plot_Behavioral_met = show_plots(1);
    plot_Cond_ent = show_plots(2);
    plot_Mutual_info = show_plots(3);
    plot_supp_probs = show_plots(4);
            
if contains(group_label,'eGFP')
    lesion_region = 'eGFP';
else
    lesion_region = group_label(1:3);
end            
            
if ~exist('lesion_colors','var')
    lesion_groups = ["ACC_eGFP", "ACC_hm4Di", ...
                 "BLA_eGFP", "BLA_hm4Di", ...
                 "OFC_eGFP", "OFC_hm4Di"]; 
    lesion_colors = {[0 0.4470 0.7410],[0.0157 0.2157 0.3490],...
                [0.6350 0.0780 0.1840],[0.3686 0.0980 0.1451],...
                [0.9290 0.6940 0.1250],[0.4706 0.3412 0.0431]};
    % select lesion_groups to plot             
    f = [find(contains(lesion_groups,lesion_region))];  
    lesion_colors = lesion_colors(f);    % These will be R1-CNO and R2-CNO colors     
end

markersSet = {'o','s'};

if plot_Behavioral_met    
    tic
    %% 1. behavioral metrics
    parameters = ["pwin","pbetter","pstay", "winstay","loseswitch","matching_measure", "RI_BW","RI_B","RI_W"];
    parameterLabels = ["P(Win)","P(Better)","P(Stay)", "Win-Stay","Lose-Switch","dev.matching", "RI_{BW}","RI_B","RI_W"];

    figure; clf
    set(gcf,'Units','normalized','Position',[0,0,.8,1],'color','w');

    for i = 1:numel(parameters)
        %% Subplot
        subplot(3,3,i);
        if contains(group_label,"eGFP") 
            groupColor = 'g';               % green for fluorescent
        end
        
        thisDat = all_output.(taskType).(group_label);
        subLabels = unique(thisDat.subjectLabel);   
        N = numel(subLabels);                % track subject N for title
        DataToPlot = thisDat.(parameters(i));
        
        numDays = 5;
        days_idx = cell(numDays,1); this_idx = []; 
        % find phase indices
        if length(phaseType)==2
            this_idx = thisDat.(phaseType(1))&thisDat.(phaseType(2));
        else
            this_idx = thisDat.(phaseType(1));
        end
        % find respective session day indices
        for d = 1:numDays
            days_idx{d} = this_idx&(thisDat.SessionDay==d);     %disp(find(days_idx{d}));
        end
        
        % divide eGFP (hm4Di) group into two: CNO first or VEH first
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
            NN(S) = numel(thisGroup);    % track subject N for legend
            if contains(group_label,"hm4Di"); groupColor = lesion_colors{S}; end
            
            % set x-axis: session day 1 through 5
            xx = [1:numDays];
            yy = {};
            clear SEM y E
            for s = 1:numel(thisGroup)
                thisSub_idx = (thisDat.subjectLabel==thisGroup(s));     %disp(sum(thisSub_idx));
                yy{s} = nan(size(xx)); t = 0; 
                for x = [xx]
                    if s==1; Y{x} = []; end
                    % set subject session day indices
                    thisDayIdx = thisSub_idx&days_idx{x};
                    
                    yy{s}(x) = nanmean(DataToPlot(thisDayIdx));
%                     Y{x} = [Y{x}; DataToPlot(thisIdx)]; % bin for Method 2
                end
            end
            y = nanmean(cell2mat(yy'),1);   % METHOD 1: Mean of subject means
            SEM = std(cell2mat(yy'),0,1,'omitnan')./sqrt(sum(~isnan(cell2mat(yy')),1));  % divide by # of non-NaN valued subjects
            
            if isempty(thisGroup); y = nan(size(xx)); SEM = y; end          % if no data, fill with NaN (for legend to work)
            
            if length(y)~=length(xx);   error("check y length: "+length(y)); end
            if length(SEM)~=length(xx);   error("check y length: "+length(SEM)); end
            
            if shadedError
            % 1. Shaded Error bar
                E = shadedErrorBar(xx,y,SEM,'lineProps',{'LineWidth',1.5,'Color',groupColor}); hold on;
                if (S==1&&strcmp(phaseType(1),"R1"))||(S==2&&strcmp(phaseType(1),"R2"))
%                     E.patch.FaceColor = 'm';     
                end
            else
                % Normal error bar
                E0 = errorbar(xx,y,SEM,'Color',groupColor,'Color',groupColor,'LineWidth',1.5); hold on;
    %             E0 = errorbar(xx,y,SEM,'Color',groupColor,'LineWidth',1.5,'LineStyle','none','Color',groupColor,'HandleVisibility','off'); hold on;
    %             P0 = plot(xx,y,'Color',groupColor,'LineWidth',2,'Marker',markersSet{S},'MarkerSize',8,'MarkerEdgeColor',groupColor,'MarkerFaceColor',groupColor,'MarkerIndices',[1+S:2:9]);
                if ~contains(group_label,"eGFP") 
                    E0.MarkerEdgeColor = 'm'; %P0.MarkerFaceColor = P0.MarkerEdgeColor;
                end
            end

        end
        
        ax = gca;
        ax.Box = 'off';
        ax.XLim = [1 numDays];
        ax.XTick = xx;
        ax.YLabel.String = parameterLabels(i);
        set(ax,'FontName','Helvetica','FontSize',15,'FontWeight','normal','LineWidth',3)
        if i==1
            L = legend({['R1-CNO (n=',num2str(NN(1)),')'],['R2-CNO (n=',num2str(NN(2)),')']},'Box','off');
            %title([group_label(1:3),'   (N = ',num2str(N),')']);
        end   
        if i==2
            L = legend({['R1-CNO (n=',num2str(NN(1)),')'],['R2-CNO (n=',num2str(NN(2)),')']},'Box','off');
            L.Location = 'northwest';
            text(1,1,phaseType(1),'units','normalized','FontWeight','bold','FontSize',20,'HorizontalAlignment','right');
        end
        if i<=5; ax.YLim = [0 1]; plot(ax.XLim, [0.5, 0.5],':k','HandleVisibility','off'); end
        if i>=7; ax.XLabel.String = "Days"; end
    
    end
    sgtitle(taskType+" task: "+ phaseType(1)+"_{100/0}");
    if length(phaseType)==2
        if strcmp(phaseType(2),"prob1000")
            schd_lbl = "100/0";
        else
            schd_lbl = convertStringsToChars(phaseType(2)); schd_lbl = schd_lbl(5:end); schd_lbl = [schd_lbl(1:2),'/',schd_lbl(3:4)];
            schd_lbl = convertCharsToStrings(schd_lbl);
        end
        sgtitle(group_label(1:3)+"   (N = "+num2str(N)+"): "+taskType+" task: "+phaseType(1)+"_{"+schd_lbl+"}");
    end
    
    toc
end
    
if plot_Cond_ent
    tic
    %% 2. Conditional Entropy
    parameters = ["ERDS", "ERDS_win", "ERDS_lose"; "EODS", "EODS_better", "EODS_worse"; "ERODS", "ERODS_winbetter", "ERODS_winworse"; "BLANK", "ERODS_losebetter", "ERODS_loseworse"];
    parameterLabels = ["ERDS", "ERDS_+", "ERDS_-"; "EODS", "EODS_B", "EODS_W"; 'ERODS', "ERODS_{B+}", "ERODS_{W+}"; 'BLANK', "ERODS_{B-}", "ERODS_{W-}"];
     
    if strcmp(taskType,"SO")
        parameterLabels(1,1) = "ERDS_{Stim}";
%         % additional plot: ERDS_loc together w/ ERDS_stim for SO task
%         parameterLabels(1,1) = "ERDS_{Stim(-), Loc(:)}";
%         parameters(4,1) = "ERDS_loc";
    elseif strcmp(taskType,"AO")
        parameterLabels(1,1) = "ERDS_{Loc}";
    end    
    
%     parameters = ["ERDS_loc", "ERDS_win", "ERDS_lose"; "EODS", "EODS_better", "EODS_worse"; "ERODS", "ERODS_winbetter", "ERODS_winworse"; "BLANK", "ERODS_losebetter", "ERODS_loseworse"];
%     parameterLabels = ["ERDS_{loc}", "ERDS_+", "ERDS_-"; "EODS", "EODS_B", "EODS_W"; 'ERODS', "ERODS_{B+}", "ERODS_{W+}"; 'BLANK', "ERODS_{B-}", "ERODS_{W-}"];
    
    figure; clf
    set(gcf,'Units','normalized','Position',[0,0,.65,1],'color','w');

for jj = 1:size(parameters,1)
    for kk = 1:size(parameters,2)
        i = kk + (jj-1)*3;
        
        if (jj==4 && kk==1)
%             if strcmp(taskType,"AO")
                continue; 
% %             else
%                 % additional plot: ERDS_loc together w/ ERDS_stim for SO task
%                 i = 1;  % plot with ERDS_stim
%             end
        end  
        
        %% Subplot
        subplot(4,3,i); hold on
        if contains(group_label,"eGFP") 
            groupColor = 'g';               % green for fluorescent
        end
        
        thisDat = all_output.(taskType).(group_label);
        subLabels = unique(thisDat.subjectLabel);   
        N = numel(subLabels);                % track subject N for title
        DataToPlot = thisDat.(parameters(jj,kk));
        
        numDays = 5;
        days_idx = cell(numDays,1); this_idx = []; 
        % find phase indices
        if length(phaseType)==2
            this_idx = thisDat.(phaseType(1))&thisDat.(phaseType(2));
        else
            this_idx = thisDat.(phaseType(1));
        end
        % find respective session day indices
        for d = 1:numDays
            days_idx{d} = this_idx&(thisDat.SessionDay==d);     % disp(find(days_idx{d}));
        end
        
        % divide eGFP (hm4Di) group into two: CNO first or VEH first
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
            NN(S) = numel(thisGroup);    % track subject N for legend
            if contains(group_label,"hm4Di"); groupColor = lesion_colors{S}; end
            
            % set x-axis: session day 1 through 5
            xx = [1:numDays];
            yy = {};
            clear SEM y E
            for s = 1:numel(thisGroup)
                thisSub_idx = (thisDat.subjectLabel==thisGroup(s));     %disp(sum(thisSub_idx));
                yy{s} = nan(size(xx)); t = 0; 
                for x = [xx]
                    if s==1; Y{x} = []; end
                    % set subject session day indices
                    thisDayIdx = thisSub_idx&days_idx{x};
                    
                    yy{s}(x) = nanmean(DataToPlot(thisDayIdx));
%                     Y{x} = [Y{x}; DataToPlot(thisIdx)]; % bin for Method 2
                end
            end
            y = nanmean(cell2mat(yy'),1);   % METHOD 1: Mean of subject means
            SEM = std(cell2mat(yy'),0,1,'omitnan')./sqrt(sum(~isnan(cell2mat(yy')),1));  % divide by # of non-NaN valued subjects
            
            if isempty(thisGroup); y = nan(size(xx)); SEM = y; end          % if no data, fill with NaN (for legend to work)
            
            if length(y)~=length(xx);   error("check y length: "+length(y)); end
            if length(SEM)~=length(xx);   error("check y length: "+length(SEM)); end
            
            if shadedError
            % 1. Shaded Error bar
                E = shadedErrorBar(xx,y,SEM,'lineProps',{'LineWidth',1.5,'Color',groupColor}); hold on;
                if (S==1&&strcmp(phaseType(1),"R1"))||(S==2&&strcmp(phaseType(1),"R2"))
%                     if ~isnan(y); E.patch.FaceColor = 'm'; end
                end
            else
            % Normal error bar
            E0 = errorbar(xx,y,SEM,'Color',groupColor,'LineWidth',1.5); hold on;
    %             E0 = errorbar(xx,y,SEM,'Color',groupColor,'LineStyle','none','Color',groupColor,'HandleVisibility','off'); hold on;
    %             P0 = plot(xx,y,'Color',groupColor,'LineWidth',2,'Marker',markersSet{S},'MarkerSize',8,'MarkerEdgeColor',groupColor,'MarkerFaceColor',groupColor,'MarkerIndices',[1+S:2:9]);
                if ~contains(group_label,"eGFP") 
                    E0.MarkerEdgeColor = 'm'; %P0.MarkerFaceColor = P0.MarkerEdgeColor;
                end
            end
        end
        
        ax = gca;
        ax.XLim = [1 numDays];
        ax.XTick = xx;
        ax.YLabel.String = parameterLabels(jj,kk);
        if i>=7; ax.XLabel.String = "Days"; end
    
        ax.Box = 'off';
        set(ax,'FontName','Helvetica','FontSize',15,'FontWeight','normal','LineWidth',3)
%         if i<=5; ax.YLim = [0 1]; plot(ax.XLim, [0.5, 0.5],':k','HandleVisibility','off'); end
        if i==1
            ylim([.5 1.0]);
            L = legend({['R1-CNO (n=',num2str(NN(1)),')'],['R2-CNO (n=',num2str(NN(2)),')']},'Box','off');
            L.Location = 'northwest';
            text(1,1,phaseType(1),'units','normalized','FontWeight','bold','FontSize',18,'HorizontalAlignment','right');
        end
    
    end
end
    sgtitle(taskType+" task: "+ phaseType(1)+"_{100/0}");
    if length(phaseType)==2
        if strcmp(phaseType(2),"prob1000")
            schd_lbl = "100/0";
        else
            schd_lbl = convertStringsToChars(phaseType(2)); schd_lbl = schd_lbl(5:end); schd_lbl = [schd_lbl(1:2),'/',schd_lbl(3:4)];
            schd_lbl = convertCharsToStrings(schd_lbl);
        end
        sgtitle(group_label(1:3)+"   (N = "+num2str(N)+"): "+taskType+" task: "+phaseType(1)+"_{"+schd_lbl+"}");
    end
toc
end

end
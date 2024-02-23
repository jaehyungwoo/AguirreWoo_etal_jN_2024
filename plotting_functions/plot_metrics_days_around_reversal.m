function plot_metrics_days_around_reversal(all_output, taskType, group_label, reversalNum, show_plots, shadedError, lesion_colors)
    plot_Behavioral_met = show_plots(1);
    plot_Cond_ent = show_plots(2);
    plot_Mutual_info = show_plots(3);
    plot_supp_probs = show_plots(4);
    
phaseSet = {["D"],["R1","prob1000"],["R2","prob1000"],["R1","prob9010"],["R2","prob9010"],...
    ["R1","prob8020"],["R2","prob8020"],["R1","prob7030"],["R2","prob7030"]};

if reversalNum>numel(phaseSet)-1
    error("Reversal number error: set an integer between [1, 8]");
end

phaseType1 = phaseSet{reversalNum};       % before rev
phaseType2 = phaseSet{reversalNum+1};     % after rev

if reversalNum<=3
    schd_lbl1 = "100/0";
    schd_lbl2 = "100/0";
    if reversalNum==3; schd_lbl2 = "90/10"; end
else
    schd_lbl1 = convertStringsToChars(phaseType1(2)); schd_lbl1 = schd_lbl1(5:end); schd_lbl1 = [schd_lbl1(1:2),'/',schd_lbl1(3:4)];
    schd_lbl1 = convertCharsToStrings(schd_lbl1);
    schd_lbl2 = convertStringsToChars(phaseType2(2)); schd_lbl2 = schd_lbl2(5:end); schd_lbl2 = [schd_lbl2(1:2),'/',schd_lbl2(3:4)];
    schd_lbl2 = convertCharsToStrings(schd_lbl2);
end
    
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
group_linetypes = {'-','-'};

if plot_Behavioral_met    
    tic
    %% 1. behavioral metrics
    parameters = ["pwin","pbetter","pstay", "winstay","loseswitch","matching_measure", "RI_BW","RI_B","RI_W"];
    parameterLabels = ["prob(Win)","prob(Better)","prob(Stay)", "Win-Stay","Lose-Switch","dev.matching", "RI_{BW}","RI_B","RI_W"];

    figure; clf
    set(gcf,'Units','normalized','Position',[0,0,.8,1],'color','w');
    sgtitle(taskType+" task, Rev#"+reversalNum+": "+ phaseType1(1)+"_{"+schd_lbl1+"}-to-" + phaseType2(1)+"_{"+schd_lbl2+"}");

    for i = 1:numel(parameters)
        %% Subplot
        subplot(3,3,i);
        if contains(group_label,"eGFP") 
            groupColor = 'g';               % green for fluorescent
        end
        
        thisDat = all_output.(taskType).(group_label);
        animalLabels = unique(thisDat.subjectLabel);   
        N = numel(animalLabels);                % track subject N for title
        DataToPlot = thisDat.(parameters(i));
        
        numDays = 3;    % last 3 days before rev, first 3+1 days after rev
        days_idx = cell(numDays*2+1,1); this_idx1 = []; this_idx2 = [];  
        % find phase indices
        if length(phaseType1)==2
            this_idx1 = thisDat.(phaseType1(1))&thisDat.(phaseType1(2));
        else
            this_idx1 = thisDat.(phaseType1(1));    % baseline D phase
        end
        this_idx2 = thisDat.(phaseType2(1))&thisDat.(phaseType2(2));
        if sum(this_idx2)==0;   error("After-reversal data is empty"); end
        % find respective session day indices
        for d = 1:numDays
            days_idx{d+numDays} = this_idx2&(thisDat.SessionDay==d);                % after reversals
            days_idx{d} = [days_idx{d+numDays}(1+numDays:end),false(1,numDays)];    % pick indices for n days earlier 
            %disp(find(days_idx{d})); disp(find(days_idx{d+numDays}));
        end
        days_idx{end} = this_idx2&(thisDat.SessionDay==d+1);   % disp(find(days_idx{end}));
        
        % divide eGFP (hm4Di) group into two: CNO first or VEH first
        CNO_first_idx = false(numel(animalLabels),1);
        VEH_first_idx = false(numel(animalLabels),1);
        for s = 1:numel(animalLabels)
            thisSub_idx = (thisDat.subjectLabel==animalLabels(s));
            if sum(thisDat.R1(thisSub_idx)==thisDat.CNO(thisSub_idx))==sum(thisSub_idx)
                if sum(thisDat.R2(thisSub_idx)==thisDat.VEH(thisSub_idx))~=sum(thisSub_idx); error(animalLabels(s)+": Indexing error!"); end
                CNO_first_idx(s) = true;    % if R1 & CNO indices numbers match, assign to R1-CNO group
            elseif sum(thisDat.R1(thisSub_idx)==thisDat.VEH(thisSub_idx))==sum(thisSub_idx)
                if sum(thisDat.R2(thisSub_idx)==thisDat.CNO(thisSub_idx))~=sum(thisSub_idx); error(animalLabels(s)+": Indexing error!"); end
                VEH_first_idx(s) = true;    % if R1 & VEH indices numbers match, assing to R1-VEH group
            end
        end
        if sum(CNO_first_idx+VEH_first_idx)~=length(CNO_first_idx); error("Group numbers does not add up"); end
        CNO12_groups = {animalLabels(CNO_first_idx), animalLabels(VEH_first_idx)};
        
        % plot data for each (two) group: CNO first & VEH first
        for S = 1:2
            thisGroup = CNO12_groups{S};
            NN{1}(S) = numel(thisGroup);    % track subject N for legend
            if contains(group_label,"hm4Di"); groupColor = lesion_colors{S}; end
            
            % set x-axis: session day 1 through 5
            xx = [1:numDays*2+1];
            yy1 = {};   % before rev
            yy2 = {};   % after rev
            clear SEM y E
            for s = 1:numel(thisGroup)
                thisSub_idx = (thisDat.subjectLabel==thisGroup(s));     %disp(sum(thisSub_idx));
                
                yy1{s} = nan(size(1:numDays)); 
                yy2{s} = nan(size(numDays+1:numDays*2+1));
                
                % loop through each set day: before rev
                for x = 1:numDays
                    % if s==1; Y{x} = []; end       % Method2: for averaging from entire subjects data combined
                    % set subject session day indices
                    thisDayIdx = thisSub_idx&days_idx{x};
                    yy1{s}(x) = nanmean(DataToPlot(thisDayIdx));
%                     Y{x} = [Y{x}; DataToPlot(thisIdx)]; % bin for Method 2
                end
                % after rev
                t = 0;
                for x = numDays+1:numDays*2+1
                    t = t + 1;
                    thisDayIdx = thisSub_idx&days_idx{x};
                    yy2{s}(t) = nanmean(DataToPlot(thisDayIdx));
                end
            end
            y1 = nanmean(cell2mat(yy1'),1);   % METHOD 1: Mean of subject means
            y2 = nanmean(cell2mat(yy2'),1);
            SEM1 = std(cell2mat(yy1'),0,1,'omitnan')./sqrt(sum(~isnan(cell2mat(yy1')),1));  % divide by # of non-NaN valued subjects
            SEM2 = std(cell2mat(yy2'),0,1,'omitnan')./sqrt(sum(~isnan(cell2mat(yy2')),1));
            
            if isempty(thisGroup)
                y1 = nan(size(y1)); SEM1 = y1;  % if no data, fill with NaN (for legend to work)
                y2 = nan(size(y2)); SEM2 = y2; 
            end          
            if length(y1)+length(y2)~=length(xx);   error("check y length: "+length(y1)+", "+length(y2)); end
            if length(SEM1)+length(SEM2)~=length(xx);   error("check y length: "+length(SEM1),", "+length(SEM1)); end
            
            if shadedError
                E = shadedErrorBar(xx,[y1,y2],[SEM1,SEM2],'lineProps',{'LineWidth', 1.5, 'Color',groupColor,'lineStyle',group_linetypes{S}}); hold on;
                % 1. Shaded Error bar
                % befor rev
%                 E1 = shadedErrorBar([xx(1:numDays),xx(numDays+1)],[y1,y2(1)],[SEM1,SEM2(1)],'lineProps',{'LineWidth', 1.5, 'Color',groupColor,'lineStyle',group_linetypes{S}}); hold on;
%                 if (S==1&&strcmp(phaseType1(1),"R1"))||(S==2&&strcmp(phaseType1(1),"R2"))
%                     E1.patch.FaceColor = 'm';                   
%                 end
                % after rev
%                 E2 = shadedErrorBar(xx(numDays+1:end),[y2],[SEM2],'lineProps',{'LineWidth', 1.5, 'Color',groupColor,'lineStyle',group_linetypes{S}}); hold on;
%                 E2.mainLine.HandleVisibility = 'off';
%                 if (S==1&&strcmp(phaseType2(1),"R1"))||(S==2&&strcmp(phaseType2(1),"R2"))
%                     E2.patch.FaceColor = 'm';        
%                 end
            else
                % 2. Normal error bar
                E0 = errorbar(xx,[y1,y2],[SEM1,SEM2],'Color',groupColor,'Color',groupColor,'LineWidth',1.5); hold on;
                % Markers for virus group
                if ~contains(group_label,"eGFP") 
                    % befor rev
                    if S==1&&strcmp(phaseType1(1),"R1")
                        P1 = plot(xx(1:numDays),y1,'LineStyle','none','Marker',markersSet{S},'MarkerSize',6,'HandleVisibility','off');
                        P1.MarkerFaceColor = groupColor;
                        P1.MarkerEdgeColor = 'm';
                    end
                    if S==2&&strcmp(phaseType1(1),"R2")
                        P1 = plot(xx(1:numDays),y1,'LineStyle','none','Marker',markersSet{S},'MarkerSize',6,'HandleVisibility','off');
                        P1.MarkerFaceColor = groupColor;
                        P1.MarkerEdgeColor = 'm';
                    end
                    % after rev
                    if S==1&&strcmp(phaseType2(1),"R1")
                        P2 = plot(xx(numDays+1:end),y2,'LineStyle','none','Marker',markersSet{S},'MarkerSize',6,'HandleVisibility','off');
                        P2.MarkerFaceColor = groupColor;
                        P2.MarkerEdgeColor = 'm';
                    end
                    if S==2&&strcmp(phaseType2(1),"R2")
                        P2 = plot(xx(numDays+1:end),y2,'LineStyle','none','Marker',markersSet{S},'MarkerSize',6,'HandleVisibility','off');
                        P2.MarkerFaceColor = groupColor;
                        P2.MarkerEdgeColor = 'm';
                    end
                end
            end
        end
        
        ax = gca;
        ax.XLim = [1 numDays*2+1];
        ax.XTick = xx;
        ax.XTickLabel = [-numDays:1:numDays];
        ax.YLabel.String = parameterLabels(i);
        if i==1
            title([lesion_region,'   (N = ',num2str(N),')']);
            L = legend({['R1-CNO (n=',num2str(NN{1}(1)),')'],['R2-CNO (n=',num2str(NN{1}(2)),')']});
            L.Box = 'off';
        end
        if i==2
            L = legend({['R1-CNO (n=',num2str(NN{1}(1)),')'],['R2-CNO (n=',num2str(NN{1}(2)),')']});
            L.Box = 'off';
            text(.25,1,phaseType1(1),'units','normalized','FontWeight','bold','FontSize',20,'HorizontalAlignment','right');
            text(.75,1,phaseType2(1),'units','normalized','FontWeight','bold','FontSize',20,'HorizontalAlignment','right');
        end
        if i>=7; ax.XLabel.String = "Days"; end
    
        ax.Box = 'off';
        set(ax,'FontName','Helvetica','FontSize',15,'FontWeight','normal','LineWidth',3)
        if i<=5; ax.YLim = [0 1]; plot(ax.XLim, [0.5, 0.5],':k','HandleVisibility','off'); end  % chance level line
        plot([numDays+1,numDays+1],ax.YLim,'--k','HandleVisibility','off'); % reversal line
        
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
%     parameterLabels = ["ERDS_{Loc}", "ERDS_+", "ERDS_-"; "EODS", "EODS_B", "EODS_W"; 'ERODS', "ERODS_{B+}", "ERODS_{W+}"; 'BLANK', "ERODS_{B-}", "ERODS_{W-}"];

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
        animalLabels = unique(thisDat.subjectLabel);   
        N = numel(animalLabels);                % track subject N for title
        DataToPlot = thisDat.(parameters(jj,kk));
        
        numDays = 3;    % last 3 days before rev, first 3+1 days after rev
        days_idx = cell(numDays*2+1,1); this_idx1 = []; this_idx2 = [];  
        % find phase indices
        if length(phaseType1)==2
            this_idx1 = thisDat.(phaseType1(1))&thisDat.(phaseType1(2));
        else
            this_idx1 = thisDat.(phaseType1(1));    % baseline D phase
        end
        this_idx2 = thisDat.(phaseType2(1))&thisDat.(phaseType2(2));
        if sum(this_idx2)==0;   error("After-reversal data is empty"); end
        % find respective session day indices
        for d = 1:numDays
            days_idx{d+numDays} = this_idx2&(thisDat.SessionDay==d);                % after reversals
            days_idx{d} = [days_idx{d+numDays}(1+numDays:end),false(1,numDays)];    % pick indices for n days earlier 
            %disp(find(days_idx{d})); disp(find(days_idx{d+numDays}));
        end
        days_idx{end} = this_idx2&(thisDat.SessionDay==d+1);   % disp(find(days_idx{end}));
        
        % divide eGFP (hm4Di) group into two: CNO first or VEH first
        CNO_first_idx = false(numel(animalLabels),1);
        VEH_first_idx = false(numel(animalLabels),1);
        for s = 1:numel(animalLabels)
            thisSub_idx = (thisDat.subjectLabel==animalLabels(s));
            if sum(thisDat.R1(thisSub_idx)==thisDat.CNO(thisSub_idx))==sum(thisSub_idx)
                if sum(thisDat.R2(thisSub_idx)==thisDat.VEH(thisSub_idx))~=sum(thisSub_idx); error(animalLabels(s)+": Indexing error!"); end
                CNO_first_idx(s) = true;    % if R1 & CNO indices numbers match, assign to R1-CNO group
            elseif sum(thisDat.R1(thisSub_idx)==thisDat.VEH(thisSub_idx))==sum(thisSub_idx)
                if sum(thisDat.R2(thisSub_idx)==thisDat.CNO(thisSub_idx))~=sum(thisSub_idx); error(animalLabels(s)+": Indexing error!"); end
                VEH_first_idx(s) = true;    % if R1 & VEH indices numbers match, assing to R1-VEH group
            end
        end
        if sum(CNO_first_idx+VEH_first_idx)~=length(CNO_first_idx); error("Group numbers does not add up"); end
        CNO12_groups = {animalLabels(CNO_first_idx), animalLabels(VEH_first_idx)};
        
        % plot data for each (two) group: CNO first & VEH first
        for S = 1:2
            thisGroup = CNO12_groups{S};
            NN{1}(S) = numel(thisGroup);    % track subject N for legend
            % if contains(group_label,"hm4Di"); groupColor = lesion_colors{S}; end
            groupColor = lesion_colors{S};
            
            % set x-axis: session day 1 through 5
            xx = [1:numDays*2+1];
            yy1 = {};   % before rev
            yy2 = {};   % after rev
            clear SEM y E
            for s = 1:numel(thisGroup)
                thisSub_idx = (thisDat.subjectLabel==thisGroup(s));     %disp(sum(thisSub_idx));
                
                yy1{s} = nan(size(1:numDays)); 
                yy2{s} = nan(size(numDays+1:numDays*2+1));
                
                % loop through each set day: before rev
                for x = 1:numDays
                    % if s==1; Y{x} = []; end       % Method2: for averaging from entire subjects data combined
                    % set subject session day indices
                    thisDayIdx = thisSub_idx&days_idx{x};
                    yy1{s}(x) = nanmean(DataToPlot(thisDayIdx));
%                     Y{x} = [Y{x}; DataToPlot(thisIdx)]; % bin for Method 2
                end
                % after rev
                t = 0;
                for x = numDays+1:numDays*2+1
                    t = t + 1;
                    thisDayIdx = thisSub_idx&days_idx{x};
                    yy2{s}(t) = nanmean(DataToPlot(thisDayIdx));
                end
            end
            y1 = nanmean(cell2mat(yy1'),1);   % METHOD 1: Mean of subject means
            y2 = nanmean(cell2mat(yy2'),1);
            SEM1 = std(cell2mat(yy1'),0,1,'omitnan')./sqrt(sum(~isnan(cell2mat(yy1')),1));  % divide by # of non-NaN valued subjects
            SEM2 = std(cell2mat(yy2'),0,1,'omitnan')./sqrt(sum(~isnan(cell2mat(yy2')),1));
            
            if isempty(thisGroup)
                y1 = nan(size(y1)); SEM1 = y1;  % if no data, fill with NaN (for legend to work)
                y2 = nan(size(y2)); SEM2 = y2; 
            end          
            if length(y1)+length(y2)~=length(xx);   error("check y length: "+length(y1)+", "+length(y2)); end
            if length(SEM1)+length(SEM2)~=length(xx);   error("check y length: "+length(SEM1),", "+length(SEM1)); end
            
            if shadedError
                % 1. Shaded Error bar
                E = shadedErrorBar(xx,[y1,y2],[SEM1,SEM2],'lineProps',{'LineWidth', 1.5, 'Color',groupColor,'lineStyle',group_linetypes{S}}); hold on;
            else
                % 2. Normal error bar
                E0 = errorbar(xx,[y1,y2],[SEM1,SEM2],'Color',groupColor,'Color',groupColor,'LineWidth',1.5); hold on;
                % Markers for virus group
                if ~contains(group_label,"eGFP") 
                    % befor rev
                    if S==1&&strcmp(phaseType1(1),"R1")
                        P1 = plot(xx(1:numDays),y1,'LineStyle','none','Marker',markersSet{S},'MarkerSize',6,'HandleVisibility','off');
                        P1.MarkerFaceColor = groupColor;
                        P1.MarkerEdgeColor = 'm';
                    end
                    if S==2&&strcmp(phaseType1(1),"R2")
                        P1 = plot(xx(1:numDays),y1,'LineStyle','none','Marker',markersSet{S},'MarkerSize',6,'HandleVisibility','off');
                        P1.MarkerFaceColor = groupColor;
                        P1.MarkerEdgeColor = 'm';
                    end
                    % after rev
                    if S==1&&strcmp(phaseType2(1),"R1")
                        P2 = plot(xx(numDays+1:end),y2,'LineStyle','none','Marker',markersSet{S},'MarkerSize',6,'HandleVisibility','off');
                        P2.MarkerFaceColor = groupColor;
                        P2.MarkerEdgeColor = 'm';
                    end
                    if S==2&&strcmp(phaseType2(1),"R2")
                        P2 = plot(xx(numDays+1:end),y2,'LineStyle','none','Marker',markersSet{S},'MarkerSize',6,'HandleVisibility','off');
                        P2.MarkerFaceColor = groupColor;
                        P2.MarkerEdgeColor = 'm';
                    end
                end
            end
        end
        
        ax = gca;
        ax.XLim = [1 numDays*2+1];
        ax.XTick = xx;
        ax.XTickLabel = [-numDays:1:numDays];
        ax.YLabel.String = parameterLabels(jj,kk);
        if i>=7; ax.XLabel.String = "Days"; end
    
        ax.Box = 'off';
        set(ax,'FontName','Helvetica','FontSize',15,'FontWeight','normal','LineWidth',3)
%         if i<=5; ax.YLim = [0 1]; plot(ax.XLim, [0.5, 0.5],':k','HandleVisibility','off'); end
        if i==1
            ylim([.5 1.0]);
%             title([group_label(1:3),'   (N = ',num2str(N),')']);
            L = legend({['R1-CNO (n=',num2str(NN{1}(1)),')'],['R2-CNO (n=',num2str(NN{1}(2)),')']});
            L.Box = 'off';
            text(.25,1,phaseType1(1),'units','normalized','FontWeight','bold','FontSize',18,'HorizontalAlignment','right');
            text(.75,1,phaseType2(1),'units','normalized','FontWeight','bold','FontSize',18,'HorizontalAlignment','right');
        end
        plot([numDays+1,numDays+1],ax.YLim,'--k','HandleVisibility','off'); % reversal line        
    
    end
end
sgtitle(taskType+" task: "+ phaseType1(1));
if length(phaseType1)==2
    sgtitle([group_label(1:3)+" (N = "+num2str(N)+"): "+taskType+" task: "+phaseType1(1)+" & "+phaseType1(2)]);
end
toc
end

end
function plot_metrics_distribution_AO_SO(all_output, task_label, phase_or_schedule_label, addstring, show_plots, plot_MI_cdf)

    plot_Behavioral_met = show_plots(1);
    plot_Cond_ent = show_plots(2);
    plot_Mutual_info = show_plots(3);
    plot_supp_probs = show_plots(4);
    
    task_output = all_output.(task_label);
    
    lesion_subsets = ["ACC_eGFP", "ACC_hm4Di", ...
                 "BLA_eGFP", "BLA_hm4Di", ...
                 "OFC_eGFP", "OFC_hm4Di"];  
    % colors: ACC, BLA, OFC         
    lesion_colors = {[0 0.4470 0.7410],[0 0.4470 0.7410],...
                [0.6350 0.0780 0.1840],[0.6350 0.0780 0.1840],...
                [0.9290 0.6940 0.1250],[0.9290 0.6940 0.1250]};
    linetypes = {':','-',':','-',':','-'};          % lesion = straight line; control = dotted lines

    if contains(addstring,"eGFP only")
        id = [1,3,5];
    elseif contains(addstring,"hm4Di only")
        id = [2,4,6];
    else 
        id = 1:numel(lesion_subsets);
    end
    lesion_subsets = lesion_subsets(id);
    lesion_colors = lesion_colors(id);
    linetypes = linetypes(id);
    
    numGroups = numel(lesion_subsets);
    
    alpha_val = 0.05;
    
    %% select phase/schedule label

    plot_data_idx = cell(1,numGroups);
    for k = 1:numGroups
        
        if length(phase_or_schedule_label)==1
            phase_and_schedule_idx = task_output.(lesion_subsets(k)).(phase_or_schedule_label);
        elseif length(phase_or_schedule_label)==2
            phase_and_schedule_idx = task_output.(lesion_subsets(k)).(phase_or_schedule_label(1))&task_output.(lesion_subsets(k)).(phase_or_schedule_label(2));
        else
            error("Select two labels only");
        end
        
        plot_data_idx{k} = logical(phase_and_schedule_idx);
    end
    
    %% 1. Behavioral metrics
if plot_Behavioral_met
    
    parameters = ["pwin","pbetter","pstay", "winstay","loseswitch","matching_measure", "RI_BW","RI_B","RI_W"];
    parameterLabels = ["prob(Win)","prob(Better)","prob(Stay)", "Win-Stay","Lose-Switch","dev.matching", "RI_{BW}","RI_B","RI_W"];

    figure;
    set(gcf,'Units','normalized','Position',[0,0,1,1],'color','w');

    for i = 1:numel(parameters)
        dataTemp = cell(numGroups,1);
        for k = 1:numGroups   
            dataTemp{k} = task_output.(lesion_subsets(k)).(parameters(i))(plot_data_idx{k});
        end
        tempMean = [];
        for cntDD = 1:numGroups
            tempMean(cntDD) = nanmean(dataTemp{cntDD});
        end

        subplot(3,3,i);
        ax = gca;   %axYLim = ax.YLim; %disp(ax.YLim);
        set(ax,'FontName','Helvetica','FontSize',14,'FontWeight','normal','LineWidth',3);
        set(ax, 'tickdir', 'out');
        ax.Box = 'off'; 
        ylabel(parameterLabels(i));
        hold on
        
        xmax = ceil(max(cell2mat(dataTemp))*10)/10;
        if xmax>0.5;    xmax = 1;   end
        xmin = 0;        
        if i>=6 
            xmin = min(cell2mat(dataTemp));
            xmin = floor(xmin*10)/10;
        end
        x_pd = linspace(xmin,xmax);
        for p = 1:numGroups
            pd{p} = fitdist(dataTemp{p},'Kernel');
            y_pd = pdf(pd{p}, x_pd);
            plot(x_pd, y_pd,'LineWidth',2,'Color', lesion_colors{p},'LineStyle',linetypes{p}); hold on;
        end

        for p = 1:numGroups
           if strcmp(linetypes{p},'-'); meanlinetype = '--';    else;   meanlinetype = linetypes{p};    end
           plot([tempMean(p) tempMean(p)],ax.YLim, 'Color',lesion_colors{p},'LineWidth',2,'LineStyle',meanlinetype); hold on  
        end
        
        % significance test from chance levels (zero or 0.5)
        if strcmp(parameters(i),"pbetter")||strcmp(parameters(i),"pstay")||strcmp(parameters(i),"winstay")||strcmp(parameters(i),"loseswitch")
            chance_lvl = 0.5;
            perf_test = 1;
        elseif contains(parameters(i),"RI")||strcmp(parameters(i),"matching_measure")
            chance_lvl = 0;
            perf_test = 1;
        else
            perf_test = 0;
        end
        if perf_test
            for p = 1:numGroups
                [h, ~] = ttest(dataTemp{p},chance_lvl,'alpha',alpha_val);
                if h
                   text(tempMean(p),ax.YLim(2),["*"],'Color', lesion_colors{p},'FontSize',20,'VerticalAlignment','bottom','HorizontalAlignment','center');
                end
            end
        end
        
        % pair-wise significance test
%         disp('----------------');
%         disp(parameterLabels(i)+" - "+phase_or_schedule_label)
        groups = nchoosek([1:numGroups],2);
        for c = 1:size(groups,1)
            gA = groups(c,1);
            gB = groups(c,2);
            groupA_dat = dataTemp{gA};      
            groupB_dat = dataTemp{gB};
            cohensD_vals(c) = my_cohensD(groupA_dat, groupB_dat);
            [t(c), p_vals(c)] = ttest2(groupA_dat, groupB_dat,'alpha',alpha_val);
%             disp("T-test with Cohen's D b.w. "+lesion_subsets{gA}+" & "+lesion_subsets{gB}+" groups: (D="+cohensD_vals(c)+", p="+p_vals(c)+ ")"); 
            [k(c), p(c), K] = kstest2(groupA_dat, groupB_dat,'alpha',alpha_val);
%             disp("Kolmogorov-Smirnov Test b.w. "+lesion_subsets{gA}+" & "+lesion_subsets{gB}+" groups: (h = "+k(c)+", D="+K+", p="+p(c)+ ")"); 
        end
        for m = 1:length(t)
            if t(m)
                text(1,1-(m-1)*.05,['$',num2str(groups(m,1)),'\ne$',num2str(groups(m,2))],'VerticalAlignment','top','HorizontalAlignment','right','Units','normalized','interpreter','latex'); 
            end
        end
        % all pairwise significant difference
        if sum(t)==numGroups; text(ax.XLim(2),ax.YLim(2),'T*','VerticalAlignment','bottom','HorizontalAlignment','right');   end
        if sum(k)==numGroups; text(ax.XLim(2),ax.YLim(2),'k*','VerticalAlignment','bottom');  end
        
        if i==1; title([task_label+" task, "+phase_or_schedule_label+" "+addstring]); end 
    end

end
    %% 2. Conditional Entropy
if plot_Cond_ent
    
    parameters = ["ERDS", "ERDS_win", "ERDS_lose"; "EODS", "EODS_better", "EODS_worse"; "ERODS", "ERODS_winbetter", "ERODS_winworse"; "BLANK", "ERODS_losebetter", "ERODS_loseworse"];
    parameterLabels = ["ERDS", "ERDS_+", "ERDS_-"; "EODS", "EODS_B", "EODS_W"; 'ERODS', "ERODS_{B+}", "ERODS_{W+}"; 'BLANK', "ERODS_{B-}", "ERODS_{W-}"];
            
    figure;
    set(gcf,'Units','normalized','Position',[0,0,1,1],'Color','w');
    
    for jj = 1:size(parameters,1)
        for kk = 1:size(parameters,2)
            i = kk + (jj-1)*3;
   
            if ~(jj==4 && kk == 1)  
                dataTemp = cell(numGroups,1);
                for k = 1:numGroups
                    dataTemp{k} = task_output.(lesion_subsets(k)).(parameters(jj,kk))(plot_data_idx{k});
                end
                
                tempMean = [];
                tempMed = [];
                for cntDD = 1:numGroups
                    tempMean(cntDD) = nanmean(dataTemp{cntDD});
                    tempMed(cntDD) = nanmedian(dataTemp{cntDD});
                end

                subplot(4,3,i);
                ax = gca;   %disp(ax.YLim);              
                ylabel(parameterLabels(jj,kk));
                set(ax,'FontName','Helvetica','FontSize',15,'FontWeight','normal','LineWidth',3);
                       
                % fix: skip empty metrics (100/0 schedule)
                if sum(isnan(dataTemp{1}))==length(dataTemp{1});   continue;   end
                if sum(isnan(dataTemp{2}))==length(dataTemp{2});   continue;   end
                if sum(isnan(dataTemp{3}))==length(dataTemp{3});   continue;   end
                
                xmax = ceil(max(cell2mat(dataTemp))*10)/10;
                for p = 1:numGroups
                    pd{p} = fitdist(dataTemp{p},'Kernel');
                    x_pd = linspace(0,xmax);   y_pd = pdf(pd{p}, x_pd);
                     plot(x_pd, y_pd,'LineWidth',2,'Color', lesion_colors{p}); hold on;
                end
                set(ax, 'tickdir', 'out');
                ax.Box = 'off';  
                ylabel(parameterLabels(jj,kk));
                set(ax,'FontName','Helvetica','FontSize',15,'FontWeight','normal','LineWidth',3);
                for p = 1:numGroups
                    if strcmp(linetypes{p},'-'); meanlinetype = '--';    else;   meanlinetype = linetypes{p};    end
                    plot([tempMean(p) tempMean(p)],ax.YLim,'Color',lesion_colors{p},'LineWidth',2,'LineStyle',meanlinetype); hold on  
                    scatter([tempMed(p)],[0],'d','filled','MarkerFaceColor',lesion_colors{p});
                end
                if i<8
                    ax.XLim = [0, 1]; 
                else
                    ax.XLim = [0, 0.65]; 
                end
                
                % pair-wise significance test
%                 disp('----------------');
%                 disp(parameterLabels(jj,kk)+" - "+phase_or_schedule_label)
                groups = nchoosek([1:numGroups],2);
                for c = 1:size(groups,1)
                    gA = groups(c,1);
                    gB = groups(c,2);
                    groupA_dat = dataTemp{gA};      
                    groupB_dat = dataTemp{gB};
                    cohensD_vals(c) = my_cohensD(groupA_dat, groupB_dat);
                    [t(c), p_vals(c)] = ttest2(groupA_dat, groupB_dat,'alpha',alpha_val);
        %             disp("T-test with Cohen's D b.w. "+lesion_subsets{gA}+" & "+lesion_subsets{gB}+" groups: (D="+cohensD_vals(c)+", p="+p_vals(c)+ ")"); 
                    [k(c), p(c), K] = kstest2(groupA_dat, groupB_dat,'alpha',alpha_val);
        %             disp("Kolmogorov-Smirnov Test b.w. "+lesion_subsets{gA}+" & "+lesion_subsets{gB}+" groups: (h = "+k(c)+", D="+K+", p="+p(c)+ ")"); 
                end
                for m = 1:length(t)
                    if t(m)
                        text(1,1-(m-1)*.05,['$',num2str(groups(m,1)),'\ne$',num2str(groups(m,2))],'VerticalAlignment','top','HorizontalAlignment','right','Units','normalized','interpreter','latex'); 
                    end
                end
                % all pairwise significant difference
                if sum(t)==numGroups; text(ax.XLim(2),ax.YLim(2),'T*','VerticalAlignment','bottom','HorizontalAlignment','right');   end
                if sum(k)==numGroups; text(ax.XLim(2),ax.YLim(2),'k*','VerticalAlignment','bottom');  end
                                 
                if i==1; title([task_label+" task, "+phase_or_schedule_label+" "+addstring]); end 
                
            end
        end
    end    
        fff = gcf;
        mysub = fff.Children(5);
        mysub.Position = mysub.Position - [0 0.12 0 0];
%}
end

    %% 3. Mutual Information
if plot_Mutual_info

    if contains(addstring,'What')||contains(addstring,'what')
        parameters = ["MIRS_stim", "MIRS_stim_win", "MIRS_stim_lose"; "MIOS_stim", "MIOS_stim_better", "MIOS_stim_worse"; "MIROS_stim", "MIROS_stim_winbetter", "MIROS_stim_winworse"; "BLANK", "MIROS_stim_losebetter", "MIROS_stim_loseworse"];
        parameterLabels = ["MIRS", "MIRS_+", "MIRS_-"; "MIOS", "MIOS_B", "MIOS_W"; "MIROS", "MIROS_{B+}", "MIROS_{W+}"; 'BLANK', "MIROS_{B-}", "MIROS_{W-}"];
    elseif contains(addstring,'Where')||contains(addstring,'where')
        parameters = ["MIRS_side", "MIRS_side_win", "MIRS_side_lose"; "MIOS_side", "MIOS_side_better", "MIOS_side_worse"; "MIROS_side", "MIROS_side_winbetter", "MIROS_side_winworse"; "BLANK", "MIROS_side_losebetter", "MIROS_side_loseworse"];
        parameterLabels = ["MIRS", "MIRS_+", "MIRS_-"; "MIOS", "MIOS_B", "MIOS_W"; "MIROS", "MIROS_{B+}", "MIROS_{W+}"; 'BLANK', "MIROS_{B-}", "MIROS_{W-}"];
    end
    
    figure('Units','normalized','Position',[0,0,1,1]);
    
    for jj = 1:size(parameters,1)           % row
        for kk = 1:size(parameters,2)       % column
            i = kk + (jj-1)*3;
            if ~(jj==4 && kk == 1)  
                dataTemp = cell(numGroups,1);
                for k = 1:numGroups
                    dataTemp{k} = task_output.(lesion_subsets(k)).(parameters(jj,kk))(plot_data_idx{k});
                end
%                 dataTemp = {all_output.control.(parameters(jj,kk))(plot_data_idx{1}); all_output.amygdala.(parameters(jj,kk))(plot_data_idx{2}); all_output.VS.(parameters(jj,kk))(plot_data_idx{3})};
                tempMean = [];
                for cntDD = 1:numGroups
                    tempMean(cntDD) = nanmean(dataTemp{cntDD});
                end
                
                subplot(4,3,i);
                ax = gca; 
                ylabel(parameterLabels(jj,kk));
                set(ax,'FontName','Helvetica','FontSize',15,'FontWeight','normal','LineWidth',3);
                
                % for deterministic task, there's no winworse & losebetter
                if strcmp(phase_or_schedule_label,"prob1000")
                   if strcmp(parameters(jj,kk),"MIROS_winworse")||strcmp(parameters(jj,kk),"MIROS_losebetter")
                      continue; 
                   end
                end
                
                
            if plot_MI_cdf
                for p = 1:numGroups
                    pd{p} = fitdist(dataTemp{p},'Kernel');
                    [fe,xe] = ecdf(dataTemp{p});
                    stairs(xe, fe, lesion_colors{p},'LineWidth',0.5); hold on;
                    xmin = min(dataTemp{p});
                    xmax = max(dataTemp{p});
                    x = linspace(xmin,xmax);
                    y = cdf(pd{p},x);
                     plot(x_pd, y_pd,'LineWidth',2,'Color', lesion_colors{p}); hold on;
                end
            else
                xmax = max(cell2mat(dataTemp));
%                 if xmax>0.5;    xmax = 1;   end
                for p = 1:numel(plot_data_idx)
                    pd{p} = fitdist(dataTemp{p},'Kernel');
                    x_pd = linspace(0,xmax);   y_pd = pdf(pd{p}, x_pd);
                     plot(x_pd, y_pd,'LineWidth',2,'Color', lesion_colors{p}); hold on;
                end
            end
            set(ax,'FontName','Helvetica','FontSize',15,'FontWeight','normal','LineWidth',3);
            set(ax, 'tickdir', 'out');
            ax.Box = 'off';
            if contains(addstring,'loc')||contains(addstring,'Loc')
                if kk==1 
                    if jj==3
                        ax.XLim = [0, 0.3];
                    else
                        ax.XLim = [0, 0.2]; 
                    end
                else
                    ax.XLim = [0, 0.1]; 
                end
            else
                if kk==1 
                    ax.XLim = [0, 0.3]; 
                else
                    ax.XLim = [0, 0.2]; 
                end
            end
            
            ylabel(parameterLabels(jj,kk));
            for p = 1:numGroups
               plot([tempMean(p) tempMean(p)],ax.YLim,'Color',lesion_colors{p},'LineWidth',2,'LineStyle',':'); hold on  
            end    
                
            % pair-wise significance test
            disp('----------------');
            disp(parameterLabels(jj,kk)+" - "+phase_or_schedule_label)
            if numGroups==3
               groups = {[1,2],[2,3],[1,3]};
            elseif numGroups==4
               groups = {[1,3],[3,4],[1,4]};
            end
            for c = 1:numel(groups)
                gA = groups{c}(1);
                gB = groups{c}(2);
                groupA_dat = dataTemp{gA};      
                groupB_dat = dataTemp{gB};
                cohensD_vals(c) = my_cohensD(groupA_dat, groupB_dat);
                [t(c), p_vals(c)] = ttest2(groupA_dat, groupB_dat,'alpha',0.001);
                disp("T-test with Cohen's D b.w. "+lesion_subsets{gA}+" & "+lesion_subsets{gB}+" groups: (D="+cohensD_vals(c)+", p="+p_vals(c)+ ")"); 
                [k(c), p(c), K] = kstest2(groupA_dat, groupB_dat,'alpha',0.001);
                disp("Kolmogorov-Smirnov Test b.w. "+lesion_subsets{gA}+" & "+lesion_subsets{gB}+" groups: (h = "+k(c)+", D="+K+", p="+p(c)+ ")"); 
            end
            if sum(t)==3; text(ax.XLim(2),ax.YLim(2),'T*','VerticalAlignment','bottom','HorizontalAlignment','right');   end
            if sum(k)==3; text(ax.XLim(2),ax.YLim(2),'k*','VerticalAlignment','bottom');  end
                

                if i==1; title([phase_or_schedule_label +" "+addstring]); end 
            end
        end
    end    
    fff = gcf;
    mysub = fff.Children(5);
    mysub.Position = mysub.Position - [0 0.12 0 0];
        
    %}
        
end

%% 4. Supplementary conditional probabilities

if plot_supp_probs
    parameters = ["betterstay", "winbetter", "winworse", "losebetter", "loseworse";"worseswitch", "winstaybetter", "winstayworse", "loseswitchbetter", "loseswitchworse"];
    parameterLabels = ["BetterStay", "WinBetter", "WinWorse", "LoseBetter", "LoseWorse";"WorseSwitch", "WinStayBetter", "WinStayWorse", "LoseSwitchBetter", "LoseSwitchWorse"];
    if contains(addstring,'loc')||contains(addstring,'Loc')
        parameters = ["winright","loseright","rightstay","winstayright","loseswitchright"; ...
                      "winleft","loseleft","leftstay","winstayleft","loseswitchleft"];
        parameterLabels = ["WinRight", "LoseRight", "RightStay", "WinStayRight", "LoseSwitchRight"; ...
                            "WinLeft", "LoseLeft", "LeftStay", "WinStayLeft", "LoseSwitchLeft"];
    end

    
    
    
    figure('Units','normalized','Position',[0 0 1.0 0.6]);
    for jj = 1:size(parameters,1)           % row
        for kk = 1:size(parameters,2)       % column
            i = kk + (jj-1)*size(parameters,2);
            dataTemp = cell(numGroups,1);
            for k = 1:numGroups
                dataTemp{k} = task_output.(lesion_subsets(k)).(parameters(jj,kk))(plot_data_idx{k});
            end
%             dataTemp = {all_output.control.(parameters(jj,kk))(plot_data_idx{1}); all_output.amygdala.(parameters(jj,kk))(plot_data_idx{2}); all_output.VS.(parameters(jj,kk))(plot_data_idx{3})};
            tempMean = [];
            for cntDD = 1:numGroups
                tempMean(cntDD) = nanmean(dataTemp{cntDD});
            end
            
            subplot(2,5,i);
            ax = gca; 
            ylabel(parameterLabels(jj,kk));
            set(ax,'FontName','Helvetica','FontSize',15,'FontWeight','normal','LineWidth',3);
    
            xmax = ceil(max(cell2mat(dataTemp))*10)/10;
%                 xmax = 1;
%                 if i>=8; xmax = 0.65;    end
%                 xmax = max(cell2mat(dataTemp));
%                 if xmax>0.5;    xmax = 1;   end
            for p = 1:numGroups
                pd{p} = fitdist(dataTemp{p},'Kernel');
                x_pd = linspace(0,xmax);   y_pd = pdf(pd{p}, x_pd);
                plot(x_pd, y_pd, lesion_colors{p},'LineWidth',2); hold on;
            end
            set(ax, 'tickdir', 'out');
            ax.Box = 'off';  
            ylabel(parameterLabels(jj,kk));
            set(ax,'FontName','Helvetica','FontSize',15,'FontWeight','normal','LineWidth',3);
            
            for p = 1:numGroups
               plot([tempMean(p) tempMean(p)],ax.YLim,'Color',lesion_colors{p},'LineWidth',2,'LineStyle',':'); hold on  
            end
            if kk<=2;   ax.XLim = [0 0.6];  end

            % T-test from mean of 0.5
            if kk>=3
                disp('----------------'); 
                for p = 1:numGroups
                    [t, pval] = ttest(dataTemp{p}, 0.5, 'alpha',0.05/length(dataTemp{p}));
                    disp("One sample T-Test from mean of 0.5 in "+lesion_subsets{p}+" group: (h = "+t+", p = "+pval+")");
                    disp("> "+parameterLabels(jj,kk)+ " "+lesion_subsets{p}+" group mean = "+tempMean(p));
                    if t==1
                        text(tempMean(p), ax.YLim(2),'*','FontSize',30,'Color',lesion_colors{p},'VerticalAlignment','baseline','HorizontalAlignment','center');
                        text(ax.XLim(1),ax.YLim(2)*(9-p)/10,num2str(tempMean(p),3),'Color',lesion_colors{p},'HorizontalAlignment','left');
                    end
                end
            end
            if jj==1
                dataTem2p = cell(numGroups,1);
                for k = 1:numGroups
                    dataTemp{k} = task_output.(lesion_subsets(k)).(parameters(2,kk))(plot_data_idx{k});
                end
%                 dataTemp2 = {all_output.control.(parameters(2,kk))(plot_data_idx{1}); all_output.amygdala.(parameters(2,kk))(plot_data_idx{2}); all_output.VS.(parameters(2,kk))(plot_data_idx{3})};
                for p = 1:numGroups
                    [t2, pval2] = ttest2(dataTemp{p},dataTemp2{p}, 'alpha',0.05/length(dataTemp{p}));
                    if t2==1
                        text(tempMean(p), ax.YLim(1),'+','FontSize',40,'Color',lesion_colors{p},'VerticalAlignment','cap','HorizontalAlignment','center');
                    end
                end
            end
            
%             % pair-wise significance test
%             disp('----------------');
%             disp(parameterLabels(jj,kk)+" - "+schedule_label)
%             if numGroups==3
%                groups = {[1,2],[2,3],[1,3]};
%             elseif numGroups==4
%                groups = {[1,3],[3,4],[1,4]};
%             end
%             for c = 1:numel(groups)
%                 gA = groups{c}(1);
%                 gB = groups{c}(2);
%                 groupA_dat = dataTemp{gA};      
%                 groupB_dat = dataTemp{gB};
%                 cohensD_vals(c) = my_cohensD(groupA_dat, groupB_dat);
%                 [t(c), p_vals(c)] = ttest2(groupA_dat, groupB_dat,'alpha',0.001);
%                 disp("T-test with Cohen's D b.w. "+lesion_subsets{gA}+" & "+lesion_subsets{gB}+" groups: (D="+cohensD_vals(c)+", p="+p_vals(c)+ ")"); 
%                 [k(c), p(c), K] = kstest2(groupA_dat, groupB_dat,'alpha',0.001);
%                 disp("Kolmogorov-Smirnov Test b.w. "+lesion_subsets{gA}+" & "+lesion_subsets{gB}+" groups: (h = "+k(c)+", D="+K+", p="+p(c)+ ")"); 
%             end
%             if sum(t)==3; text(ax.XLim(2),ax.YLim(2),'T*','VerticalAlignment','bottom','HorizontalAlignment','right');   end
%             if sum(k)==3; text(ax.XLim(2),ax.YLim(2),'k*','VerticalAlignment','bottom');  end

            if i==1; title([phase_or_schedule_label +" "+addstring]); end 
            
        end
    end

end

end
% Analyze models for rats dataset
clearvars; clc; close all
if contains(pwd,'f004p51')
    cd('/Users/f004p51/Dropbox (Dartmouth College)/CCNL/Alicia_Data/Claudia_SO_AO_data')
else
    cd('C:\Users\wojh1\Dropbox (Dartmouth College)\CCNL\Alicia_Data\Claudia_SO_AO_data');
end

task_label = "AO";

% load dataset
include_NoChoiceTrials = 1;     % dataset including uncommitted choice trials
if include_NoChoiceTrials
    dataset_label = task_label+"_IncludeNoChoice";
else
    dataset_label = task_label;
end

% load models
[models] = initialize_RL_models(dataset_label, 1, 1);
winModel = models{2};   % RL1_decay

drugOrderSet = ["CNO_first","VEH_first"];
drugOrderLabel = ["CNO1","VEH1"];
drugOrderColors = {{'#FF6F6F','c'},{'#6666FF','g'}};

allPhases = {["D","prob1000"],["R1","prob1000"],["R2","prob1000"],["R1","prob9010"],["R2","prob9010"],["R1","prob8020"],["R2","prob8020"],["R1","prob7030"],["R2","prob7030"]};
allPhasesLegends = ["D","R1 (100/0)","R2 (100/0)","R3 (90/10)","R4 (90/10)","R5(80/20)","R6(80/20)","R7(70/30)","R8(70/30)"];
sexSet = ["Female","Male"];

%% Figure 3-1. Plot parameters by group*sex
clc
groups.set = {'eGFP','OFC_hm4Di','BLA_hm4Di'};
groups.labels = {'eGFP','vlOFC hm4Di','BLA hm4Di'};
groups.colors = {'#32CC32','#FF2E92','#FF2E92'};
% Row1: eGFP
% Row2: vlOFC hm4Di
% Row3: BLA hm4Di
% Column 1: Female R1
% Column 2: Male R1
% Column 3: Female R3
% Column 4: Male R3

PhasesToShow = [2,4];   % R1 & R3

allpars = cell2mat(winModel.Fit.fitpar');
%     allpars(:,2) = 1./allpars(:,2);   % temperature sigma instead of beta
%     thisModel.ub(2) = 1;
    
figure(31); clf
set(gcf,'Units','normalized','Position',[0,0,0.8,0.6], 'Color','w');
numfig = 0;
% loop through each lesion group
for g = 1:numel(groups.set)
    % R1/R3
    for i = 1:numel(PhasesToShow)
        phaseNum = PhasesToShow(i);
        phase_label = allPhases{phaseNum}(1);
        schedule_label = allPhases{phaseNum}(2);
        phase_idx = strcmp(winModel.SessInfo.phase,phase_label)&strcmp(winModel.SessInfo.schedules,schedule_label);
        disp(">>> "+i+". "+allPhasesLegends(phaseNum)+" = "+sum(phase_idx));

        % Female/Male
        for k = 1:numel(sexSet)
            if contains(groups.set{g},'hm4Di')
                group_idx = phase_idx&strcmp(winModel.SessInfo.condition, groups.set{g});
            elseif strcmp(groups.set{g},'eGFP')
                group_idx = phase_idx&contains(winModel.SessInfo.condition, groups.set{g});
            end
            disp("  ("+g+") "+groups.set{g}+": "+sum(group_idx)+" sessions");
            
            % one panel per each parameter
            for param = 1:numel(winModel.initpar)
                numfig = numfig + 1;
                if mod(numfig,12)==1
                    xPos = 0.0; 
                    x_offset = 0.0;
                else
                    if mod(numfig,3)==1
                        x_offset = 0.015;
                        if mod(numfig,6)==1
                            x_offset = 0.03;
                        end
                    else
                        x_offset = 0.0;
                    end
                end
                SP = subplot(numel(groups.set),numel(PhasesToShow)*2*numel(winModel.initpar),numfig);              
                xPos = xPos + 0.065 + x_offset;
                SP.Position(1) = xPos;
                SP.Position(3) = 0.06;
                
                % Compile data & plot bar graph for each group: CNO/VEH-first
                thisParams = struct;
                MeanVals = nan(1,2);
                SEMvals = nan(1,2);                
                for j = 1:numel(drugOrderSet)
                    sess_idx = group_idx&strcmp(winModel.SessInfo.drugOrder, drugOrderSet(j))&strcmp(winModel.SessInfo.sex, sexSet(k));
                    thisGroup_label = groups.set{g}+"_"+drugOrderSet(j)+"_"+sexSet(k);
                    if param==1; disp("     "+drugOrderSet(j)+"("+k+") "+sexSet(k)+": "+sum(sess_idx)+" params"); end
                    thisParams.(thisGroup_label) = allpars(sess_idx,param);
                    
                    MeanVals(j) = mean(thisParams.(thisGroup_label));
                    SEMvals(j) = std(thisParams.(thisGroup_label))/sqrt(numel(thisParams.(thisGroup_label)));
                   
                    hold on;
                    bar(j,MeanVals(j),'FaceColor',drugOrderColors{i}{j},'LineWidth',1);
                    errorbar(j,MeanVals(j),SEMvals(j),'k','LineStyle','none');
                    scatter(j*ones(numel(thisParams.(thisGroup_label)),1),thisParams.(thisGroup_label),15,'MarkerFaceColor',drugOrderColors{i}{j},'MarkerFaceAlpha',0.4,'MarkerEdgeColor','k','LineWidth',0.5);
                end
                xticks(1:2);
                xticklabels(drugOrderLabel);
                xtickangle(30);
                
                % sig-test between two groups
                all_groups = fieldnames(thisParams);
                [pval, Wilcoxon_r, asterisk] = Wilcoxon_ranksum_test(thisParams.(all_groups{1}), thisParams.(all_groups{2}), 1);

                if pval<.075
                    % yPosFact = 0.6;
                    yPosFact = ceil(max(MeanVals+SEMvals)*10)/10*2;
                    plot([1 2], repmat(yPosFact*winModel.ub(param),1,2),'-k', 'LineWidth',1.5);
                    plot([1 1], [yPosFact-.09, yPosFact].*winModel.ub(param),'-k', 'LineWidth',1.5);
                    plot([2 2], [yPosFact-.09 yPosFact].*winModel.ub(param),'-k', 'LineWidth',1.5);
                    if strcmp(asterisk,"n.s.")
                        asterisk = ["\it{p} = \rm"+num2str(pval,3)]; 
                    else
                        asterisk = [asterisk+"\it{p} = \rm"+num2str(pval,3)]; 
                    end
                    text(1.5, yPosFact*winModel.ub(param), [asterisk], 'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',10);
                end
                
                text(0.075, 1.0, winModel.plabels(param),'Units','normalized','VerticalAlignment','middle','FontSize',14);
                if param==1&&k==1
%                     text(0.15, 1.0, groups.labels{g},'Units','normalized','VerticalAlignment','top');
                    ylabel(groups.labels{g},'Color',groups.colors{g});
                end
                
                yLimFact = 1;
                ylim([winModel.lb(param), winModel.ub(param)*yLimFact]);
                set(gca,'FontName','Helvetica','FontSize',10,'FontWeight','normal','LineWidth',1, 'tickdir','out','Box','off');
                if g==1&&param==2
                    title(sexSet(k)+" "+allPhasesLegends(phaseNum),'FontSize',14,'Position',[1, 110],'HorizontalAlignment','center'); 
                end
            end
        end
    end
end

%% Alternative version1: Female + Male in the same panel
clc
groups.set = {'eGFP','OFC_hm4Di','BLA_hm4Di'};
groups.labels = {'eGFP','vlOFC hm4Di','BLA hm4Di'};
groups.colors = {'#32CC32','#FF2E92','#FF2E92'};
% Row1: eGFP
% Row2: vlOFC hm4Di
% Row3: BLA hm4Di
% Column 1: Female R1
% Column 2: Male R1
% Column 3: Female R3
% Column 4: Male R3

PhasesToShow = [2,4];   % R1 & R3

allpars = cell2mat(winModel.Fit.fitpar');
%     allpars(:,2) = 1./allpars(:,2);   % temperature sigma instead of beta?
%     thisModel.ub(2) = 1;
    
figure(310); clf
set(gcf,'Units','normalized','Position',[0,0,0.6,0.8], 'Color','w');
numfig = 0;
% loop through each lesion group
for g = 1:numel(groups.set)
    % R1/R3
    for i = 1:numel(PhasesToShow)
        phaseNum = PhasesToShow(i);
        phase_label = allPhases{phaseNum}(1);
        schedule_label = allPhases{phaseNum}(2);
        phase_idx = strcmp(winModel.SessInfo.phase,phase_label)&strcmp(winModel.SessInfo.schedules,schedule_label);
        disp(">>> "+i+". "+allPhasesLegends(phaseNum)+" = "+sum(phase_idx));
            
        % one panel per each parameter
        for param = 1:numel(winModel.initpar)
            disp("-----------"+winModel.plabels(param)+"-----------");
            numfig = numfig + 1;
            if mod(numfig,6)==1
                xPos = 0.0; 
                x_offset = 0.0;
            else
                if mod(numfig,3)==1
                    x_offset = 0.05;
                else
                    x_offset = 0.0;
                end
            end
            SP = subplot(numel(groups.set),numel(PhasesToShow)*numel(winModel.initpar),numfig);              
            xPos = xPos + 0.125 + x_offset;
            SP.Position(1) = xPos;
            SP.Position(3) = 0.11;
                
            % Compile data & plot bar graph for each group: Female/Male
            thisParams = struct;
            MeanVals = nan(1,4);
            SEMvals = nan(1,4);
            barNum = 0;
            for k = 1:numel(sexSet)
                if contains(groups.set{g},'hm4Di')
                    group_idx = phase_idx&strcmp(winModel.SessInfo.condition, groups.set{g});
                elseif strcmp(groups.set{g},'eGFP')
                    group_idx = phase_idx&contains(winModel.SessInfo.condition, groups.set{g});
                end
                disp("  ("+g+") "+groups.set{g}+": "+sum(group_idx)+" sessions");
                for j = 1:numel(drugOrderSet)
                    barNum = barNum + 1;
                    
                    sess_idx = group_idx&strcmp(winModel.SessInfo.drugOrder, drugOrderSet(j))&strcmp(winModel.SessInfo.sex, sexSet(k));
                    thisGroup_label = groups.set{g}+"_"+drugOrderSet(j)+"_"+sexSet(k);
                    if strcmp(groups.set{g},'eGFP')
                        sess_idx = group_idx&strcmp(winModel.SessInfo.sex, sexSet(k));
                        thisGroup_label = groups.set{g}+"_"+sexSet(k);
                    end
                    
                    if param==1; disp("     "+drugOrderSet(j)+"("+k+") "+sexSet(k)+": "+sum(sess_idx)+" params"); end
                    thisParams.(thisGroup_label) = allpars(sess_idx,param);
                    
                    MeanVals(barNum) = mean(thisParams.(thisGroup_label));
                    SEMvals(barNum) = std(thisParams.(thisGroup_label))/sqrt(numel(thisParams.(thisGroup_label)));
                    hold on;
                    B = bar(barNum,MeanVals(barNum),'FaceColor',drugOrderColors{i}{j},'LineWidth',1);
                    if barNum>2; B.HandleVisibility = 'off'; end
                    errorbar(barNum,MeanVals(barNum),SEMvals(barNum),'k','LineStyle','none','HandleVisibility','off');
                    scatter(barNum*ones(numel(thisParams.(thisGroup_label)),1),thisParams.(thisGroup_label),15,'MarkerFaceColor',drugOrderColors{i}{j},'MarkerFaceAlpha',0.4,'MarkerEdgeColor','k','LineWidth',0.5,'HandleVisibility','off');
                end
                xticks([1.5, 3.5]);
                xlim([0.25, 4.75]);
                xticklabels(sexSet); % xtickangle(30);
            end
            
            % sig-test between two groups
            all_groups = fieldnames(thisParams);
            for j = 1:1+(~strcmp(groups.set{g},'eGFP'))
                [pval, ~, asterisk, stats] = Wilcoxon_ranksum_test(thisParams.(all_groups{2*j-1}), thisParams.(all_groups{2*j}), 1);
                if 1%pval<.075
                    disp("Mean diff = "+ [mean(thisParams.(all_groups{2*j-1}))-mean(thisParams.(all_groups{2*j}))]+", pval "+pval);
                    yPosFact = 0.4;
                    plot([2*j-1 2*j], repmat(yPosFact*winModel.ub(param),1,2),'-k', 'LineWidth',1.5);
                    plot([2*j-1 2*j-1], [yPosFact-.05, yPosFact].*winModel.ub(param),'-k', 'LineWidth',1.5);
                    plot([2*j 2*j], [yPosFact-.05 yPosFact].*winModel.ub(param),'-k', 'LineWidth',1.5);
                    if strcmp(asterisk,"n.s.")
                        asterisk = ["\it{p} = \rm"+num2str(pval,3)]; 
                    else
                        asterisk = [asterisk+"\it{p} = \rm"+num2str(pval,3)]; 
                    end
                    text(2*j-0.5, yPosFact*winModel.ub(param), [asterisk], 'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',10);
                end
                if ~strcmp(groups.set{g},'eGFP')
                    % Female vs. Male
                    [pval, ~, asterisk] = Wilcoxon_ranksum_test(thisParams.(all_groups{j}), thisParams.(all_groups{j+2}), 1);
                    if pval<.075
                        yPosFact = 0.4 + j*0.1;
                        plot([j j+2], repmat(yPosFact*winModel.ub(param),1,2),'-k', 'LineWidth',1.5);
                        plot([j j], [yPosFact-.05, yPosFact].*winModel.ub(param),'-k', 'LineWidth',1.5);
                        plot([j+2 j+2], [yPosFact-.05 yPosFact].*winModel.ub(param),'-k', 'LineWidth',1.5);
                        if strcmp(asterisk,"n.s.")
                            asterisk = ["\it{p} = \rm"+num2str(pval,3)]; 
                        else
                            asterisk = [asterisk+"\it{p} = \rm"+num2str(pval,3)]; 
                        end
                        text(j+1, yPosFact*winModel.ub(param), [asterisk], 'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',10);
                    end
                end
            end   

            text(0.05, 1.05, winModel.plabels(param),'Units','normalized','VerticalAlignment','middle','FontSize',14);
            if param==1
                legend(["CNO1","VEH1"],'box','off','fontsize',8,'Location','northwest');
                ylabel(groups.labels{g},'Color',groups.colors{g});
                if g==1
                    text(-.5,1.1,char(i+'A'-1),'Units','normalized','FontWeight','bold','FontSize',20,'FontName','helvetica','VerticalAlignment','bottom');
                end
            end

            yLimFact = 1;
            ylim([winModel.lb(param), winModel.ub(param)*yLimFact]);
            yticks([winModel.lb(param):winModel.ub(param)/4:winModel.ub(param)]);
            SP.YTickLabel{2} = '';
            SP.YTickLabel{4} = '';
            set(gca,'FontName','Helvetica','FontSize',12,'FontWeight','normal','LineWidth',1, 'tickdir','out','Box','off');
            if g==1&&param==2
                title(allPhasesLegends(phaseNum),'FontSize',14,'Position',[2.5, 110],'HorizontalAlignment','center'); 
            end
            fprintf('\n');
        end
    end
end

%% Alternative version 2: 
clc
% Row1: vlOFC hm4Di vs. eGFP
% Row2: BLA hm4Di vs.eGFP
% Column 1: Female R1
% Column 2: Male R1
% Column 3: Female R3
% Column 4: Male R3

groups.set = {'OFC_hm4Di','BLA_hm4Di'};
groups.labels = {'vlOFC','BLA'};
groups.colors = {'#FF2E92','#FF2E92'};

PhasesToShow = [2,4];   % R1 & R3

allpars = cell2mat(winModel.Fit.fitpar');
%     allpars(:,2) = 1./allpars(:,2);   % temperature sigma instead of beta
%     thisModel.ub(2) = 1;
    
figure(32); clf
set(gcf,'Units','normalized','Position',[0,0,1,0.55], 'Color','w');
numfig = 0;
% loop through each lesion group
for g = 1:numel(groups.set)
    % R1/R3
    for i = 1:numel(PhasesToShow)
        phaseNum = PhasesToShow(i);
        phase_label = allPhases{phaseNum}(1);
        schedule_label = allPhases{phaseNum}(2);
        phase_idx = strcmp(winModel.SessInfo.phase,phase_label)&strcmp(winModel.SessInfo.schedules,schedule_label);
        disp(">>> "+i+". "+allPhasesLegends(phaseNum)+" = "+sum(phase_idx));
        % Female/Male
        for k = 1:numel(sexSet)
            group_idx = phase_idx&strcmp(winModel.SessInfo.condition, groups.set{g});
            disp("  ("+g+") "+groups.set{g}+": "+sum(group_idx)+" sessions");
            
            % one panel per each parameter
            for param = 1:numel(winModel.initpar)
                numfig = numfig + 1;
                if mod(numfig,12)==1
                    xPos = 0.0; 
                    x_offset = 0.0;
                else
                    if mod(numfig,3)==1
                        x_offset = 0.015;
                        if mod(numfig,6)==1
                            x_offset = 0.03;
                        end
                    else
                        x_offset = 0.0;
                    end
                end
                SP = subplot(numel(groups.set),numel(PhasesToShow)*2*numel(winModel.initpar),numfig);              
                xPos = xPos + 0.065 + x_offset;
                SP.Position(1) = xPos;
                SP.Position(3) = 0.06;
                
                % Compile data & plot bar graph for each group: CNO/VEH-first
                thisParams = struct;
                MeanVals = nan(1,4);
                SEMvals = nan(1,4);  
                % hm4Di group
                for j = 1:numel(drugOrderSet)
                    sess_idx = group_idx&strcmp(winModel.SessInfo.drugOrder, drugOrderSet(j))&strcmp(winModel.SessInfo.sex, sexSet(k));
                    thisGroup_label = groups.set{g}+"_"+drugOrderSet(j)+"_"+sexSet(k);
                    if param==1; disp("     "+drugOrderSet(j)+"("+k+") "+sexSet(k)+": "+sum(sess_idx)+" params"); end
                    thisParams.(thisGroup_label) = allpars(sess_idx,param);
                    MeanVals(j) = mean(thisParams.(thisGroup_label));
                    SEMvals(j) = std(thisParams.(thisGroup_label))/sqrt(numel(thisParams.(thisGroup_label)));
                    hold on;
                    bar(j,MeanVals(j),'FaceColor',drugOrderColors{i}{j},'LineWidth',1);
                    errorbar(j,MeanVals(j),SEMvals(j),'k','LineStyle','none','HandleVisibility','off');
                    scatter(j*ones(numel(thisParams.(thisGroup_label)),1),thisParams.(thisGroup_label),15,'MarkerFaceColor',drugOrderColors{i}{j},'MarkerFaceAlpha',0.4,'MarkerEdgeColor','k','LineWidth',0.5,'HandleVisibility','off');
                end
                % Place eGFP group next to hm4Di
                for j = 1:numel(drugOrderSet)
                    egfp_idx = phase_idx&contains(winModel.SessInfo.condition, 'eGFP'); % eGFP
                    sess2_idx = egfp_idx&strcmp(winModel.SessInfo.drugOrder, drugOrderSet(j))&strcmp(winModel.SessInfo.sex, sexSet(k));
                    thisGroup_label = "eGFP_"+drugOrderSet(j)+"_"+sexSet(k);
                    if param==1; disp(" vs. eGFP: "+drugOrderSet(j)+"("+k+") "+sexSet(k)+": "+sum(sess2_idx)+" params"); end
                    thisParams.(thisGroup_label) = allpars(sess2_idx,param);
                    MeanVals(j+2) = mean(thisParams.(thisGroup_label));
                    SEMvals(j+2) = std(thisParams.(thisGroup_label))/sqrt(numel(thisParams.(thisGroup_label)));
                    hold on;
                    bar(j+2,MeanVals(j+2),'FaceColor',drugOrderColors{i}{j},'LineWidth',1);
                    errorbar(j+2,MeanVals(j+2),SEMvals(j+2),'k','LineStyle','none','HandleVisibility','off');
                    scatter((j+2)*ones(numel(thisParams.(thisGroup_label)),1),thisParams.(thisGroup_label),15,'MarkerFaceColor',drugOrderColors{i}{j},'MarkerFaceAlpha',0.4,'MarkerEdgeColor','k','LineWidth',0.5,'HandleVisibility','off');
                end
                xticks([1.5, 3.5]);
                %xlim([0.5, 4.5]);
                xticklabels({['\color[rgb]{1 .18 .5725}',groups.labels{g}],'\color[rgb]{.196 .8 .196}eGFP'});
                % xtickangle(30);
                
                % sig-test between two groups
                all_groups = fieldnames(thisParams);
                for j = 1:2
                    [pval, ~, asterisk] = Wilcoxon_ranksum_test(thisParams.(all_groups{2*j-1}), thisParams.(all_groups{2*j}), 1);
                    if pval<.075
                        yPosFact = 0.5;
%                         yPosFact = ceil(max(MeanVals+SEMvals)*10)/10*1.5;
                        plot([2*j-1 2*j], repmat(yPosFact*winModel.ub(param),1,2),'-k', 'LineWidth',1.5);
                        plot([2*j-1 2*j-1], [yPosFact-.05, yPosFact].*winModel.ub(param),'-k', 'LineWidth',1.5);
                        plot([2*j 2*j], [yPosFact-.05 yPosFact].*winModel.ub(param),'-k', 'LineWidth',1.5);
                        if strcmp(asterisk,"n.s.")
                            asterisk = ["\it{p} = \rm"+num2str(pval,3)]; 
                        else
                            asterisk = [asterisk+"\it{p} = \rm"+num2str(pval,3)]; 
                        end
                        text(2*j-0.5, yPosFact*winModel.ub(param), [asterisk], 'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',10);
                    end
                    
                    % hm4Di vs. eGFP
                    [pval, ~, asterisk] = Wilcoxon_ranksum_test(thisParams.(all_groups{j}), thisParams.(all_groups{j+2}), 1);
                    if pval<.075
                        yPosFact = 0.5 + j*0.1;
%                         yPosFact = min(ceil(max(MeanVals+SEMvals)*10)/10*1.5, winModel.ub(param)*.9) + j*0.1*winModel.ub(param);
                        plot([j j+2], repmat(yPosFact*winModel.ub(param),1,2),'-k', 'LineWidth',1.5);
                        plot([j j], [yPosFact-.05, yPosFact].*winModel.ub(param),'-k', 'LineWidth',1.5);
                        plot([j+2 j+2], [yPosFact-.05 yPosFact].*winModel.ub(param),'-k', 'LineWidth',1.5);
                        if strcmp(asterisk,"n.s.")
                            asterisk = ["\it{p} = \rm"+num2str(pval,3)]; 
                        else
                            asterisk = [asterisk+"\it{p} = \rm"+num2str(pval,3)]; 
                        end
                        text(j+1, yPosFact*winModel.ub(param), [asterisk], 'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',10);
                    end
                end            
                              
                text(0.075, 0.95, winModel.plabels(param),'Units','normalized','VerticalAlignment','middle','FontSize',14,'FontWeight','bold');
                if param==1&&k==1
                    legend("CNO1","VEH1",'box','off','fontsize',7);
%                     text(0.15, 1.0, groups.labels{g},'Units','normalized','VerticalAlignment','top');
%                     ylabel(groups.labels{g},'Color',groups.colors{g});
                end
                
                yLimFact = 1;
                ylim([winModel.lb(param), winModel.ub(param)*yLimFact]);
                set(gca,'FontName','Helvetica','FontSize',11,'FontWeight','normal','LineWidth',1, 'tickdir','out','Box','off');
                if g==1&&param==2
                    title(sexSet(k)+" "+allPhasesLegends(phaseNum),'FontSize',14); 
                end
            end
        end
    end
end
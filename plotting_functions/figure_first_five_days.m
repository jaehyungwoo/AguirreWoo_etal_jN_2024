function [outVals, subjectLabels] = figure_first_five_days(block_output, taskType, MET, group_colors, font_size)
outVals = struct;
subjectLabels = struct;

allPhases = {["D","prob1000"],["R1","prob1000"],["R2","prob1000"],["R1","prob9010"],["R2","prob9010"],...
    ["R1","prob8020"],["R2","prob8020"],["R1","prob7030"],["R2","prob7030"]};

allPhasesLabels = ["D","R1","R2","R3","R4","R5","R6","R7","R8"];
allPhasesLegends = ["D","R1 (100/0)","R2 (100/0)","R3 (90/10)","R4 (90/10)",...
    "R5 (80/20)","R6 (80/20)","R7 (70/30)","R8 (70/30)"];

Phases_to_plot = allPhases(1:5);
allPhasesLabels = allPhasesLabels(1:5);

% Groups, including eGFP
lesion_groups = {'all_eGFP','OFC_hm4Di','BLA_hm4Di','ACC_hm4Di'};
lesion_labels = {'eGFP','vOFC','BLA','ACC'};
legend_labels = {'eGFP','OFC hm4Di','BLA hm4Di','ACC hm4Di'};

% first N days to plot
numDays = 5;

%% 1. First 5 days from reversal
    figure;
    set(gcf,'Units','normalized','Position',[0,0,.85,.9],'color','w');
    
    tic
    for g = 1:numel(lesion_groups)
        group_label = lesion_groups{g};
        
        subLabels = fieldnames(block_output.(taskType).(group_label));   
        subN = numel(subLabels);                % track subject N
        
        for i = 1:numel(Phases_to_plot)
            phaseType = Phases_to_plot{i};
            
            %% Subplot
            subplot(numel(lesion_groups),numel(Phases_to_plot),i+(g-1)*numel(Phases_to_plot));            
            
            DataToPlot = struct;
            if i==1
                DataToPlot.D = [];      % for Discrim. phase only
                DataToPlot.D_subjects = [];      % for Discrim. phase only
            else
                DataToPlot.CNO1 = [];   % CNO-first group
                DataToPlot.VEH1 = [];   % VEH-first group    
                DataToPlot.CNO1_subjects = [];
                DataToPlot.VEH1_subjects = [];
            end
            
            % loop through each subject
            for sub = 1:subN
                subjectLabel = subLabels{sub};
                thisDat = block_output.(taskType).(group_label).(subjectLabel);
                
                yy = nan(1,numDays);    % bin for one subject
                
                if ~isfield(thisDat,'pbetter')
                    continue;
                else
                    % find indices for each phase
                    if ~strcmp(phaseType(1),"D")
                        phase_idx = thisDat.(phaseType(1))&thisDat.(phaseType(2));
                    else
                        phase_idx = logical(thisDat.(phaseType(1)));
                    end
                    
                    % assign group: CNO first, or VEH first?
                    if strcmp(phaseType(1),"D")
                        group_name = "D";
                        groups_set = ["D"];
                    else
                        group_name = thisDat.group;
                        groups_set = ["CNO1","VEH1"];
                    end

                    % find respective session day indices
                    days_idx = {};
                    for d = 1:numDays
                        days_idx{d} = phase_idx&(thisDat.SessionDay==d);     %disp(find(days_idx{d}));
                        if sum(days_idx{d})==0
                            yy(d) = NaN;
                        else
                            yy(d) = thisDat.(MET.name)(days_idx{d});
                        end
                    end                    
                end
                
                DataToPlot.(group_name) = [DataToPlot.(group_name); yy];
                DataToPlot.(group_name+"_subjects") = [DataToPlot.(group_name+"_subjects"); convertCharsToStrings(subjectLabel)];
            end
            
            % plot data for each (two) group
            NN = [];
            for S = 1:numel(groups_set)
                thisGroupDat = DataToPlot.(groups_set(S));
                NN(S) = sum(~isnan(sum(thisGroupDat,2,'omitnan')));    % track subject N for legend

                % x-axis: session day 1 through 5
                xx = [1:numDays];
                YY = mean(thisGroupDat,1,'omitnan');   % METHOD 1: Mean of subject means
                SEM = std(thisGroupDat,1,'omitnan')./sqrt(sum(~isnan(thisGroupDat),1));
                              
                outVals.(allPhasesLabels(i)).(lesion_labels{g}).(MET.name).(groups_set(S)) = thisGroupDat;
                subjectLabels.(allPhasesLabels(i)).(lesion_labels{g}).(MET.name).(groups_set(S)) = DataToPlot.(groups_set(S)+"_subjects");
                
                % if no data, fill with NaN
                if isempty(thisGroupDat); YY = nan(size(xx)); SEM = YY; end          
                if length(YY)~=numDays;   error("check y length: "+length(YY)); end
                if length(SEM)~=numDays;   error("check y length: "+length(SEM)); end

                % Plot: Shaded Error bar
                shadedErrorBar(xx,YY,SEM,'lineProps',{'LineWidth',1,'Color',group_colors{S}}); hold on;
            end

            % X-axis setting
            ax = gca;
            ax.Box = 'off';
            ax.XLim = [1 numDays];
            ax.XTick = xx;
            if g==numel(lesion_groups); ax.XLabel.String = "Days"; end
            
            % Y-axis settings
            ax.YLabel.String = MET.label;
            if contains(MET.label,"ERDS")
                ax.YLim = [0.4 1];
                ax.YTick = [0.4:0.2:1];
                ax.YTickLabel = {'.4','.6','.8','1'};
            else
                ax.YLim = [0 1]; 
                ax.YTick = [0:.25:1];
                ax.YTickLabel = {'0','.25','.5','.75','1'};
                plot(ax.XLim, [0.5, 0.5],':k','HandleVisibility','off'); % chance level line
            end  
            
            % Legends
            text(1,1,allPhasesLegends(i),'units','normalized','FontWeight','bold','FontSize',font_size+4,'HorizontalAlignment','right','VerticalAlignment','bottom');
            if i==1
                L = legend([legend_labels{g},' (n=',num2str(sum(NN)),')']);
                L.FontSize = font_size - 3;
                L.Box = 'off';
            else
                L = legend({['CNO1 (n=',num2str(NN(1)),')'],['VEH1 (n=',num2str(NN(2)),')']});
                L.FontSize = font_size - 3;
                L.Box = 'off';
                L.Location = 'southwest';
            end
            
            % Figure panel labels
            if i==1
                text(-0.7,1.,['\bf',char(g+'A'-1),'   \rm',lesion_labels{g}],'FontName','Helvetica','FontSize',font_size+8,'FontWeight','bold',...
                    'Units','normalized','HorizontalAlignment','left','VerticalAlignment','bottom');
            end
            
            % set fonts, axis line width, tick directions
            set(ax,'FontName','Helvetica','FontSize',font_size,'FontWeight','normal','LineWidth',1,'tickdir', 'out','Box','off')
        end
    end
    toc
    
end
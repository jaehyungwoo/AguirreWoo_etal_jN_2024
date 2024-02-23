function [outVals, subjectLabels] = sigTest_days_around_rev(block_output, taskType, allMet, group_colors, font_size)
outVals = struct;
subjectLabels = struct;
MET = struct;

allPhases = {["D","prob1000"],["R1","prob1000"],["R2","prob1000"],["R1","prob9010"],["R2","prob9010"],...
    ["R1","prob8020"],["R2","prob8020"],["R1","prob7030"],["R2","prob7030"]};
allPhasesLabels = ["D","R1","R2","R3","R4","R5","R6","R7","R8"];
allPhasesLegends = ["D","R1(100/0)","R2(100/0)","R3(90/10)","R4(90/10)",...
    "R5(80/20)","R6(80/20)","R7(70/30)","R8(70/30)"];

Revs_to_test = [1:5];

% Including eGFP
lesion_groups = {'all_eGFP','OFC_hm4Di','BLA_hm4Di','ACC_hm4Di'};
lesion_labels = {'eGFP','vOFC','BLA','ACC'};
legend_labels = {'eGFP','OFC hm4Di','BLA hm4Di','ACC hm4Di'};

numDays = 3;    % last N-1 days before rev, first N days after rev

%% 2. Days around reversal
%     figure;
%     set(gcf,'Units','normalized','Position',[0,0,.85,0.9],'color','w');

for m = 1:numel(allMet)   
    MET.name = allMet(m);
    tic
    for g = 1:numel(lesion_groups)
        group_label = lesion_groups{g}; 
        disp("----------------------------"+group_label+"----------------------------");
        
        subLabels = fieldnames(block_output.(taskType).(group_label));   
        subN = numel(subLabels);                % track subject N
        
        for i = 1:numel(Revs_to_test)
            beforeR = allPhases{i};       % before rev
            afterR = allPhases{i+1};      % after rev
            
            %% Subplot
%             subplot(numel(lesion_groups),numel(Revs_to_test),i+(g-1)*numel(Revs_to_test));

            DataToPlot = struct;
            DataToPlot.CNO1 = [];   % CNO-first group
            DataToPlot.VEH1 = [];   % VEH-first group    
            DataToPlot.CNO1_subjects = [];
            DataToPlot.VEH1_subjects = [];
            
            % loop through each subject
            for sub = 1:subN
                subjectLabel = subLabels{sub};
                thisDat = block_output.(taskType).(group_label).(subjectLabel);
                
                yy = nan(1,numDays*2+1);    % bin for one subject
                if ~isfield(thisDat,'pbetter')
                    continue;
                else
                    % assign group: CNO first, or VEH first?
                    group_name = thisDat.group;
                    groups_set = ["CNO1","VEH1"];

                    % find indices for each phase
                    phase_idx1 = thisDat.(beforeR(1))&thisDat.(beforeR(2));
                    phase_idx2 = thisDat.(afterR(1))&thisDat.(afterR(2));
                    
                    % find respective session day indices:                    
                    % before Rev data
                    last_idx = find(phase_idx1);
                    if length(last_idx)>=numDays
                        yy(1:1+numDays-1) = thisDat.(MET.name)(last_idx(end-numDays+1:end));
                    else
                        numEmpt = numDays - length(last_idx);
                        yy(1:1+numDays-1) = [nan(numEmpt,1); thisDat.(MET.name)(last_idx(end-length(last_idx)+1:end))];
                    end
                    % after Rev data
                    first_idx = find(phase_idx2);
                    if length(first_idx)>=numDays+1
                        yy(1+numDays:numDays*2+1) = thisDat.(MET.name)(first_idx(1:numDays+1));
                    else
                        numEmpt = numDays+1 - length(first_idx);
                        yy(1+numDays:numDays*2+1) = [thisDat.(MET.name)(first_idx(1:length(first_idx))); nan(numEmpt,1)];
                    end
                end
                
                DataToPlot.(group_name) = [DataToPlot.(group_name); yy];
                DataToPlot.(group_name+"_subjects") = [DataToPlot.(group_name+"_subjects"); convertCharsToStrings(subjectLabel)];
            end
            
            %% Compare data between two group
            NN = [];
            
            Dat1_before = DataToPlot.CNO1(:,1:numDays);
            Dat2_before = DataToPlot.VEH1(:,1:numDays);
            [pval, cohensD, asterisk, stats] = two_sample_T_test(Dat1_before(:), Dat2_after(:), 1, 'both');

            
            for S = 1:numel(groups_set)
                thisGroupDat = DataToPlot.(groups_set(S));
                NN(S) = sum(~isnan(sum(thisGroupDat,2,'omitnan')));
                
                % set x-axis: N+(N+1)
                xx = [1:numDays*2+1];
                YY = mean(thisGroupDat,1,'omitnan');   % METHOD 1: Mean of subject means
                SEM = std(thisGroupDat,1,'omitnan')./sqrt(sum(~isnan(thisGroupDat),1));
                
                outVals.(allPhasesLabels(i)).(lesion_labels{g}).(MET.name).(groups_set(S)) = thisGroupDat;
                subjectLabels.(allPhasesLabels(i)).(lesion_labels{g}).(MET.name).(groups_set(S)) = DataToPlot.(groups_set(S)+"_subjects");
                
                % if no data, fill with NaN
                if isempty(thisGroupDat); YY = nan(size(xx)); SEM = YY; end          
                if length(YY)~=numDays*2+1;   error("check y length: "+length(YY)); end
                if length(SEM)~=numDays*2+1;   error("check y length: "+length(SEM)); end

                % Plot: Shaded Error bar
                shadedErrorBar(xx,YY,SEM,'lineProps',{'LineWidth',1,'Color',group_colors{S}}); hold on;
            end

            % X-axis setting
            ax = gca;
            ax.XLim = [1 numDays*2+1];
            ax.XTick = xx;
            ax.XTickLabel = [-numDays:1:numDays];
            if g==numel(lesion_groups); ax.XLabel.String = "Days"; end
            % Y-axis settings
            ax.YLabel.String = MET.label;
            if contains(MET.label,"ERDS")||contains(MET.label,"H_")
                ax.YLim = [0.4 1];
                ax.YTick = [0.4:0.2:1];
                ax.YTickLabel = {'.4','.6','.8','1'};
            else
                ax.YLim = [0 1]; 
                ax.YTick = [0:.25:1];
                ax.YTickLabel = {'0','.25','.5','.75','1'};
                plot(ax.XLim, [0.5, 0.5],':k','HandleVisibility','off'); % chance level line
            end  
            % reversal line & labels
            plot([numDays,numDays]+1,ax.YLim,'--k','LineWidth',1,'HandleVisibility','off'); 
            text(0,1,allPhasesLegends(i),'units','normalized','FontWeight','bold','FontSize',font_size+3,'HorizontalAlignment','left','VerticalAlignment','bottom');
            text(1,1,allPhasesLegends(i+1),'units','normalized','FontWeight','bold','FontSize',font_size+3,'HorizontalAlignment','right','VerticalAlignment','bottom');            
            % Legends
            L = legend({['CNO1 (n=',num2str(NN(1)),')'],['VEH1 (n=',num2str(NN(2)),')']});  % legend({['R1-CNO'],['R2-CNO']});
            L.Location = 'southwest';
            L.FontSize = font_size - 3;
            L.Box = 'off';      
            % Figure panel labels
            if i==1
                text(-0.7,1.,['\bf',char(g+'A'-1),'   \rm',lesion_labels{g}],'FontName','Helvetica','FontSize',font_size+8,'FontWeight','bold',...
                    'Units','normalized','HorizontalAlignment','left','VerticalAlignment','bottom');
            end
            
            % set fonts, axis line width, tick directions
            set(ax,'FontName','Helvetica','FontSize',font_size,'FontWeight','normal','LineWidth',1,'tickdir', 'out','Box','off'); %,'FontWeight','bold'
        end 
    end
    toc
    
end
end
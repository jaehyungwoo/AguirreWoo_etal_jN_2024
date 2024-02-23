function [outVal, subjectLabels] = figure_trials_around_rev(trial_output, plot_type, taskType, MET, group_colors, font_size, smooth_window)
outVal = struct;
subjectLabels = struct;

if ~exist('smooth_window','var')
    smooth_window = 1; 
end

allPhases = {["D","prob1000"],["R1","prob1000"],["R2","prob1000"],["R1","prob9010"],["R2","prob9010"],...
    ["R1","prob8020"],["R2","prob8020"],["R1","prob7030"],["R2","prob7030"]};
allPhasesLabels = ["D","R1","R2","R3","R4","R5","R6","R7","R8"];
allPhasesLegends = ["D","R1(100/0)","R2(100/0)","R3(90/10)","R4(90/10)",...
    "R5(80/20)","R6(80/20)","R7(70/30)","R8(70/30)"];

% Revs_to_plot = [1:5];
Revs_to_plot = [1:5];

% Including eGFP
% lesion_groups = {'all_eGFP','OFC_hm4Di','BLA_hm4Di','ACC_hm4Di'};
% lesion_labels = {'eGFP','vOFC','BLA','ACC'};
% legend_labels = {'eGFP','OFC hm4Di','BLA hm4Di','ACC hm4Di'};

lesion_groups = {'all_eGFP','OFC_hm4Di','BLA_hm4Di'};
lesion_labels = {'eGFP','vOFC','BLA'};
legend_labels = {'eGFP','OFC hm4Di','BLA hm4Di'};

% N trials to plot around reversal
numTrials = trial_output.N;

    %% 4. Trials around reversal
    figure; 
    set(gcf,'Units','normalized','Position',[0,0,.85,.90],'color','w');
    
    tic
    for g = 1:numel(lesion_groups)
        group_label = lesion_groups{g};     %disp(g+". "+group_label);
        
        subLabels = fieldnames(trial_output.(taskType).(group_label));   
        subN = numel(subLabels);
        
        for i = 1:numel(Revs_to_plot)
            beforeR = allPhases{i};       % before rev
            afterR = allPhases{i+1};      % after rev
            
            %% Subplot
            subplot(numel(lesion_groups),numel(Revs_to_plot),i+(g-1)*numel(Revs_to_plot));
            
            DataToPlot = struct;
            DataToPlot.CNO1 = [];   % CNO-first group
            DataToPlot.VEH1 = [];   % VEH-first group    
            DataToPlot.CNO1_subjects = [];
            DataToPlot.VEH1_subjects = [];
            
            % loop through each subject
            for sub = 1:subN
                subjectLabel = subLabels{sub};
                thisDat = trial_output.(taskType).(group_label).(subjectLabel);
                
                yy = nan(1,numTrials*2);    % bin for one subject
                if ~isfield(thisDat,'firstN')||~isfield(thisDat.firstN.(plot_type),'pbetter')
                    continue;
                else
                    % assign group: CNO first, or VEH first?
                    group_name = thisDat.group;
                    groups_set = ["CNO1","VEH1"];
                    
                    % find indices for each phase
                    phase_idx1 = thisDat.beforeRev.(beforeR(1))&thisDat.beforeRev.(beforeR(2));
                    phase_idx2 = thisDat.firstN.(afterR(1))&thisDat.firstN.(afterR(2));
                    
                    % find respective session day indices:                    
                    % before Rev data
                    if sum(phase_idx1)~=0
                        yy(1:numTrials) = thisDat.beforeRev.(plot_type).(MET.name)(phase_idx1,:);
                    end
                    % after Rev data
                    if sum(phase_idx2)~=0
                        yy(numTrials+1:numTrials*2) = thisDat.firstN.(plot_type).(MET.name)(phase_idx2,:);
                    end
                end
                
                DataToPlot.(group_name) = [DataToPlot.(group_name); yy];
                DataToPlot.(group_name+"_subjects") = [DataToPlot.(group_name+"_subjects"); convertCharsToStrings(subjectLabel)];
            end
            
            % plot data for each (two) group
            NN = [];
            for S = 1:numel(groups_set)
                thisGroupDat = DataToPlot.(groups_set(S));
                NN(S) = sum(~isnan(sum(thisGroupDat,2,'omitnan')));
                
                % set x-axis: N+(N+1)
                xx = [1:numTrials*2];
                YY = mean(thisGroupDat,1,'omitnan');   % METHOD 1: Mean of subject means
                SEM = std(thisGroupDat,1,'omitnan')./sqrt(sum(~isnan(thisGroupDat),1));
                
                % Additional smoothing: default is none (i.e., window of 1)
                % Caution: Smooth Before/After rev separately (o/w non-existent NaN values gets smoothed)
                y1 = smooth(YY(1:numTrials), smooth_window);
                y2 = smooth(YY(numTrials+1:end), smooth_window);
                YY = [y1; y2]';
                sem1 = smooth(SEM(1:numTrials), smooth_window);
                sem2 = smooth(SEM(numTrials+1:end), smooth_window);
                SEM = [sem1; sem2]';
                
                % output raw values
                outVal.(allPhasesLabels(i)).(lesion_labels{g}).(MET.name).(groups_set(S)) = thisGroupDat;
                subjectLabels.(allPhasesLabels(i)).(lesion_labels{g}).(MET.name).(groups_set(S)) = DataToPlot.(groups_set(S)+"_subjects");
                
                % if no data, fill with NaN
                if isempty(thisGroupDat); YY = nan(size(xx)); SEM = YY; end          
                if length(YY)~=numTrials*2;   error("check y length: "+length(YY)); end
                if length(SEM)~=numTrials*2;   error("check y length: "+length(SEM)); end
                
                % Plot: Shaded Error bar
                shadedErrorBar(xx,YY,SEM,'lineProps',{'LineWidth',1,'Color',group_colors{S}}); hold on;
            end

            % X-axis setting
            ax = gca;
            ax.XLim = [1 numTrials*2+1];
            ax.XTick = [1:50:numTrials*2+1];
            ax.XTickLabel = [-numTrials:50:numTrials];
            if g==numel(lesion_groups); ax.XLabel.String = "Trials"; end
            
            % Y-axis settings
            ax.YLabel.String = MET.label;
            if contains(MET.name,"RI_")
                ax.YLim = [-.5 .5];
                ax.YTick = [ax.YLim(1):ax.YLim(2)/2:ax.YLim(2)];
                plot(ax.XLim, [0, 0],':k','HandleVisibility','off'); % baseline
            else
                ax.YLim = [0 1]; 
                ax.YTick = [0:.25:1];
                ax.YTickLabel = {'0','.25','.5','.75','1'};
                if ~(contains(MET.label,"ERDS")||contains(MET.label,"H_"))
                    plot(ax.XLim, [0.5, 0.5],':k','HandleVisibility','off'); % chance level line
                end
            end
            % reversal line & labels
            plot([numTrials,numTrials]+1,ax.YLim,'--k','LineWidth',1,'HandleVisibility','off'); 
            text(0,1,allPhasesLegends(i),'units','normalized','FontWeight','bold','FontSize',font_size+3,'HorizontalAlignment','left','VerticalAlignment','bottom');
            text(1,1,allPhasesLegends(i+1),'units','normalized','FontWeight','bold','FontSize',font_size+3,'HorizontalAlignment','right','VerticalAlignment','bottom');            

            
            if 1%i==2
            % Legends    
                L = legend({['CNO1 (n=',num2str(NN(1)),')'],['VEH1 (n=',num2str(NN(2)),')']});  % legend({['R1-CNO'],['R2-CNO']});
                L.Location = 'southwest';
                L.FontSize = font_size - 3;
                L.Box = 'off';
            
            % Figure panel labels
%                 text(-0.7,1.,['\bf',char(g+'A'-1),'   \rm',lesion_labels{g}],'FontName','Helvetica','FontSize',font_size+8,'FontWeight','bold',...
%                     'Units','normalized','HorizontalAlignment','left','VerticalAlignment','bottom');
            end
            
            % set fonts, axis line width, tick directions
            set(ax,'FontName','Helvetica','FontSize',font_size,'FontWeight','normal','LineWidth',1,'tickdir', 'out','Box','off'); %,'FontWeight','bold'
        end 
    end  
%     sgtitle(taskType+" task");
    toc
end
 % analyze rat behavior
close all; clc; clearvars

% Set task info
task_subsets = ["SO","AO"];
schedule_subsets = ["prob1000", "prob9010", "prob8020", "prob7030"];
stage_subsets = ["D","R1","R2"];
injection_subsets = ["none","CNO","VEH"];

% load preprocessed data 
load(["dataset/preprocessed/allAnimalsDat.mat"],'allDat','conditions_set');   % all groups
lesion_groups = unique(conditions_set);
animalLabels = fieldnames(allDat);

% Set file names
outputFname1 = "output/plot_data/block_output";
outputFname2 = "output/plot_data/trial_output";

% compute settings:
calculate_whole_block_stats = 1;    % whole day stats
calculate_trials_beforeafter = 1;   % calculate across trials stats

if calculate_trials_beforeafter
    beforeAfterN = 100;            % N trials before/after reversal
end

%% Calculate output    
tic

if 1
    % intialize data bin
    block_output = struct;
    trial_output = struct;
    trial_output.N = beforeAfterN;
    
    numGroups = length(lesion_groups);
    for i = 1:numGroups
        group_label = lesion_groups(i);
        if contains(group_label,"eGFP"); group_label = "all_eGFP"; end
        for task = task_subsets
            block_output.(task).(group_label) = struct;
            trial_output.(task).(group_label) = struct;
        end
    end
    
    %% loop through each group/animal and compute metrics
    for tt = 1:numel(task_subsets)
        disp("Task "+tt+": "+task_subsets(tt));  fprintf('\n');
        all_blockcnt = 0;
        
        % loop through each lesion group
        for g = 1:numGroups
            %%
            lesionGroup_label = lesion_groups(g);
            group_idx = find(strcmp(conditions_set,lesionGroup_label))';  % Subjects ID for this lesion group
            disp("========="+lesionGroup_label+": "+length(group_idx)+" animals=========");
            
            group_label = lesionGroup_label;
            if contains(group_label,"eGFP"); group_label = "all_eGFP"; end  % save all eGFP animals into one group
            
            % loop through each subject
            groupSub_cnt = 0;
            for n = [group_idx]
                subjectLabel = animalLabels{n};
                groupSub_cnt = groupSub_cnt + 1;
                
                tempSubBlock = struct;
                tempSubTrials = struct;
                tempSubTrials.beforeRev.Cum = struct;
                tempSubTrials.firstN.Cum = struct;
                tempSubTrials.beforeRev.Run = struct;
                tempSubTrials.firstN.Run = struct;
                
                % group assignment
                if strcmp(allDat.(subjectLabel).R1,'CNO')
                    tempSubBlock.group = "CNO1";
                    tempSubTrials.group = "CNO1";
                elseif strcmp(allDat.(subjectLabel).R1,'VEH')
                    tempSubBlock.group = "VEH1";
                    tempSubTrials.group = "VEH1";
                end
                
                %% One subject data
                fprintf("%d. Subject %d (%s): ",[groupSub_cnt,n,subjectLabel]);
                
                this_animal_dat = allDat.(subjectLabel).(task_subsets(tt));
                lesion_cond = allDat.(animalLabels{n}).condition;
                if ~strcmp(lesionGroup_label, lesion_cond); error('lesion label error'); end
                
                % loop through each phase
                stages = fieldnames(this_animal_dat);
                for i = 1:numel(stages)
                    this_stage_dat = this_animal_dat.(stages{i});

                    %% loop through each session day
                    for j = 1:numel(this_stage_dat)
                        all_blockcnt = all_blockcnt + 1;
                        
                        block_stats = this_stage_dat{j};
                        if isempty(block_stats); continue; end  
                        
                        % whole-day metrics
                        if calculate_whole_block_stats
                            output = ComputeWholeBlockStats(block_stats, allDat.(subjectLabel), task_subsets(tt));
                            output.subjectLabel = string(subjectLabel);
                            output.SessionDay = j;
                            tempSubBlock = append_to_fields(tempSubBlock, {output});
                        end
                        
                        % Across-trial plot: first & last N trials data                        
                        if calculate_trials_beforeafter
                            if j==numel(this_stage_dat)
                                trialsType = "beforeRev";      % Last N trials from the last day before reversal
                                prevBlock_temp = block_stats; % keep record of block stats before reversal 
                                prevBlock_stats = [];
                                compute_trials_flag = 1;
                            elseif j==1%&&i~=1
                                trialsType = "firstN";      % Frist N tirals (after reversal, except for D phase)
                                if i==1
                                    prevBlock_stats = [];
                                else
                                    prevBlock_stats = prevBlock_temp;
                                end 
                                compute_trials_flag = 1;
                            else
                                compute_trials_flag = 0;
                            end
                            
                            if compute_trials_flag
                                [outputCum, outputRun, output_idx] = ComputeAcrossTrialsStats(block_stats, allDat.(subjectLabel), beforeAfterN, trialsType, prevBlock_stats);
                                output_idx.subjectLabel = string(subjectLabel);
                                output_idx.SessionDay = j;
                                tempSubTrials.(trialsType).Cum = append_to_fields(tempSubTrials.(trialsType).Cum, {outputCum});
                                tempSubTrials.(trialsType).Run = append_to_fields(tempSubTrials.(trialsType).Run, {outputRun});
                                tempSubTrials.(trialsType) = append_to_fields(tempSubTrials.(trialsType), {output_idx});
                            end
                        end

                        display_counter(j);
                    end
                end
                fprintf('\n');
                block_output.(task_subsets(tt)).(group_label).(subjectLabel) = tempSubBlock;
                trial_output.(task_subsets(tt)).(group_label).(subjectLabel) = tempSubTrials;
            end
            fprintf('\n');
        end
        disp("Total block count: "+all_blockcnt);
        fprintf('\n\n');
    end
    
    % save output
    if calculate_whole_block_stats
        save(outputFname1+".mat",'block_output'); 
    end
    if calculate_trials_beforeafter 
        save(outputFname2+beforeAfterN+".mat",'trial_output'); 
    end
    disp("File saved!");
end 

ET = toc;
disp("Elaspsed time is "+ET/60+" minutes.");


function [block_stats] = preprocessBlockData(rawData, sched_set, taskType, phase_label, sched_label, include_noChoiceTrials)
    if ~exist('include_noChoiceTrials','var')
        include_noChoiceTrials = false;
    end

    % obtain raw trial vectors from table
    committed_choice_idx =  logical(table2array(rawData(:,'PerTrial_CommitedChoices')));    % trials where a choice has been made
    if include_noChoiceTrials
        reward_vec = table2array(rawData(:,'PerTrial_Rewarded'));
        reward_vec(~committed_choice_idx) = NaN;
        correct_vec = table2array(rawData(:,'PerTrial_Correct'));
        correct_vec(~committed_choice_idx) = NaN;
        correctPosition = table2array(rawData(:,'PerTrial_PositionOfCorrectImage')); % 1 (=left) or 3 (=right)
        correctPosition = double(correctPosition==3);               % 0 (=left) or 1 (=right)
        numTrials = length(correctPosition);
    else
        reward_vec = table2array(rawData(committed_choice_idx,'PerTrial_Rewarded'));
        correct_vec = table2array(rawData(committed_choice_idx,'PerTrial_Correct'));
        correctPosition = table2array(rawData(committed_choice_idx,'PerTrial_PositionOfCorrectImage')); % 1 (=left) or 3 (=right)
        correctPosition = double(correctPosition==3);               % 0 (=left) or 1 (=right)
        numTrials = length(correctPosition);
    end

    % initialize data bin for the whole block
    block_stats = struct;
    block_stats.r = reward_vec;                             % reward(=1) or not(=0)
    block_stats.cloc = zeros(numTrials,1);                    % choice location: 0=Left, 1=Right
    block_stats.hr_side = correctPosition;                  % higher reward ("correct") side: 0=Left, 1=Right  
    block_stats.hr_side(block_stats.hr_side==0) = -1;       % -1=Left, 1=Right  

    % assign choice location (L/R)
    block_stats.cloc(correctPosition==1&correct_vec==1) = 1;       % if correct side was Right and animal was correct, chosen side is Right(=1)
    block_stats.cloc(correctPosition==1&correct_vec==0) = -1;     % if correct side was Right and animal was wrong, chosen side is Left(=-1)
    block_stats.cloc(correctPosition==0&correct_vec==1) = -1;     % if correct side was Left and animal was correct, chosen side is Left(=1)
    block_stats.cloc(correctPosition==0&correct_vec==0) = 1;     % if correct side was Left and animal was wrong, chosen side is Right(=1)

    if strcmp(taskType,"SO")
        % Stimulus-Outcome task: note AO task does not have stimulus choice                  
        % assign choice stimulus (A/B)
        block_stats.cstim = zeros(numTrials,1);      % choice stimulus: 0=Left, 1=Right
        switch phase_label
            case "D"
                better_stim = -1;       % assume stimA (arbitrary) was a better option at D
            case "R1"
                better_stim = 1;        % this switches to stimB at R1
            case "R2"
                better_stim = -1;       % back to stimA at R2
        end
        block_stats.cstim(correct_vec==1) = better_stim;              
        block_stats.cstim(correct_vec==0) = -better_stim;
        block_stats.hr_stim = better_stim*ones(numTrials,1); 
    end

    % assign reward schedule & task info
    for s = 1:numel(sched_set) 
        block_stats.(sched_set(s))= false;
    end
    block_stats.(sched_label) = true;       % indicate reward schedule for this block
    block_stats.stage = phase_label;
    block_stats.task = taskType;

end
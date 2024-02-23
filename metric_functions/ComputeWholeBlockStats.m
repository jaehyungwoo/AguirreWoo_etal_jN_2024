function [output] = ComputeWholeBlockStats(block_stats, animalInfo, task)
    output = struct;
    
    reward = block_stats.r;
    choice_loc = block_stats.cloc;
    hr_side = block_stats.hr_side;
    
    if strcmp(task,"SO")
        choice_vec = block_stats.cstim;
        better_opt = block_stats.hr_stim;
        choseBetter = choice_vec==better_opt;
        str_loc = choice_loc(1:end-1)==choice_loc(2:end);
        LocBased = Conditional_Entropy(str_loc, reward(1:end-1),"ERDS_loc");
        LocBased.H_str_loc = Shannon_Entropy(str_loc);
        
        chooseLeft = (choice_loc==-1);
        LocBased.pLeft = mean(chooseLeft);
        LocBased.pstay_loc = mean(str_loc);
        chooseLeft = chooseLeft(1:end-1);
        LocBased.RI_LR = LocBased.pstay_loc - (mean(chooseLeft)^2+(1-mean(chooseLeft))^2);
        LocBased.RI_L = mean(str_loc&chooseLeft) - LocBased.pLeft^2;
        LocBased.RI_R = mean(str_loc&~chooseLeft) - (1-LocBased.pLeft)^2;
        
        % compute metrics: HR based on stimulus
        output = append_to_fields(output, {compute_metrics_rats(choice_vec, reward, better_opt), LocBased});            
        
    elseif strcmp(task,"AO")
        choice_vec = choice_loc;
        better_opt = hr_side;
        
        % compute metrics: HR based on location
        output = compute_metrics_rats(choice_vec, reward, better_opt);                        
    else
        error("Set correct task name");        
    end

    %% mark schedule label, phase, injection type, and subject name
    
    schedule_subsets = ["prob1000", "prob9010", "prob8020", "prob7030"];
    stage_subsets = ["D","R1","R2"];
    injection_subsets = ["none","CNO","VEH"];
    
    for s = 1:numel(schedule_subsets)
        output.(schedule_subsets(s)) = block_stats.(schedule_subsets(s));
    end
    for s = 1:3
        output.(stage_subsets(s)) = 0;
        output.(injection_subsets(s)) = 0;
    end
    output.(block_stats.stage) = 1;
    
    switch block_stats.stage
        case "D"
            output.none = 1;
        case "R1"
            output.(animalInfo.R1) = 1;
        case "R2"
            output.(animalInfo.R2) = 1;
    end
    output.sex = string(animalInfo.sex);

end
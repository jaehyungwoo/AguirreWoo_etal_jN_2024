function [idx_chr,sched_labels,stage_labels] = SessionsChronologicalSort(subFolderNames, Animal_Data)
% Sort phases by session order chronologically
% output returns:
%   idx_chr: chrnological phase index
%   sched_labels: schedule labels in a chronological order
%   stage_labels: phase labels in a chronological order
    
    idx_chr = nan(1,numel(subFolderNames));
    sched_labels = strings(1,numel(subFolderNames));
    stage_labels = strings(1,numel(subFolderNames));
   
    idx_chr(1) = find(strcmp(subFolderNames,'D'));    % D: baseline discrimination task, w/o any injections
    sched_labels(1) = "prob1000";
    stage_labels(1) = "D";

    if length(idx_chr)>=2
        idx_chr(2) = find(strcmp(subFolderNames,[Animal_Data.R1,'100']));   % R1 deterministic (100/0)
        idx_chr(3) = find(strcmp(subFolderNames,[Animal_Data.R2,'100']));   % R2 deterministic (100/0)
        sched_labels(2:3) = "prob1000";
        stage_labels(2) = "R1";
        stage_labels(3) = "R2";
        if length(idx_chr)>=4
            idx_chr(4) = find(strcmp(subFolderNames,[Animal_Data.R1,'90']));
            sched_labels(4) = "prob9010";
            stage_labels(4) = "R1";
            if length(idx_chr)>=5   % fix: some subjects did not perfom second reversals
                idx_chr(5) = find(strcmp(subFolderNames,[Animal_Data.R2,'90']));
                sched_labels(5) = "prob9010";
                stage_labels(5) = "R2";
                if length(idx_chr)>=6
                    idx_chr(6) = find(strcmp(subFolderNames,[Animal_Data.R1,'80']));
                    sched_labels(6) = "prob8020";    
                    stage_labels(6) = "R1";
                    if length(idx_chr)>=7
                        idx_chr(7) = find(strcmp(subFolderNames,[Animal_Data.R2,'80']));
                        sched_labels(7) = "prob8020";                        
                        stage_labels(7) = "R2";
                        if length(idx_chr)>=8
                            idx_chr(8) = find(strcmp(subFolderNames,[Animal_Data.R1,'70']));
                            sched_labels(8) = "prob7030";
                            stage_labels(8) = "R1";
                            if length(idx_chr)>=9     % fix: some subjects did not perfom second reversals for the same reward condition
                                idx_chr(9) = find(strcmp(subFolderNames,[Animal_Data.R2,'70']));
                                sched_labels(9) = "prob7030";
                                stage_labels(9) = "R2";
                            end
                        end
                    end
                end
            end
        end
    end

end
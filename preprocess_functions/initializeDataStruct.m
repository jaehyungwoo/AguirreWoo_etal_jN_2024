function [allDat, sched_set, conditions_set] = initializeDataStruct(CohortStatus, task_set, unqualifiedSub)
    % Obtain subject label numbers from existing folder names
    cd('Animal_Data');

    allDat = struct;                        % initialize data structure
    sched_set = ["prob1000","prob9010","prob8020","prob7030"];

    subj = dir;
    AnimalNum = [];     
    for i = 1:length(subj)
       this_name = subj(i).name;
       if contains(this_name,'.'); continue; end
       if contains(this_name,'CA')
           AnimalNum(end+1,1) = str2double(this_name(3:end));
       end
    end
    
    % Subjects to be excluded from analysis
%     % unqualifiedSub = [14, 16, 23, 29, 35:39, 41, 43, 51:53];   % based on excel file
    CA_excluded = [unqualifiedSub];     
    
    AnimalNum = setdiff(AnimalNum, CA_excluded);      
    AnimalNum = sort(AnimalNum,'ascend');   % sort subject numbers
    
    for n = 1:length(AnimalNum); allDat.("CA"+AnimalNum(n)) = struct; end

    % Assign subject info from the excel file
    TotalSubjLabels = CohortStatus.Properties.RowNames;
    AnimalLabel = fieldnames(allDat);
    disp("Excluded subject: "+cell2mat(setdiff(TotalSubjLabels,AnimalLabel)));   % reporting excluded subjects
    conditions_set = strings(numel(AnimalLabel),1);
    for i = 1:numel(AnimalLabel)
        this_idx = find(strcmp(TotalSubjLabels,AnimalLabel{i}));
        if isempty(this_idx); error(i+": "+AnimalLabel{i}+" data missing"); end
        this_cond = table2array(CohortStatus(this_idx,2:3));                % brain region and virus group

        if contains(this_cond{2},'4Di'); this_cond{2} = 'hm4Di'; end        % correcting for upper/lower case
    %     if strcmp(this_cond{2},'hM4Di'); this_cond{2} = 'hm4Di'; end        

        allDat.(AnimalLabel{i}).condition = [this_cond{1},'_',this_cond{2}]; % disp(allDat.(AnimalLabel{i}).condition);
        allDat.(AnimalLabel{i}).sex = char(table2array(CohortStatus(this_idx,1)));

        injections = table2array(CohortStatus(this_idx,4:5));
        allDat.(AnimalLabel{i}).R1 = injections{1};
        allDat.(AnimalLabel{i}).R2 = injections{2};
        
        for task = [task_set]
            allDat.(AnimalLabel{i}).(task) = struct;
        end
%         allDat.(AnimalLabel{i}).SO = struct;
%         allDat.(AnimalLabel{i}).AO = struct;
        
        conditions_set(i) = string([this_cond{1},'_',this_cond{2}]);
    end

    
    

end
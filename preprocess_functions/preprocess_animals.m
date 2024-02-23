%% preprocess animal data
close all; clc; clearvars

% Root to dataset directory and load the excel file first
if contains(pwd,'wojh')
    cd('C:\Users\wojh1\Dropbox (Dartmouth College)\CCNL\Alicia_Data\Claudia_SO_AO_data');
else
    cd('~'); cd('Dropbox (Dartmouth College)/CCNL/Alicia_Data/Claudia_SO_AO_data');
end

cd('dataset');
statusExcelFname = 'Cohort_Status_Sep2023.xlsx';
CohortStatus = readtable(statusExcelFname,'ReadRowNames',true);
disp("Loaded Cohort Status excel file");

task_set = ["SO","AO"];     % 'Stimulus-Outcome' / 'Action-Outcome'
Sub_to_exclude = [14, 16, 23, 29, 35:39, 41, 43, 51:53, 66, 68, 69, 70, 72:81, 87];   % based on excel file

[allDat, sched_set, conditions_set] = initializeDataStruct(CohortStatus, task_set, Sub_to_exclude);     % specify excluded subjects in this function

include_noChoiceTrials = 0;     % include data for uncommitted trials (for model fitting purposes)

%% loop through each animal and load dataset
tic

disp("Preprocessing starting...");
AnimalLabel = fieldnames(allDat);
disp("Total "+numel(AnimalLabel)+" subjects");
for k = 1:numel(AnimalLabel)
    fprintf('> Subject %d (%s)', [k, AnimalLabel{k}]);
    
    % loop through each task
    for t = 1:numel(task_set)
        fprintf('\n '+task_set(t));
        dataset_address = [string(AnimalLabel{k})+"/"+task_set(t)];
        
        % Skip if the corresponding task folder does not exist
        if ~exist(dataset_address,'dir');  continue;   end
        
        sessfolders = dir(dataset_address);
        dirFlags = [sessfolders.isdir];
        subFolders = sessfolders(dirFlags);
        subFolderNames = {subFolders(3:end).name};
        
        [idx_chr,sched_labels,stage_labels] = SessionsChronologicalSort(subFolderNames, allDat.(AnimalLabel{k}));
        
        %% loop through each phase (folder)
        for i = 1:numel(subFolderNames)
            sessfiles = dir([dataset_address+"/"+subFolderNames{idx_chr(i)}]);
            dirFlags = ~[sessfiles.isdir];
            subFiles = sessfiles(dirFlags);
            subFileNames = {subFiles(1:end).name};
            
            if isempty(subFileNames)
                % phase with no data
                allDat.(AnimalLabel{k}).(task_set(t)).(subFolderNames{idx_chr(i)}) = {};
                continue;
            end
            
            % Edit: extract day number and sort accordingly (w/o this step, day #10 comes in second)
            error_flag = 0;
            % Also check for missing datafiles
            dayNum = zeros(numel(subFileNames),1);
            for j = 1:numel(subFileNames)
                str = subFileNames{j};
                numStart = regexp(str,'D\d.');
                numEnd = regexp(str,'.xls');
                dayNum(j) = str2double(str(numStart+1:numEnd-1));
            end
            [dayNum,idx] = sort(dayNum, 'ascend');
            subFileNames = subFileNames(idx)';
            if dayNum(1)~=1
                disp(" (Missing datafile for the first day)");
                subFileNames(2:end+1) = subFileNames;
                subFileNames{1} = {};
            end
            if sum(diff(dayNum)~=1)>0
                error_flag = 1;
                disp(subFileNames);
                disp("Day number is not consecutive");                
                fname = subFileNames{end};
                numStart = regexp(fname,'D\d.');
                numEnd = regexp(fname,'.xls');
                endDay = str2double(fname(numStart+1:numEnd-1));
            end
            
            allDat.(AnimalLabel{k}).(task_set(t)).(subFolderNames{idx_chr(i)}) = cell(1,numel(subFileNames));
            if error_flag
                allDat.(AnimalLabel{k}).(task_set(t)).(subFolderNames{idx_chr(i)}) = cell(1,endDay);
            end
            
            %% loop through each session day: load excel data file and preprocess into dataStruct
            for j = 1:numel(allDat.(AnimalLabel{k}).(task_set(t)).(subFolderNames{idx_chr(i)}))
                display_counter(j);
                
                %fname = convertStringsToChars([dataset_address+"/"+subFolderNames{idx_chr(i)}+'/'+subFileNames{j}]);
                fname0 = AnimalLabel{k}+"_"+subFolderNames{idx_chr(i)}+"_"+task_set(t)+"_D";
                fname = convertStringsToChars([dataset_address+"/"+subFolderNames{idx_chr(i)}+"/"+fname0+j+".xlsx"]);
                
                % Edit: make sure if the day number corresponds to the indexed num
                numStart = regexp(fname,'D\d.');
                numEnd = regexp(fname,'.xls');
                if ~exist(fname,'file')
                    % set missing file as empty cell and skip
                    allDat.(AnimalLabel{k}).(task_set(t)).(subFolderNames{idx_chr(i)}){j} = {};
                    continue;
                end
                dayID = str2double(fname(numStart+1:numEnd-1));
                if j~=dayID
                    disp(fname+": Make sure if the day number ("+dayID+") corresponds to the indexing num ("+j+")"); 
                end
                
                % Load excel data & preprocess raw data
                rawData = readtable(fname,'VariableNamingRule','modify');
                block_stats = preprocessBlockData(rawData, sched_set, task_set(t), stage_labels(i), sched_labels(i), include_noChoiceTrials);
                
                allDat.(AnimalLabel{k}).(task_set(t)).(subFolderNames{idx_chr(i)}){j} = block_stats;
                
                % for shutting off annoying warning messages
                [w, MSGID] =lastwarn();
                warning('off', MSGID);
            end
        end
    end
    fprintf('\n');
end
ET = toc;
disp("Elapsed time is "+ET/60+" minutes.");

%% save preprocessed output struct
cd('../../dataset/preprocessed');

fname = "allAnimalsDat";
if numel(task_set)==1; fname = fname + "_"+task_set; end
if include_noChoiceTrials; fname = fname + "_IncludeNoChoiceTrials"; end

save([fname+".mat"],'allDat','conditions_set');
disp("Preprocessed file saved: "+fname+".mat");


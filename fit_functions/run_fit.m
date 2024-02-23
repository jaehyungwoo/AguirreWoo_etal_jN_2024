% run fit for rats dataset
clearvars; clc; close all
if contains(pwd,'f004p51')
    cd('/Users/f004p51/Dropbox (Dartmouth College)/CCNL/Alicia_Data/Claudia_SO_AO_data')
else
    cd('C:\Users\wojh1\Dropbox (Dartmouth College)\CCNL\Alicia_Data\Claudia_SO_AO_data');
end

% task_label = "SO";
task_label = "AO";

% load dataset
include_NoChoiceTrials = 1;     % dataset including uncommitted choice trials

if include_NoChoiceTrials
    dataset_label = task_label+"_IncludeNoChoice";
    load("dataset/preprocessed/allAnimalsDat_IncludeNoChoiceTrials.mat");
else
    dataset_label = task_label;
    load("dataset/preprocessed/allAnimalsDat.mat");
end

% load models
numFitRun = 50;
[models] = initialize_RL_models(dataset_label, 1, 0, numFitRun);

ModelsArray = string;
for m = 1:length(models)
    ModelsArray(m) = models{m}.name;
end
%% Run model fitting

use_parfor = 1;     % run w/ parallel computing?

animalIDs = fieldnames(allDat);
for m = 1:length(models)
    this_model = models{m}; disp("Fitting Model #"+m+". "+this_model.name+" for data: "+dataset_label);    
    
    if ~this_model.fit_exists
        tic
        % initialize
        Fit = struct;
        SessInfo = struct;
        all_sess_cnt = 0;
                
        if ~strcmp(this_model.extract_initpar_from,'none')
            extractedParModel = models{ModelsArray==this_model.extract_initpar_from};
        else
            extractedParModel = {};
        end
        
        for i = 1:numel(animalIDs)
            animal_label = animalIDs{i};    disp(i+"/"+numel(animalIDs)+". "+animal_label);
            thisDat = allDat.(animal_label).(task_label);
            if isempty(thisDat)
               continue;
            else
                %% fit this animal data
                all_phases = fieldnames(thisDat);
                for j = 1:numel(all_phases)
                    thisPhase = all_phases{j};
                    for k = 1:numel(thisDat.(thisPhase))
                        thisSess = thisDat.(thisPhase){k};
                        all_sess_cnt = all_sess_cnt + 1;
                        
                        % compile choice data to fit
                        stats = struct;
                        stats.r = thisSess.r;
                        if strcmp(task_label,"AO")
                            stats.c = thisSess.cloc;
                        elseif strcmp(task_label,"SO")
                            stats = thisSess.cstim;
                        end
                        
                        if use_parfor
                            allSessStats{all_sess_cnt} = stats;
                        else
                            % fit using random initpar
                            minBIC = intmax;
                            for n = 1:numFitRun
                                if ~isempty(extractedParModel)&&n>numFitRun/2
                                    this_model.initparpool{n}(1:numel(extractedParModel.initpar)) = extractedParModel.Fit.fitpar{all_sess_cnt};
                                end
                                [qpar, negloglike, bic, nlike, aic] = ...
                                    fit_fun(stats, this_model.fun, this_model.initparpool{n},this_model.lb,this_model.ub);
                                if bic<minBIC
                                    minBIC = bic;
                                    fitpar0 = qpar;
                                    LL0 = negloglike;
                                    nlike0 = nlike;
                                    aic0 = aic;
                                end
                            end 
                            Fit.fitpar{all_sess_cnt} = fitpar0;
                            Fit.LL(all_sess_cnt) = LL0;
                            Fit.AIC(all_sess_cnt) = aic0;
                            Fit.NLike(all_sess_cnt) = nlike0;
                        end
                        
                        % compile subject & session info
                        SessInfo.numTrials(all_sess_cnt) = sum(~isnan(stats.r));
                        SessInfo.animal_ids(all_sess_cnt) = convertCharsToStrings(animal_label);
                        SessInfo.sex(all_sess_cnt) = convertCharsToStrings(allDat.(animal_label).sex);
                        
                        if strcmp(thisPhase,'D')
                            SessInfo.phase(all_sess_cnt) = convertCharsToStrings(thisPhase);
                        else
                            if contains(thisPhase, allDat.(animal_label).R1)
                                SessInfo.phase(all_sess_cnt) = "R1";
                            elseif contains(thisPhase, allDat.(animal_label).R2)
                                SessInfo.phase(all_sess_cnt) = "R2";
                            else
                               error(all_sess_cnt+": Set correct phase label."); 
                            end
                        end
                        
                        if thisSess.prob1000
                            SessInfo.schedules(all_sess_cnt) = "prob1000";
                        elseif thisSess.prob9010
                            SessInfo.schedules(all_sess_cnt) = "prob9010";
                        elseif thisSess.prob8020
                            SessInfo.schedules(all_sess_cnt) = "prob8020";
                        elseif thisSess.prob7030
                            SessInfo.schedules(all_sess_cnt) = "prob7030";
                        end
                        SessInfo.condition(all_sess_cnt) = convertCharsToStrings(allDat.(animal_label).condition);
                        SessInfo.drugOrder(all_sess_cnt) = convertCharsToStrings(allDat.(animal_label).R1)+"_first";
                        
                        if mod(all_sess_cnt,100)==0; disp(all_sess_cnt); end   
                    end
                end
            end
        end
        
        % if fitting in parallel
        if use_parfor
            Fitpar = cell(1,all_sess_cnt);    
            LL = nan(1,all_sess_cnt);    AIC = LL;   NLike = LL; 
            parfor i = 1:all_sess_cnt
                InitPool = this_model.initparpool;
                stats = allSessStats{i};
                minBIC = intmax;
                fitpar0 = [];   LL0 = NaN;  aic0 = NaN; nlike0 = NaN;
                for n = 1:numFitRun 
                    if ~isempty(extractedParModel)&&n>numFitRun/2
                        InitPool{n}(1:numel(extractedParModel.initpar)) = extractedParModel.Fit.fitpar{i};
                    end
                    [qpar, negloglike, bic, nlike, aic] = ...
                        fit_fun(stats, this_model.fun,this_model.initparpool{n},this_model.lb,this_model.ub);
                    if bic<minBIC
                        minBIC = bic;
                        fitpar0 = qpar;
                        LL0 = negloglike;
                        nlike0 = nlike;
                        aic0 = aic;
                    end
                end 
                Fitpar{i} = fitpar0;
                LL(i) = LL0;
                AIC(i) = aic0;
                NLike(i) = nlike0;                
                if mod(i,100)==0; disp(i); end
            end
            Fit.fitpar = Fitpar;
            Fit.LL = LL;
            Fit.AIC = AIC;
            Fit.NLike = NLike;
        end
        this_model.Fit = Fit;
        this_model.Fit.numFitRun = numFitRun;
        this_model.SessInfo = SessInfo;        
        
        % save output
        model_struct = this_model;
        fname = "output/model/fitNew/"+dataset_label+"/"+this_model.name+".mat";
        save(fname, 'model_struct');    disp("File saved!");
        ET = toc;
        disp("Elaspsed time is "+ET/60+" minutes");
        
        fprintf('\n\n');
        models{m} = this_model;
    end
end

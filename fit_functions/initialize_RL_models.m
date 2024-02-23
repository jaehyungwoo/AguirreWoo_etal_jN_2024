function [models] = initialize_RL_models(dataset_label, load_fit, load_fitN, numfitruns)
    if ~exist("load_fit",'var')
        load_fit = true;
    end
    if ~exist("load_fitN",'var')
        load_fitN = true;
    end
    if ~exist("numfitruns",'var')
        numfitruns = 10;
    end

    single_alpha = 1;   
    Side_bias = 0;
    separate_alpha = 0;
    decay_chosen = 0;
    hybrid_mod = 0;

models = {};
%% RL2: Single-component
if single_alpha
    m = length(models) + 1;
    models{m}.name = 'RL1_nondecay';    
    models{m}.fun = 'RL1_nondecay';     
    models{m}.initpar=[.5  5];   
    models{m}.lb     =[ 0  0];   
    models{m}.ub     =[ 1 100];
    models{m}.label = "RL1";
    models{m}.plabels = ["\alpha", "\beta"];
    models{m}.extract_initpar_from = 'none';

    m = length(models) + 1;
    models{m}.name = 'RL1_decay';    
    models{m}.fun = 'RL1_decay';     
    models{m}.initpar=[.5  5  .5];   
    models{m}.lb     =[ 0  0   0];   
    models{m}.ub     =[ 1 100  1];
    models{m}.label = "RL1_{decay}";
    models{m}.plabels = ["\alpha", "\beta","\gamma_{\it{d}}"];
    models{m}.extract_initpar_from = 'RL1_nondecay';

    if Side_bias
        m = length(models) + 1;
        models{m}.name = 'RL15_nondecay';    
        models{m}.fun = 'RL15_nondecay';     
        models{m}.initpar=[.5   5  0];   
        models{m}.lb     =[ 0   0 -1];   
        models{m}.ub     =[ 1 100  1];
        models{m}.label = "RL1+\beta_0";
        models{m}.plabels = ["\alpha", "\beta", "\beta_0"];
        models{m}.extract_initpar_from = 'RL1_nondecay';
        
        m = length(models) + 1;
        models{m}.name = 'RL15_decay';    
        models{m}.fun = 'RL15_decay';     
        models{m}.initpar=[.5  5  .5  0];   
        models{m}.lb     =[ 0  0   0 -1];   
        models{m}.ub     =[ 1 100  1  1];
        models{m}.label = "RL1_{decay}+\beta_0";
        models{m}.plabels = ["\alpha", "\beta_V","\gamma_{\it{d}}", "\beta_0"];
        models{m}.extract_initpar_from = 'RL1_decay';
    end    


    if decay_chosen
        m = length(models) + 1;
        models{m}.name = 'RL1_decayChosen';    
        models{m}.fun = 'RL1_decayChosen';     
        models{m}.initpar=[.5  5  .5  1];   
        models{m}.lb     =[ 0  0   0  0];   
        models{m}.ub     =[ 1 100  1  1];
        models{m}.label = "RL1_{decay+chosen}";
        models{m}.plabels = ["\alpha", "\beta","\gamma_{decay}","\gamma_C"];
        models{m}.extract_initpar_from = 'RL1_decay';
        
        m = length(models) + 1;
        models{m}.name = 'RL1_decayChosen05';    
        models{m}.fun = 'RL1_decayChosen05';     
        models{m}.initpar=[.5  5  .5  1];   
        models{m}.lb     =[ 0  0   0  0];   
        models{m}.ub     =[ 1 100  1  1];
        models{m}.label = "RL1_{decay+chosen\rightarrow0.5}";
        models{m}.plabels = ["\alpha", "\beta","\gamma_{decay}","\gamma_{C,0.5}"];
        models{m}.extract_initpar_from = 'RL1_decay';
    end
end

if separate_alpha
    m = length(models) + 1;
    models{m}.name = 'RL2_nondecay';    
    models{m}.fun = 'RL2_nondecay';     
    models{m}.initpar=[.5  5  .5];   
    models{m}.lb     =[ 0  0   0];   
    models{m}.ub     =[ 1 100  1];
    models{m}.label = "RL2";
    models{m}.plabels = ["\alpha_{+}", "\beta", "\alpha_{-}"];
    models{m}.extract_initpar_from = 'RL1_nondecay';

    m = length(models) + 1;
    models{m}.name = 'RL2_decay';     
    models{m}.fun = 'RL2_decay';     
    models{m}.initpar=[.5  5  .5 .5];   
    models{m}.lb     =[ 0  0   0  0];   
    models{m}.ub     =[ 1 100  1  1];
    models{m}.label = "RL2_{decay}";
    models{m}.plabels = ["\alpha_{+}", "\beta", "\alpha_{-}","\gamma_{decay}"];
    models{m}.extract_initpar_from = 'RL2_nondecay';
end

if hybrid_mod
    m = length(models) + 1;
    models{m}.name = 'Hyb_RL2_decay';     
    models{m}.fun = 'Hybrid_RL_decay';     
    models{m}.initpar=[.5  5  .5 .5 .5];   
    models{m}.lb     =[ 0  0   0  0  0];   
    models{m}.ub     =[ 1 100  1  1  1];
    models{m}.label = "RL2-Hybrid_{decay}";
    models{m}.plabels = ["\alpha_{rew}", "\beta", "\alpha_{unrew}","decay","\omega_V"];
    models{m}.extract_initpar_from = 'RL1_nondecay';
end
%% check for errors

numfields = numel(fieldnames(models{1}));
Names_set = {};

for m = 1:length(models)
    disp(models{m}.name);

    if numel(models{m}.initpar)~=numel(models{m}.plabels)
        disp(models{m}.name);
        error('Error in number of parameters');
    end
    if numel(fieldnames(models{m}))~=numfields
        error("Check number of fields : Model "+m);
    end
    Names_set{m} = models{m}.name;    
    ffunc = models{m}.fun;
    if ~(exist(ffunc,'file'))
        error(ffunc+": corresponding model fit function does not exist")
    end

end
if numel(unique(Names_set))~=numel(models)
   error('Error: should assign unique labels to each model');
end

%% load existing output
% if the saved output file exists, flag and load the data

output_dir = "output/model/";

for m = 1:length(models)
    % fitting data
    fitfname = output_dir+"fitNew/"+dataset_label+"/"+models{m}.name+".mat";    
    if exist(fitfname, 'file') && load_fit
        load(fitfname, 'model_struct'); 
        orig_struct = model_struct;
        %models{m} = model_struct;
        models{m}.fit_exists = 1;
        field_names = fieldnames(orig_struct);
        for cnt = 1:length(field_names)
            if ~isfield(models{m}, field_names{cnt})
                models{m}.(field_names{cnt}) = orig_struct.(field_names{cnt});
            end
        end
    else
        models{m}.fit_exists = 0;
    end
    
    % Alternative fitting data: fit 2 params per session
    fitfname = output_dir+"fit2/"+dataset_label+"/"+models{m}.name+".mat";    
    if exist(fitfname, 'file') && load_fitN
        load(fitfname, 'model_struct'); 
        orig_struct = model_struct;
        models{m}.fit2_exists = 1;
        field_names = fieldnames(orig_struct);
        for cnt = 1:length(field_names)
            if ~isfield(models{m}, field_names{cnt})
                models{m}.(field_names{cnt}) = orig_struct.(field_names{cnt});
            end
        end
    else
        models{m}.fit2_exists = 0;
    end
    
    % Alternative fitting data: fit 2 params per session
    fitfname = output_dir+"fit3/"+dataset_label+"/"+models{m}.name+".mat";    
    if exist(fitfname, 'file') && load_fitN
        load(fitfname, 'model_struct'); 
        orig_struct = model_struct;
        models{m}.fit3_exists = 1;
        field_names = fieldnames(orig_struct);
        for cnt = 1:length(field_names)
            if ~isfield(models{m}, field_names{cnt})
                models{m}.(field_names{cnt}) = orig_struct.(field_names{cnt});
            end
        end
    else
        models{m}.fit3_exists = 0;
    end
    
    
    % cross-validation data
    cvfname = output_dir+"cv/"+dataset_label+"/"+models{m}.name+".mat";  
    if exist(cvfname, 'file')
        models{m}.cv_exists = 1;
        load(cvfname, 'model_struct');
        field_names = fieldnames(model_struct);
        for cnt = 1:length(field_names)
            if ~isfield(models{m}, field_names{cnt})
                models{m}.(field_names{cnt}) = model_struct.(field_names{cnt});
            end
        end
    else
        models{m}.cv_exists = 0;
    end
    
end

%% generate initial value pool for evenly spaced search space
for m = 1:length(models)
    initparpool = nan(numfitruns-1,length(models{m}.initpar));     
    
    for p = 1:length(models{m}.initpar)
        initpars = linspace(models{m}.lb(p),models{m}.ub(p),numfitruns-1);
        initpars = initpars(randperm(length(initpars)));    %randomize order
        initparpool(:,p) = initpars;
    end
    
    models{m}.initparpool = mat2cell(initparpool,ones(numfitruns-1,1),length(models{m}.initpar));
    models{m}.initparpool{numfitruns} = models{m}.initpar;
end

%% display info
disp(">> "+ numel(models)+" models initialized:"); 
disp('-----------------');

end
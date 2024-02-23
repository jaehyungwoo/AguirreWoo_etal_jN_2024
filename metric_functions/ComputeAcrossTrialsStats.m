% Compute cumulative stats (or runnig window stats) across trials within a block
function [outputCum, outputRun, output_idx] = ComputeAcrossTrialsStats(block_stats, animalInfo, beforeAfterN, trialsType, prevBlock_stats)
    numTrial = length(block_stats.r);
    
    outputCum = struct;   % cumulative stats
    outputRun = struct;   % running window stats
    outputRun.N = 10;     % running window size N
    output_idx = struct;    % index for schedule label, phase, injection type
    
    % initialize output bins
    beh_metrics = ["pwin","pbetter","pstay","winstay"',"loseswitch","RI_BW","RI_B","RI_W"];
    ent_metrics = ["ERDS","EODS","ERODS","MIRS","MIOS","MIORS","H_str"];

    varNames = [beh_metrics, ent_metrics];
    varNames_Loc = ["pstay","winstay"',"loseswitch","ERDS","H_str"]+"_Loc";
    varNames_Loc = [varNames_Loc, "RI_LR", "RI_L", "RI_R"];
    
    for v = 1:numel(varNames)
        outputCum.(varNames(v)) = nan(1,numTrial);      % cumulative ERDS using all previous trials        
        outputRun.(varNames(v)) = nan(1,numTrial);      % ERDS from a running window of arbitrary size
    end    
    
    if strcmp(block_stats.task,"SO")
        for v = 1:numel(varNames_Loc)
            outputCum.(varNames_Loc{v}) = nan(1,numTrial);      % location-based ERDS for SO task
            outputRun.(varNames_Loc{v}) = nan(1,numTrial);
        end
    end

    %% Behavioral metrics
    
    shortL_flag = 0;
    if isempty(prevBlock_stats)||strcmp(trialsType,"beforeRev")
        combined_stats = block_stats;
        t_offset = 0;    
    elseif strcmp(trialsType,"firstN")
        t_offset = outputRun.N;
        combined_stats = struct;        % combined info from previous block before reversal
        if ~(length(prevBlock_stats.r)<t_offset)
            combined_stats.r = [prevBlock_stats.r(end-t_offset+1:end); block_stats.r];
            combined_stats.cloc = [prevBlock_stats.cloc(end-t_offset+1:end); block_stats.cloc];
            combined_stats.hr_side = [prevBlock_stats.hr_side(end-t_offset+1:end); block_stats.hr_side];
            if isfield(block_stats,'cstim')
                combined_stats.cstim = [prevBlock_stats.cstim(end-t_offset+1:end); block_stats.cstim];
                combined_stats.hr_stim = [prevBlock_stats.hr_stim(end-t_offset+1:end); block_stats.hr_stim];
            end
        else
            numFill = t_offset - length(prevBlock_stats.r);
            combined_stats.r = [nan(numFill,1); prevBlock_stats.r; block_stats.r];
            combined_stats.cloc = [nan(numFill,1); prevBlock_stats.cloc; block_stats.cloc];
            combined_stats.hr_side = [nan(numFill,1); prevBlock_stats.hr_side; block_stats.hr_side];
            if isfield(block_stats,'cstim')
                combined_stats.cstim = [nan(numFill,1); prevBlock_stats.cstim; block_stats.cstim];
                combined_stats.hr_stim = [nan(numFill,1); prevBlock_stats.hr_stim; block_stats.hr_stim];
            end
            shortL_flag = 1;
            nan_idx = isnan(combined_stats.r);
        end
    end
    
    if strcmp(block_stats.task,"SO")
        choice_vec = combined_stats.cstim;     % SO task
    else
        choice_vec = combined_stats.cloc;      % AO task
    end
    
    % trial stats vectors
    choseBetter = (combined_stats.cloc==combined_stats.hr_side);
    stay = choice_vec(1:end-1)==choice_vec(2:end);
    stayVec = [NaN; stay];
    if shortL_flag; combined_stats.r(nan_idx) = 0; end
    WinStayVec = [NaN; stay&combined_stats.r(1:end-1)];
    if shortL_flag; combined_stats.r(nan_idx) = 1; end
    LoseSwitchVec = [NaN; ~stay&~combined_stats.r(1:end-1)];
    BetterStay = [NaN; stay&choseBetter(1:end-1)];
    WorseStay = [NaN; stay&~choseBetter(1:end-1)];
    if strcmp(block_stats.task,"SO")
        choseLeft = (combined_stats.cloc==-1);
        stayLoc = combined_stats.cloc(1:end-1)==combined_stats.cloc(2:end);
        stayLocVec = [NaN; stayLoc];
        LeftStay = [NaN; stayLoc&choseLeft(1:end-1)];
        RightStay = [NaN; stayLoc&~choseLeft(1:end-1)];
    end
    
    % 1. compute cumulative stats: behavioral
    for t = 1:numTrial
        cumWindow_idx = [1+t_offset:t+t_offset];
        outputCum.pwin(t) = mean(combined_stats.r(cumWindow_idx));
        outputCum.pbetter(t) = mean(choseBetter(cumWindow_idx));
        outputCum.pstay(t) = mean(stayVec(cumWindow_idx),'omitnan');
            pBetterStay(t) = mean(BetterStay(cumWindow_idx),'omitnan');
            pWorseStay(t) = mean(WorseStay(cumWindow_idx),'omitnan');
        outputCum.winstay(t) = mean(WinStayVec(cumWindow_idx),'omitnan')/mean(combined_stats.r(cumWindow_idx));
        outputCum.loseswitch(t) = mean(LoseSwitchVec(cumWindow_idx),'omitnan')/mean(~combined_stats.r(cumWindow_idx));
        if strcmp(block_stats.task,"SO")
            outputCum.pstay_Loc(t) = mean(stayLocVec(cumWindow_idx),'omitnan');
            pChooseLeft(t) = mean(choseLeft(cumWindow_idx),'omitnan');
            pLeftStay(t) = mean(LeftStay(cumWindow_idx),'omitnan');
            pRightStay(t) = mean(RightStay(cumWindow_idx),'omitnan');
        end
    end
    outputCum.RI_BW = outputCum.pstay - (outputCum.pbetter.^2+(1-outputCum.pbetter).^2);
    outputCum.RI_B = pBetterStay - outputCum.pbetter.^2;
    outputCum.RI_W = pWorseStay - (1-outputCum.pbetter).^2;  
    if strcmp(block_stats.task,"SO")
        outputCum.RI_LR = outputCum.pstay_Loc - (pChooseLeft.^2+(1-pChooseLeft).^2);
        outputCum.RI_L = pLeftStay - pChooseLeft.^2;
        outputCum.RI_R = pRightStay - (1-pChooseLeft).^2;  
    end
    
    % 2. running window stats: behavioral
    % for after reversal stats, use carry over information from previous block
    for t = 1:numTrial
        runWindow_idx = [max(t+t_offset-outputRun.N+1,1):t+t_offset];
        outputRun.pwin(t) = mean(combined_stats.r(runWindow_idx));
        outputRun.pbetter(t) = mean(choseBetter(runWindow_idx),'omitnan');
        outputRun.pstay(t) = mean(stayVec(runWindow_idx),'omitnan');
            pBetterStay(t) = mean(BetterStay(runWindow_idx),'omitnan');
            pWorseStay(t) = mean(WorseStay(runWindow_idx),'omitnan');
        outputRun.winstay(t) = mean(WinStayVec(runWindow_idx),'omitnan')/mean(combined_stats.r(runWindow_idx));
        outputRun.loseswitch(t) = mean(LoseSwitchVec(runWindow_idx),'omitnan')/mean(~combined_stats.r(runWindow_idx));
        if strcmp(block_stats.task,"SO")
            outputRun.pstay_Loc(t) = mean(stayLocVec(runWindow_idx),'omitnan');
            pChooseLeft(t) = mean(choseLeft(runWindow_idx),'omitnan');
            pLeftStay(t) = mean(LeftStay(runWindow_idx),'omitnan');
            pRightStay(t) = mean(RightStay(runWindow_idx),'omitnan');
        end
    end
    outputRun.RI_BW = outputRun.pstay - (outputRun.pbetter.^2+(1-outputRun.pbetter).^2);
    outputRun.RI_B = pBetterStay - outputRun.pbetter.^2;
    outputRun.RI_W = pWorseStay - (1-outputRun.pbetter).^2;
    if strcmp(block_stats.task,"SO")
        outputRun.RI_LR = outputRun.pstay_Loc - (pChooseLeft.^2+(1-pChooseLeft).^2);
        outputRun.RI_L = pLeftStay - pChooseLeft.^2;
        outputRun.RI_R = pRightStay - (1-pChooseLeft).^2;  
    end
    
    %% entropy metrics

    if strcmp(block_stats.task,"SO")
        choice = combined_stats.cstim;
        choiceLoc = combined_stats.cloc;
        strLoc = choiceLoc(1:end-1)==choiceLoc(2:end);
        hr_opt = combined_stats.hr_stim;
    else
        choice = combined_stats.cloc;
        hr_opt = combined_stats.hr_side;
    end
    str = choice(1:end-1)==choice(2:end);
    rew = combined_stats.r(1:end-1);
    opt = choice==hr_opt;
    
    % 1. cumulative stats: entropy-based
    for t = 1:numTrial
        cumWindow_idx = [1+t_offset:t+t_offset-1];
        CumStr = str(cumWindow_idx);
        CumRew = rew(cumWindow_idx);
        CumOpt = opt(cumWindow_idx);        
        CumRewOpt = binary_to_decimal([CumRew, CumOpt]);
        
        Eoutput = copy_field_names(struct, ...
                    {Conditional_Entropy(CumStr, CumRew, "ERDS"), ...
                     Mutual_Information(CumStr, CumRew, "MIRS"), ...
                     Conditional_Entropy(CumStr, CumOpt, "EODS"), ...
                     Mutual_Information(CumStr, CumOpt, "MIOS"), ...
                     Conditional_Entropy(CumStr, CumRewOpt, "ERODS"), ...
                     Mutual_Information(CumStr, CumRewOpt, "MIROS"), ...
                    });
        Eoutput.H_str = Shannon_Entropy(CumStr);
        if strcmp(block_stats.task,"SO")
            CumStrLoc = strLoc(cumWindow_idx);
            Eoutput = copy_field_names(Eoutput, {Conditional_Entropy(CumStrLoc, CumRew, "ERDS_Loc"),...
                                                 Mutual_Information(CumStrLoc, CumRew, "MIRS_Loc") });
            Eoutput.H_str_Loc = Shannon_Entropy(CumStrLoc);
        end
        Eset = fieldnames(Eoutput);
        for e = 1:numel(Eset)
            outputCum.(Eset{e})(t) = Eoutput.(Eset{e});
        end
    end
    
    % 2. running window stats: entropy-based
    for t = 1:numTrial
        runWindow_idx = [max(t+t_offset-outputRun.N,1):t+t_offset-1];
        RunStr = str(runWindow_idx);
        RunRew = rew(runWindow_idx);
        RunOpt = opt(runWindow_idx);
        RunRewOpt = binary_to_decimal([RunRew, RunOpt]);
        Eoutput = copy_field_names(struct, ...
                {Conditional_Entropy(RunStr, RunRew, "ERDS"), ...
                 Mutual_Information(RunStr, RunRew, "MIRS"), ...
                 Conditional_Entropy(RunStr, RunOpt, "EODS"), ...
                 Mutual_Information(RunStr, RunOpt, "MIOS"), ...
                 Conditional_Entropy(RunStr, RunRewOpt, "ERODS"), ...
                 Mutual_Information(RunStr, RunRewOpt, "MIROS"), ...
                });
        Eoutput.H_str = Shannon_Entropy(RunStr);
        if strcmp(block_stats.task,"SO")
            RunStrLoc = strLoc(runWindow_idx);
            Eoutput = copy_field_names(Eoutput, {Conditional_Entropy(RunStrLoc, RunRew, "ERDS_Loc"),...
                                                 Mutual_Information(RunStrLoc, RunRew, "MIRS_Loc") });
            Eoutput.H_str_Loc = Shannon_Entropy(RunStrLoc);
        end
        Eset = fieldnames(Eoutput);
        for e = 1:numel(Eset)
            outputRun.(Eset{e})(t) = Eoutput.(Eset{e});
        end
    end
    
    %% First/Last N trials data
    varNames = string(fieldnames(outputCum));

    if numTrial>=beforeAfterN
        if strcmp(trialsType,"firstN")
            for v = 1:numel(varNames)
                outputCum.(varNames(v)) = outputCum.(varNames(v))(1:beforeAfterN);      % cumulative ERDS using all previous trials        
                outputRun.(varNames(v)) = outputRun.(varNames(v))(1:beforeAfterN);      % ERDS from a running window of arbitrary size
            end    
        elseif strcmp(trialsType,"beforeRev")
            for v = 1:numel(varNames)
                outputCum.(varNames(v)) = outputCum.(varNames(v))(end-beforeAfterN+1:end);      % cumulative ERDS using all previous trials        
                outputRun.(varNames(v)) = outputRun.(varNames(v))(end-beforeAfterN+1:end);      % ERDS from a running window of arbitrary size
            end    
        end
    else
        % if # of trials are less than N, fill in with NaN
        numFill = beforeAfterN - numTrial;
        if strcmp(trialsType,"firstN")
            for v = 1:numel(varNames)
                outputCum.(varNames(v)) = [outputCum.(varNames(v)), nan(1,numFill)];
                outputRun.(varNames(v)) = [outputRun.(varNames(v)), nan(1,numFill)];
            end
        elseif strcmp(trialsType,"beforeRev")
            for v = 1:numel(varNames)
                outputCum.(varNames(v)) = [nan(1,numFill), outputCum.(varNames(v))];
                outputRun.(varNames(v)) = [nan(1,numFill), outputRun.(varNames(v))];
            end
        end
    end

    %% for stats after reversal, should use information from previous block before reversal
    if strcmp(trialsType,"firstN")
        if ~(strcmp(block_stats.stage,"D")&&block_stats.prob1000)
            if isempty(prevBlock_stats)
                error("Previous block before reversal info is missing"); 
            end
        end
    end
    %% mark schedule label, phase, injection type, and subject name
    schedule_subsets = ["prob1000", "prob9010", "prob8020", "prob7030"];
    stage_subsets = ["D","R1","R2"];
    injection_subsets = ["none","CNO","VEH"];
    
    for s = 1:numel(schedule_subsets)
        output_idx.(schedule_subsets(s)) = block_stats.(schedule_subsets(s));
    end
    for s = 1:3
        output_idx.(stage_subsets(s)) = 0;
        output_idx.(injection_subsets(s)) = 0;
    end
    output_idx.(block_stats.stage) = 1;
    
    switch block_stats.stage
        case "D"
            output_idx.none = 1;
        case "R1"
            output_idx.(animalInfo.R1) = 1;
        case "R2"
            output_idx.(animalInfo.R2) = 1;
    end
    output_idx.sex = string(animalInfo.sex);
    
end

function dec_array = bin_to_dec(binary_mtx)
    dec_array = zeros(size(binary_mtx, 1), 1);
    for i=1:size(binary_mtx,1)
        dec = 0;
        len = size(binary_mtx, 2);
        for p = 1:len
            dec = dec + binary_mtx(i,p)*(2^(len-p));
        end
        dec_array(i) = dec;
    end
end
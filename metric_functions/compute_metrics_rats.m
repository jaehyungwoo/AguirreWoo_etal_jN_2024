function allOutput = compute_metrics_rats(choice, reward, better_opt)
% % behavioral_metrics %
%PURPOSE:   Compute behavioral & entropy metrics
%AUTHORS:   Jae Hyung Woo 1/16/2023
%
%INPUT ARGUMENTS
%   choice:   choice vector for block/session where -1= choose left, 1=
%       choose right. choice should not include NaN trials.
%   reward:   reward vector for block/session where 1 = reward, 
%       0 = no reward
%   better_opt:  vector of "better" option (higher reward probability) in each 
%       trial. Should be the same length as choice and reward vectors. 

%OUTPUT ARGUMENTS
%   allOutput: range of behavioral metrics, see code

    % initialize 
    allOutput = struct;
    choseBetter = choice==better_opt;

    % N-1 length vectors:
    stay = choice(1:end-1)==choice(2:end);
    prevRew = reward(1:end-1);
    prevOpt = choseBetter(1:end-1);
    rew_and_opt = bin_to_dec([prevRew, prevOpt]);

    % behavioral metrics
    allOutput.pwin = mean(reward);
    allOutput.pbetter = mean(choseBetter);
    allOutput.pstay = mean(stay);
    allOutput.winstay = mean(stay&prevRew)/mean(prevRew);
    allOutput.loseswitch = mean(~stay&~prevRew)/mean(~prevRew);
    allOutput.betterstay = mean(stay&prevOpt)/mean(prevOpt);
    allOutput.RI_BW = allOutput.pstay - (allOutput.pbetter^2+(1-allOutput.pbetter)^2);
    allOutput.RI_B = mean(stay&prevOpt)- allOutput.pbetter^2;
    allOutput.RI_W = mean(stay&~prevOpt)- (1-allOutput.pbetter)^2;
    
    % entropy metrics
    allOutput.H_str = Shannon_Entropy(stay);
    allOutput = append_to_fields(allOutput, ...
        {Conditional_Entropy(stay, prevRew, "ERDS"), ...
         Mutual_Information(stay, prevRew, "MIRS"), ...
         Conditional_Entropy(stay, prevOpt, "EODS"), ...
         Mutual_Information(stay, prevOpt, "MIOS"), ...
         Conditional_Entropy(stay, rew_and_opt, "ERODS"), ...
         Mutual_Information(stay, rew_and_opt, "MIROS"), ...
         });
    
end


%% subfunctions

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
function output = Conditional_Entropy(y, x, metric_name, decomp_map)
% % conditional_entropy %
%PURPOSE:   Compute conditional entropy H(y|x) from two logical vectors 
%AUTHORS:   Ethan Trepka 10/05/2020; edits by Jae Hyung Woo
%
%INPUT ARGUMENTS
%   y:  vector of conditioned variable, each uq value is an event, e.g. 1 = stay, 0 = switch
%   x:  vector of conditioning variable, each value is an event 
%   metric_name: string, name of conditional entropy metric
%   decomp_map: map from "events" in x to strings to label to decompositions that corresponds with the label, 
%               e.g., {1,0 -> "win","lose"}

%OUTPUT ARGUMENTS
%   output: 
%       (metric_name)
%       (metric_name + _ + decomp_map values)
    
    % exclude NaN entries if any
    if sum(isnan(y))>0||sum(isnan(x))>0
        nanIDX = isnan(y)|isnan(x);
        x(nanIDX) = [];
        y(nanIDX) = [];
    end

    y_unique = unique(y);
    x_unique = unique(x);

    if ~exist('decomp_map', 'var')
        decompose_flag = false;
    else
        decompose_flag = true;
    end
    
    if decompose_flag
        decomp_vals = values(decomp_map);
        for i=1:decomp_map.Count
            output.(strcat(metric_name, "_", decomp_vals(i))) = NaN;
        end
    end
    x_entropies_storage = [NaN];
    
    for x_num = 1:length(x_unique)
        x_entropy = 0; %initial value of current entropy decomposition
        x_signal = x==x_unique(x_num);
        prob_x = mean(x_signal);
        
        for y_num = 1:length(y_unique)
            y_signal = y == y_unique(y_num);
            prob_x_given_y = Conditional_Probability(y_signal,x_signal);
            x_entropy = nansum_zero_helper([x_entropy, prob_x_given_y*prob_x*log2(prob_x_given_y)], 'all');
        end
        if decompose_flag
            output.(strcat(metric_name, "_", decomp_map(x_unique(x_num)))) = -x_entropy;
        end
        x_entropies_storage = [x_entropies_storage, -x_entropy];
    end

    output.(metric_name) = nansum_zero_helper(x_entropies_storage, 'all');

end



%%

function prob = Conditional_Probability(y,x)
    prob = mean(y&x)/mean(x);
end

function mySum = nansum_zero_helper(myArr, dimension)
    mySum = nansum(myArr, dimension);
    shouldBeNaN = sum(isnan(myArr), 2) == size(myArr, 2);
    mySum(shouldBeNaN) = NaN; 
end
function output = Mutual_Information(y, x, metric_name, decomp_map)
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
        for i = 1:decomp_map.Count
            %initialize entropy decompositions to NaN
            output.(strcat(metric_name, "_", decomp_vals(i))) = NaN;
        end
    end
    x_entropies_storage = [NaN];

    for x_num = 1:length(x_unique)
        x_entropy = NaN; %initial value of current entropy decomposition
        x_signal = x==x_unique(x_num);
        prob_x = mean(x_signal);
        for y_num = 1:length(y_unique)
            y_signal = y == y_unique(y_num);
            prob_y = mean(y_signal);
            prob_x_and_y = mean(y_signal&x_signal);
            x_entropy = nansum_zero_helper([x_entropy, prob_x_and_y*log2((prob_x_and_y)/(prob_x*prob_y))], 'all');
        end
        if decompose_flag
            output.(strcat(metric_name, "_", decomp_map(x_unique(x_num)))) = x_entropy;
        end
        x_entropies_storage = [x_entropies_storage, x_entropy];
    end
    output.(metric_name) = nansum_zero_helper(x_entropies_storage, 'all');
end

%%

function mySum = nansum_zero_helper(myArr, dimension)
    mySum = nansum(myArr, dimension);
    shouldBeNaN = sum(isnan(myArr), 2) == size(myArr, 2);
    mySum(shouldBeNaN) = NaN; 
end
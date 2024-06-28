function best_temp_cl = update_temp_cl(temp_cl, possible_remove_points, wrong_lab,sdm, n)
% randomly drop some centers doesn't satisfy the flexible min strategy
% until it fits the requirement
if nargin < 5
    n = length(wrong_lab);
end
minpp = 100;
best_min_sim = 0;
%min_sim_save = [];

while minpp ~= 0
    comb_list = nchoosek(possible_remove_points, n);
    [nRows, ~] = size(comb_list);

    TF_value = arrayfun(@(rowIdx) ...
        if_allin(comb_list(rowIdx, :), temp_cl, wrong_lab), 1:nRows);
    comb_list = comb_list(TF_value,:);
    [nRows, ~] = size(comb_list);
    for i = 1:nRows
        test_cl = temp_cl;
        for j = 1:length(wrong_lab)
            test_cl{wrong_lab(j)} = setdiff(temp_cl{wrong_lab(j)}, comb_list(i,:));
        end
        [~,difwb2]=score_computer(test_cl,sdm, 'flexmin');

        pp=cellfun(@(x) sum(x(:) < 0), difwb2);

        if pp == 0
            minpp = 0;
            min_sim = cellfun(@(x) min(min(sdm(x, x))), test_cl);
            sum_min_sim = sum(min_sim);

            if sum_min_sim > best_min_sim
                best_temp_cl = test_cl;
                best_min_sim = sum_min_sim;
                % min_sim_save = min_sim;
            end
        end
    end
    n = n + 1;
end
end

function TF = if_allin(x, temp_cl, wrong_lab)
n = length(wrong_lab);
TF_lab = false(1, n);

for i = 1:n
    TF_lab(i) = any(ismember(x, temp_cl{wrong_lab(i)}));
end

TF = all(TF_lab);
end

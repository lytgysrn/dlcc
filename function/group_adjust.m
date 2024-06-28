function [temp_cl] = group_adjust(temp_cl, sdm,spectral_method, initial_n)


[difwb,difwb2]=score_computer(temp_cl,sdm);
ppmax=cellfun(@(x) sum(x(:) < 0), difwb);
ppmin=cellfun(@(x) sum(x(:) < 0), difwb2);
pp_min=sum(ppmin);
% if any max within less than max between group similarity try split all
% other groups to find if any group fusion.
while sum(ppmax) > 0
    suspect_label = find(ppmax > 0);


    [~, best_split] = try_split(temp_cl, sdm, spectral_method, pp_min, suspect_label); 
    if ~isempty(best_split)
    [score_now,~] = score_computer(best_split, sdm, 'score'); 


    worst_label = find(score_now == max(score_now), 1);
    int_worst = intersect(best_split{worst_label}, cell2mat(temp_cl(suspect_label)'));
    if ~isempty(int_worst)
        temp_cl = best_split;
        temp_cl(worst_label) = [];
        [difwb,difwb2]=score_computer(temp_cl,sdm);
        ppmax=cellfun(@(x) sum(x(:) < 0), difwb);
        ppmin=cellfun(@(x) sum(x(:) < 0), difwb2);
        pp_min=sum(ppmin);
    else
        ppmax = [];
    end
    else
        ppmax=[];
    end
end

% randomly drop centers for handling group contamination
if pp_min>0
    wrong_lab = find(ppmin > 0);
    possible_remove_points = [];
    lwr = length(wrong_lab);
    for v = 1:lwr
        rcg_lab = wrong_lab(v);

        if_less_0 = difwb2{rcg_lab} < 0;

        add_points = temp_cl{rcg_lab}(if_less_0);
        possible_remove_points = [possible_remove_points; add_points];
    end

    if ~isempty(possible_remove_points)
        if initial_n == 1
            temp_cl = update_temp_cl(temp_cl, possible_remove_points, wrong_lab,sdm);
        else
            temp_cl = update_temp_cl(temp_cl, possible_remove_points, wrong_lab,sdm, initial_n);
        end
    end
end
end

function [avg_increase, best_split] = try_split(temp_cl, sdm, spectral_method, pp_cut, notcut)
%try split all groups except groups in notcut, and pick up the one increase
%the flexible threshold the most
    if nargin < 4
        pp_cut = 0;
    end
    if nargin < 5
        notcut = [];
    end
    % initial flexible threshold for each group
    k = length(temp_cl);
    initial_within_sim = cellfun(@(x) min(min(sdm(x, x))), temp_cl, 'UniformOutput', true);
    avg_increase = zeros(1, k);
    best_split = {};
    
    for i = 1:k
        if ~ismember(i, notcut)
            if length(temp_cl{i}) > 2
                tmp_g = temp_cl;
                split2 = spec_clus_withsim(sdm(temp_cl{i},temp_cl{i}), 2, spectral_method);
                split2=split2{1};
                tmp_g2 = cell(1, 2);
                for j = 1:2
                    tmp_g2{j} = temp_cl{i}(split2 == j);
                end
                tmp_g(i) = [];
                tmp_g = [tmp_g, tmp_g2];
                [~,difwb2]= score_computer(tmp_g,sdm);
                pp=sum(cellfun(@(x) sum(x(:) < 0), difwb2));
                if pp <= pp_cut
                    update_sim = zeros(1, 2);
                    for x = 1:2
                        sdm_subset = sdm(tmp_g2{x}, tmp_g2{x});
                        update_sim(x) = min(min(sdm_subset)); 
                    end
                    avg_increase(i) = mean(update_sim - initial_within_sim(i));
                    if i == find(avg_increase == max(avg_increase), 1)
                        best_split = tmp_g;
                    end
                end
            end
        end
    end
end




function [temp_cl] = group_adjust(temp_cl, sdm)
%adjust temp_cl to satisfy min strategy

check=score_computer(temp_cl,sdm);
ppmin=cellfun(@(x) nnz(x(:) < 0), check.difwb2);
pp_min=sum(ppmin);
pp_max=sum(cellfun(@(x) nnz(x(:) < 0), check.difwb));
pp_counts=pp_min+pp_max;
if (pp_counts)>0
    wrong_lab = find(ppmin > 0);
    while pp_counts ~= 0
        if pp_min>0
            lwr=length(wrong_lab);
            for j = 1:lwr
                group_points = temp_cl{wrong_lab(j)};
                sub_sim = sdm(group_points, group_points);
                original_betsim = check.min_bet{wrong_lab(j)};

                min_sim_growth = arrayfun(@(v) mean(min(sub_sim(setdiff(1:size(sub_sim, 1), v), setdiff(1:size(sub_sim, 1), v))) - original_betsim(setdiff(1:length(original_betsim), v))), 1:length(group_points));
                [~, max_idx] = max(min_sim_growth);
                temp_cl{wrong_lab(j)} = setdiff(group_points, group_points(max_idx));
            end
            % Recalculate scores
            check = score_computer(temp_cl, sdm);
            pp = cellfun(@(x) sum(x < 0), check.difwb2);
            pp_min = sum(pp);
        end


        if pp_min ~= 0
            pp_counts=pp_min;
            wrong_lab = find(pp ~= 0);
        else
            ppmax=cellfun(@(x) sum(x < 0), check.difwb);
            pp_max=sum(ppmax);
            if pp_max>0
                wrong_lab = find(ppmax > 0);
                lwr = length(wrong_lab);
                for v = 1:lwr
                    temp_cl{wrong_lab(v)} = temp_cl{wrong_lab(v)}(check.difwb{wrong_lab(v)} >= 0);
                end
                temp_cl = temp_cl(~cellfun('isempty', temp_cl));
                check = score_computer(temp_cl, sdm);
                pp = cellfun(@(x) sum(x < 0), check.difwb2);
                pp_min = sum(pp);
                if pp_min ~= 0
                    pp_counts=pp_min;
                    wrong_lab = find(pp ~= 0);
                else
                    pp_counts=0;
                end
            else
                pp_counts=0;
            end
        end
    end
end

end




function result = DBCAmerge(dbca, dm0_save)
    % Check if there are elements in dbca that are 0
    if any(dbca == 0)
        dbca(dbca == 0) = 100:1:(100 + sum(dbca == 0) - 1);
    end

    ab = unique(dbca);
    temp_K = length(ab);
    gap_bet_tk = zeros(temp_K, temp_K);
    all_pos = arrayfun(@(x) find(dbca == ab(x)), 1:temp_K, 'UniformOutput', false);

    for y = 1:(temp_K - 1)
        for u = (y + 1):temp_K
            c1 = all_pos{y};
            c2 = all_pos{u};
            max_sim = max(max(dm0_save(c2, c1)));
            gap_bet_tk(u, y) = max_sim;
            gap_bet_tk(y, u) = max_sim;
        end
    end
    gap_bet_tk(eye(temp_K) == 1) = 1;  % Setting diagonal values to 1

    % Convert matrix to vector, sort and get unique values
    gapset = sort(unique(gap_bet_tk(gap_bet_tk ~= 0)), 'descend');
    gapset = gapset(2:end);

    m = 1;
    gs = [];

    while isempty(gs)
        v = 1;
        while v <= temp_K
            if any(gapset(m) == gap_bet_tk(:, v))
                comb_lab = find(gap_bet_tk(:, v) == gapset(m));
                for x = 1:temp_K
                    gap_bet_tk(x, v) = max(gap_bet_tk(x, [v; comb_lab]));
                end
                gap_bet_tk(v, :) = gap_bet_tk(:, v);
                gap_bet_tk(comb_lab, :) = [];
                gap_bet_tk(:, comb_lab) = [];
                temp_K = temp_K - length(comb_lab);
            else
                v = v + 1;
            end
        end
        if length(gap_bet_tk) == 1
            gs = gapset(m);
        else
            m = m + 1;
        end
    end

    result.gs = gs;
    result.gapset = gapset;
end

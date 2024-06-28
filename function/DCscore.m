function result = DCscore(dm0, cl_label, out, mn)
if nargin < 4
    mn = 3;
end
[N_total,~]=size(dm0);
cl_label_no_zero = cl_label(cl_label ~= 0);
Kest = length(unique(cl_label_no_zero));
mins = zeros(1, Kest);
ss = zeros(1, Kest);
os = zeros(1, Kest);
Nw = zeros(1, Kest);

dm0_save = dm0;
dm0_sort = sort(dm0_save, 2, 'descend');

dbcaot = DBCA(dm0_save, min(dm0_sort(:, 2)));

temp_K = length(unique(dbcaot));

if temp_K == 1
    gs = min(dm0_sort(:, 2));
else
    gs_struct = DBCAmerge(dbcaot, dm0_save);
    gs = gs_struct.gs;
end

for w = 1:Kest
    lc = find(cl_label == w);
    Nw(w) = length(lc);

    if Nw(w) == 1
        ms = 1;
    elseif Nw(w) == 2
        ms = min(min(dm0_save(lc, lc)));
    else
        dm_c = dm0_save(lc, lc);
        dm_c_order = sort(dm_c, 2, 'descend');
        ms = min(dm_c_order(:, 2));
        DBCAclr = DBCA(dm_c, ms, mn);

        ldbca = length(unique(DBCAclr));

        if ldbca == 1
            DBCAclr = DBCA(dm_c, ms + 1e-10, mn);
        end

        if out
            while isequal(sort(unique(DBCAclr)), [0; 1])
                ms = min(dm_c_order(dm_c_order(:, 2) > ms, 2));
                DBCAclr = DBCA(dm_c, ms, mn);
                ldbca = length(unique(DBCAclr));
            end

            % if ldbca == 1
            %     ldbca = 2;
            % end
        end

        if ldbca ~= 1
            if out
                % have 2 or more clusters, then DBCAmerge can perform
                save_cut = min(dm_c_order(dm_c_order(:, 2) > ms, 2));
                DBCAclr = DBCA(dm_c, save_cut);
            end

            dbcams = DBCAmerge(DBCAclr, dm_c);
            ms = dbcams.gs;

            if out

                new_pos = find(dbcams.gapset == dbcams.gs) - 1;
                DBCAclr = DBCA(dm_c, dbcams.gapset(new_pos), mn);
                while isequal(sort(unique(DBCAclr)), [0; 1])
                    new_pos = new_pos - 1;
                    ms = dbcams.gapset(new_pos);
                    DBCAclr = DBCA(dm_c, ms, mn);
                end %try ms value until there is another cluster that is not 0 (outlier) or 1

                ms = dbcams.gapset(new_pos + 1); %use the previous ms value
                DBCAclr = DBCA(dm_c, ms, mn);
            end
        end
    end

    mins(w) = ms;
    dm0(eye(size(dm0)) == 1) = 0;

    if out
        outlier = find(DBCAclr == 0);
        lengtho = length(outlier);
    else
        lengtho = 0;
    end

    if lengtho > 0
        x_length = Nw(w) - lengtho;
        lc_ro = lc;
        lc_ro(outlier) = [];
        wc1 = arrayfun(@(x) dm0(lc_ro(x), lc(dm0(lc_ro(x), lc) >= ms)), 1:x_length, 'UniformOutput', false);
        o_lab = lc(outlier);
        wc_o = arrayfun(@(x) max(dm0(o_lab(x), lc_ro)), 1:lengtho);
    else
        wc1 = arrayfun(@(x) dm0(lc(x), lc(dm0(lc(x), lc) >= ms)), 1:Nw(w), 'UniformOutput', false);
        wc_o = [];
    end

    ss(w) = mean([cell2mat(wc1), wc_o]);

    wc2 = arrayfun(@(x) max(dm0(lc(x), setdiff(1:end, lc))), 1:Nw(w));
    maxbc = max(wc2);

    if maxbc >= gs
        overms_lab = find(wc2 >= gs);
        col_indices = setdiff(1:size(dm0, 2), lc);
        wc2_s = arrayfun(@(x) dm0(lc(x), col_indices(dm0(lc(x), col_indices) >= gs)), overms_lab, 'UniformOutput', false);
        oth_s = repmat(gs, 1, Nw(w) - length(wc2_s));
        os(w) = mean([cell2mat(wc2_s), oth_s]);
    else
        os(w) = gs;
    end
end

fs = ss - os;
score = sum((Nw/N_total) .* fs);

result = struct('score', score, 'fs', fs, 'iclus_s', ss, 'bclus_s', os, 'ms', mins, 'gs', gs, 'Nw', Nw);
end

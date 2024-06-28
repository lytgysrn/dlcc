function cluster = DBCA(dm, theta, mn)
    if nargin < 3
        mn = 3;
    end

    n = size(dm, 1);
    adj = cell(n, 1);

    for i = 1:n
        label = find(dm(i, :) >= theta);
        label(label == i) = [];
        adj{i} = label;
    end

    cluster = AssignCluster(adj);

    cls_n = unique(cluster);
    cls_n(cls_n == 0) = [];

    for j = 1:length(cls_n)
        pos = find(cluster == cls_n(j));
        if length(pos) < mn
            cluster(pos) = 0;
        end
    end
end

function cluster = AssignCluster(adj)
    n = length(adj);
    cluster = -1 * ones(n, 1);
    cl = 0;

    for i = 1:n
        if cluster(i) == -1
            if isempty(adj{i})
                cluster(i) = 0;
            else
                cl = cl + 1;
                cluster = cl_loop(i, cluster, adj, cl);
            end
        end
    end
end

function cluster = cl_loop(k, cluster, adj, cl)
    if cluster(k) == -1
        cluster(k) = cl;
        if ~isempty(adj{k})
            for w = 1:length(adj{k})
                new_k = adj{k}(w);
                cluster = cl_loop(new_k, cluster, adj, cl);
            end
        end
    end
end
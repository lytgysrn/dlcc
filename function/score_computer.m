
function sc = score_computer(temp_cl, sdm, need_to_compute)
%compute with-in/between group minimum/maximum similarity and related gaps or scores
if iscolumn(temp_cl) 
    temp_cl = temp_cl'; 
    temp_cl = cellfun(@(x) x(:), temp_cl, 'UniformOutput', false);
end

if nargin <3
    k = length(temp_cl);
    bml_v = cell(1,k);
    bml_vm = bml_v;

    % minimum within group sim
    wml_v = cellfun(@(x) min(sdm(x, x)), temp_cl, 'UniformOutput', false);
    % maximum within group sim
    wml_m = cellfun(@(x) ...
    (isscalar(x)) * 1 + (length(x) > 1) * max(sdm(x, x) - diag(diag(sdm(x, x)))), ...
    temp_cl, 'UniformOutput', false);


    % maximum between group sim & largest minimum between group sim
    for i = 1:k
        idx=setdiff(1:k,i);
        others = vertcat(temp_cl{idx});
   
        bml_v{i} = max(sdm(temp_cl{i}, others),[],2)';

        sim_v = zeros(k-1, length(temp_cl{i}));
        for x=1:(k-1)
            sim_v(x,:) = min(sdm(temp_cl{i}, temp_cl{idx(x)}),[],2);
        end
        bml_vm{i}=max(sim_v);
    end

    % difference within vs. between
    difwb = cellfun(@(x, y) x - y, wml_m, bml_v, 'UniformOutput', false);
    difwb2 = cellfun(@(x, y) x - y, wml_v, bml_vm, 'UniformOutput', false);
    sc = struct('difwb', {difwb}, 'difwb2', {difwb2}, 'min_within',{wml_v},'max_within',{wml_m},'max_bet',{bml_v},'min_bet',{bml_vm});
elseif strcmp(need_to_compute, 'score')
    k = length(temp_cl);
    bml_v = cell(1,k);
    % minimum within group sim
    wml_v = cellfun(@(x) min(sdm(x, x)), temp_cl, 'UniformOutput', false);
    for i = 1:k
        idx=setdiff(1:k,i);
        others = vertcat(temp_cl{idx});
        bml_v{i} = max(sdm(temp_cl{i}, others),[],2)';
    end
    t1 = sqrt(cellfun(@sum, wml_v));
    t2 = cellfun(@mean, bml_v);
    sc = t2./t1;


elseif strcmp(need_to_compute, 'flexmin')
    k = length(temp_cl);
    bml_vm = cell(1,k);
    % minimum within group sim
    wml_v = cellfun(@(x) min(sdm(x, x)), temp_cl, 'UniformOutput', false);
    for i = 1:k
        idx=setdiff(1:k,i);
        sim_v = zeros(k-1, length(temp_cl{i}));
        for x=1:(k-1)
            sim_v(x,:) = min(sdm(temp_cl{i}, temp_cl{idx(x)}),[],2);
        end
        bml_vm{i}=max(sim_v);
    end
    difwb2 = cellfun(@(x, y) x - y, wml_v, bml_vm, 'UniformOutput', false);

    sc = struct('difwb2', {difwb2}, 'min_within',{wml_v},'min_bet',{bml_vm});
elseif strcmp(need_to_compute, 'mw')
        % minimum within group sim
     wml_v = cellfun(@(x) min(sdm(x, x)), temp_cl, 'UniformOutput', false);
     sc = wml_v;
end
end
% sc = struct('score', difwb, 'wml', wml_v, 'bml', bml_v, 'bml2', bml_vm, 'bc', difwb2);

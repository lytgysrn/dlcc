function [remove_point, sub_tcl] = max_sim_remove(tcl, sim_matrix, to_remove)
    
    sub_tcl = tcl;
    sub_tcl(to_remove) = [];
    remove_point=tcl{to_remove};
    if length(remove_point)>1
   
    th_bar = max(cellfun(@(x) ...
        max(max(sim_matrix(x, x) - diag(diag(sim_matrix(x, x))))), ...
        tcl(to_remove)));

 
    other_points = unique([sub_tcl{:}]);

    to_remove_sim = max(sim_matrix(other_points, tcl{to_remove}), [], 2);

    new_sim = sim_matrix(other_points, other_points);
    expanded_sim = ones(size(new_sim, 1) + 1, size(new_sim, 2) + 1);
    expanded_sim(2:end, 2:end) = new_sim;
    expanded_sim(1, 2:end) = to_remove_sim';
    expanded_sim(2:end, 1) = to_remove_sim;

    %Find other centers too close to the dropping group
    max_sim_group = DBCA(expanded_sim, th_bar, 1);

    % update tcl
    if max_sim_group(1) == 1
        add_remove = other_points(max_sim_group(2:end) == 1);
        sub_tcl = cellfun(@(group) setdiff(group, add_remove), sub_tcl, 'UniformOutput', false);
    else
        add_remove = [];
    end
       remove_point=[remove_point, add_remove];
    end


end
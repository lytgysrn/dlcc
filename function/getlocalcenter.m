function [rank_mat,rank_topobs,dm0_order] = getlocalcenter(data,dm0,s,depth,Lmatrix)
%find the deepest point for each subset, named as local centers.

%Input:
% data: data set
% dm0: depth-based similarity matrix
% s: size of each subset
% depth: type of depth (default: spatial depth)
% Lmatrix (optional): norm matrix, storing the lengths of all vectors between points. If
% provided, spatial depth can be computed more efficiently.

%Output:
  %rank_topobs: depth center of each subset
  %rank_mat: the depth rank of the obs in its own subset.
  %dm0_order: the ordering matrix of dm0. Each column i stores the labels of points, starting from the nearest point to the farthest point from point i.

if (nargin<4)
    depth='spatial';
end

Nobs = size(data,1);
[~,dm0_order] = sort(dm0,2,'descend');
dm0_order=dm0_order';
rank_mat = zeros(Nobs,1);
rank_topobs = zeros(1,Nobs);

% parfor_progress(Nobs/10);
for j = 1:Nobs
    idx=dm0_order(1:s,j);
    if strcmp(depth,'spatial')
        if nargin<5
            local_d = spatial_d(data(idx,:),data(idx,:));
        else
            local_d = spatial_d(data(idx,:),data(idx,:),Lmatrix(idx,idx));
        end
    elseif strcmp(depth,'mahalanobis')
        local_d=Maha_d(data(idx,:),data(idx,:),'MCD');
    end
    names = dm0_order(1:s,j);
    local_rank = tiedrank(-local_d,'min');
    rank_mat(j) = local_rank(1);
    rank_topobs(j) = names(local_rank==1);
    % if mod(j,10)==0
    % parfor_progress;
    % end
end

end

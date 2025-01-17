function result = left_class(X, temp_clus, method, varargin)
%classify left observations based on the information in temporary clusters
%Input:
  %X: the whole data set
  %temp_clus: temporary clustering results
  %method: classification methods (maxdep, rf, knn)
%Output
  %final clustering results

% Create an input parser object
p = inputParser;

% Add required arguments
addRequired(p, 'X');
addRequired(p, 'temp_clus');
addRequired(p, 'method');

% Define default values for optional arguments
defaultMaxdepth = false;
defaultDepth = 'spatial';
defaultLeafsize= 0;
defaultNtrees = 100;
defaultK_knn = 5;
defaultDm0 = [];

% Add optional arguments to the parser
addParameter(p, 'maxdepth', defaultMaxdepth);
addParameter(p, 'depth', defaultDepth);
addParameter(p, 'ntrees', defaultNtrees);
addParameter(p, 'leaf_size', defaultLeafsize);
addParameter(p, 'K_knn', defaultK_knn);
addParameter(p, 'dm0', defaultDm0);

% Parse input arguments
parse(p, X, temp_clus, method, varargin{:});

% Extract parsed data from the parser object
X = p.Results.X;
temp_clus = p.Results.temp_clus;
method = p.Results.method;
maxdepth = p.Results.maxdepth;
depth = p.Results.depth;
ntrees = p.Results.ntrees;
K_knn = p.Results.K_knn;
dm0 = p.Results.dm0;
leaf_size=p.Results.leaf_size;

% begin the function
Kclus = length(temp_clus);
labelledobs = cell2mat(temp_clus);
N_labelled=length(labelledobs);
left_obs_label = find(~ismember(1:size(X,1),labelledobs));
left_obs = X(left_obs_label,:);

if strcmp(method, 'maxdep')

    if strcmp(depth, 'mahalanobis')
        depth_mat = arrayfun(@(c) Maha_d(left_obs, X(temp_clus{c}, :),'moment'), 1:Kclus, 'UniformOutput', false);
        depth_mat = horzcat(depth_mat{:});
    elseif strcmp(depth, 'spatial')
        depth_mat = arrayfun(@(c) spatial_d(left_obs, X(temp_clus{c}, :)), 1:Kclus, 'UniformOutput', false);
        depth_mat = horzcat(depth_mat{:});
    end
    [~, depth_order] = max(depth_mat, [], 2);

    cluster = cell(1, Kclus);
    for i = 1:Kclus
        obs_label=left_obs_label(depth_order==i);
        cluster{i} = [temp_clus{i}, obs_label];
    end
    cluster_vector = zeros(size(X, 1), 1);
    for i = 1:length(cluster)
        cluster_vector(cluster{i}) = i;
    end
elseif strcmp(method, 'rf')

    current_label=cluster2cv(X, temp_clus);

    %build random forest model

    %set leaf_size proportional to N to avoid overfitting, and for classification, typically min_leaf_size will not be greater than 10
    if leaf_size==0
        leaf_size=ceil(0.0025*N_labelled);
        if leaf_size>10
            leaf_size=10;
        end
    end
    t = templateTree('MinLeafSize', leaf_size);
    RF_model = fitcensemble(X(current_label~=0, :),current_label(current_label~=0),'Method','Bag','NumLearningCycles',ntrees,'Learners',t);
    pred=predict(RF_model,X(current_label==0, :));

    % update labels
    current_label(current_label == 0) = pred;
    cluster_vector = current_label;

    % update clusters
    cluster = cell(1, Kclus);
    for i = 1:Kclus
        obs_label=left_obs_label(pred==i);
        cluster{i} = [temp_clus{i}, obs_label];
    end
elseif strcmp(method, 'knn')
    current_label=cluster2cv(X, temp_clus);

    classes = current_label(current_label~=0);
    knnc = KNNdep(Kclus,K_knn,dm0(current_label==0,current_label~=0),classes);
    current_label(current_label==0) = knnc.class;
    cluster_vector = current_label;
    % update clusters
    cluster = cell(1, Kclus);
    for i = 1:Kclus
        obs_label=left_obs_label(knnc.class==i);
        cluster{i} = [temp_clus{i}, obs_label];
    end
end

if maxdepth
    nelabel = 0;
    while ~isempty(nelabel)
        if strcmp(depth, 'mahalanobis')
            DD = arrayfun(@(c) Maha_d(X, X(cluster{c}, :),'moment'), 1:Kclus, 'UniformOutput', false);
            DD = horzcat(DD{:});
        elseif strcmp(depth, 'spatial')
            DD = arrayfun(@(c) spatial_d(X, X(cluster{c}, :)), 1:Kclus, 'UniformOutput', false);
            DD = horzcat(DD{:});
        end

        [~, DDmat] = max(DD, [], 2);
        nelabel = find(DDmat ~= cluster_vector);
        cluster_vector = DDmat;
        if ~isempty(nelabel)
            for i = 1:length(cluster)
                cluster{i} = find(cluster_vector == i);
            end
        end
    end
end
result = struct('cluster', {cluster}, 'cluster_vector', cluster_vector);
end

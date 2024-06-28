function [output,depth_clus]= DLCC(data,dm0,dm0_order,ranktopobs,rank_mat,s,Th,method,class_method,varargin)
%main function of DLCC algorithm
%input:
  %data:data set
  %dm0: depth-based similarity matrix
  %dm0_order: the corresponding ordering matrix of dm0
  %ranktopobs: depth center of each subset
  %rank_mat: the depth rank of the obs in its own subset.
  %s: size of each subset
  %Th: threshold for grouping filtered centers
  %method: min or max strategy
  %class_method: the classification method for left obs, three options are
  %max depth classifier (maxdep), random forest (rf), and K-nearest neighbor (knn)
%optional input
  %maxdepth: logical parameter controls whether DLCC loops until all
  %observations in the allocated cluster have the highest depth values.(Default: false)
  %depth: type of depth (Default: spatial depth)
  %ntrees: only when the class_method is random forest, the number of trees (Default: 100)
  %leaf_size: only when the class_method is random forest, the minimum leaf size
  %K_knn:only when the class_method is KNN, # nearest neighbors (Default:5)
  %ifloop: logical parameter under the min strategy, if only remain one
  %center in each cluster, loop until the result doesn't change. (Default: false)
  %k: the number of cluster (if known), for the min strategy to give a
  %better intial grouping of centers.
  %initial_n: only for flexible min strategy, set a start number for
  %dropping centers in update_temp_cl.m
%Output:
  %temp_center: filtered centers in final clustering
  %temp_cluster: temporary clusters before classification and maxdep
  %depth_clus: final clustering results

% Create an input parser object
p = inputParser;

% Add required arguments
addRequired(p, 'data');
addRequired(p, 'dm0');
addRequired(p, 'dm0_order');
addRequired(p, 'ranktopobs');
addRequired(p, 'rank_mat');
addRequired(p, 's');
addRequired(p, 'Th');
addRequired(p, 'method');
addRequired(p, 'class_method');

% Define default values for optional arguments
defaultMaxdepth = false;
defaultDepth = 'spatial';
defaultLeafsize= 0;
defaultNtrees = 100;
defaultK_knn = 5;
defaultifloop=false;
defaultk=0;
defaultinitial_n=1;

% Add optional arguments to the parser
addParameter(p, 'maxdepth', defaultMaxdepth);
addParameter(p, 'depth', defaultDepth);
addParameter(p, 'ntrees', defaultNtrees);
addParameter(p, 'leaf_size', defaultLeafsize);
addParameter(p, 'K_knn', defaultK_knn);
addParameter(p, 'ifloop', defaultifloop);
addParameter(p,'k',defaultk)
addParameter(p,'initial_n',defaultinitial_n)


% Parse input arguments
parse(p, data,dm0,dm0_order,ranktopobs,rank_mat,s,Th,method,class_method, varargin{:});

% Extract parsed data from the parser object
data = p.Results.data;
dm0 = p.Results.dm0;
dm0_order=p.Results.dm0_order;
s=p.Results.s;
Th=p.Results.Th;
method = p.Results.method;
class_method=p.Results.class_method;
maxdepth = p.Results.maxdepth;
depth = p.Results.depth;
ntrees = p.Results.ntrees;
K_knn = p.Results.K_knn;
leaf_size=p.Results.leaf_size;
ifloop=p.Results.ifloop;
k=p.Results.k;
initial_n=p.Results.initial_n;
%begin the function

%filtering local centers
a=filter_center(dm0,dm0_order,s,ranktopobs,rank_mat,method);

if strcmp(method, 'min')
    %obtain groups of filtered centers
    if k~=0
        temp_cl=get_temp_cl_WK(a,s,dm0_order,k,initial_n);
    else
        temp_cl=get_temp_cl(a,dm0,s,Th,'min',dm0_order);
    end
    %update centers a
    a = cell2mat(temp_cl);

    if ifloop
        stop = 0;
        a_save = a;
        loop_i = 0;
        while stop ~= 1
            Kclus = length(temp_cl);
            temp_clus = cell(1,Kclus);
            for i = 1:Kclus
                [C,~,~] = unique(dm0_order(1:s,temp_cl{i}));
                temp_clus{i} = C;
            end
            temp_clus=get_temp_cluster(a,dm0,temp_clus,Kclus, temp_cl,s, 'min');
            depth_clus = left_class(data,temp_clus,class_method,'maxdepth',maxdepth,'depth',depth,'ntrees',ntrees,'leaf_size',leaf_size,'K_knn',K_knn,'dm0',dm0);
            Kclus = length(temp_clus);
            d=depth_by_cluster(data, Kclus, depth_clus.cluster, a_save, depth);
            [~,idx]=max(d);
            new_a = a_save(idx);

            % check any duplicate elements
            [new_a, ~, ic] = unique(new_a);
            counts = accumarray(ic, 1);
            has_duplicates = any(counts > 1);
            % if yes, randomly choose 1 other center instead of that one
            if has_duplicates
                duplicated_indices = counts > 1;
                duplicate = new_a(duplicated_indices);
                for v=1:length(duplicate)
                    diff_elements = setdiff(a_save, new_a);
                    new_a = [new_a, diff_elements(1)];
                end
            end

            if ~isequal(new_a, a)
                temp_cl = arrayfun(@(x) new_a(x), 1:Kclus, 'UniformOutput', false);
                a = new_a;
                loop_i = loop_i + 1;
                if loop_i > 10
                    warning('The result does not converge. Stopping to prevent possible infinite loops.');
                    stop = 1;
                end
            else
                %stop if new_a=a, it converges
                stop = 1;
            end
        end
    else
        Kclus = length(temp_cl);
        temp_clus = cell(1,Kclus);
        for i = 1:Kclus
            [C,~,~] = unique(dm0_order(1:s,temp_cl{i}));
            temp_clus{i} = C;
        end
        temp_clus=get_temp_cluster(a,dm0,temp_clus,Kclus, temp_cl,s, 'min');
        depth_clus = left_class(data,temp_clus,class_method,'maxdepth',maxdepth,'depth',depth,'ntrees',ntrees,'leaf_size',leaf_size,'K_knn',K_knn,'dm0',dm0);
    end
elseif strcmp(method, 'max')
    temp_cl=get_temp_cl(a,dm0,s,Th,'max',dm0_order);
    %update centers a
    a = cell2mat(temp_cl);

    Kclus = length(temp_cl);
    temp_clus = cell(1,Kclus);
    count_list = cell(1,Kclus);
    for i = 1:Kclus
        [C,~,groups] = unique(dm0_order(1:s,temp_cl{i}));
        counts = accumarray(groups(:),1);
        temp_clus{i} = C;
        count_list{i}=counts/length(temp_cl{i});
    end
    temp_clus=get_temp_cluster(a,dm0,temp_clus,Kclus, temp_cl,s, 'max',count_list);
    depth_clus = left_class(data,temp_clus,class_method,'maxdepth',maxdepth,'depth',depth,'ntrees',ntrees,'leaf_size',leaf_size,'K_knn',K_knn,'dm0',dm0);
end
output = struct('temp_center', {temp_cl}, 'temp_clus', {temp_clus});
end


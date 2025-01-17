function dm = Maha_dmcreator(data,cov_strategy,class,cov_mat_set)
%Built mahalanobis depth based similarity matrix
%Input: 
  %data: data set
  %cov_strategy: different methods in calculating covariance matrix
  %classes and cov_mat_set are used specifically for the "provided" strategy, indicating which points use which covariance matrix.

if nargin<3
    class=[];
    cov_mat_set=[];
end

n=size(data,1);
dm = zeros(n, n);

if strcmp(cov_strategy, 'moment')
    cov_matrix = cov(data);
    for i=1:n
        centered_x = data - data(i,:);
        md2 = sum(centered_x / cov_matrix .* centered_x, 2);
        dm(i,:) = 1 ./ (1 + md2);
    end
elseif strcmp(cov_strategy, 'MCD')
    cov_matrix = robustcov(data,'Method','olivehawkins');
    for i=1:n
        centered_x = data - data(i,:);
        md2 = sum(centered_x / cov_matrix .* centered_x, 2);
        dm(i,:) = 1 ./ (1 + md2);
    end
elseif strcmp(cov_strategy,'provided')
    for i=1:n-1
        centered_x = data - data(i,:);
        md2=sum(centered_x / cov_mat_set(:,:,class(i)) .* centered_x, 2);
        dm(i,:) = 1 ./ (1 + md2);
    end
elseif strcmp(cov_strategy,'Euclidean')
    D=size(data,2);
    cov_matrix=eye(D);
    for i=1:n
        centered_x = data - data(i,:);
        md2 = sum(centered_x / cov_matrix .* centered_x, 2);
        dm(i,:) = 1 ./ (1 + md2);
    end
% elseif strcmp(cov_strategy,'Reflec')
% for i=1:n
%     idx = ~(1:n==i);
%     ref_pot = 2*data(i,:) - data(idx,:);
%     data_ref=[data;ref_pot];
%     %cov_matrix = robustcov(data_ref,'Method','olivehawkins');
%     cov_matrix = cov(data_ref);
%     centered_x = data - data(i,:);
%     md2 = sum(centered_x * pinv(cov_matrix) .* centered_x, 2);
%     dm(i,:) = 1 ./ (1 + md2);
% end
end
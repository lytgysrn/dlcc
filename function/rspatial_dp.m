function [dm, Lmatrix_save] = rspatial_dp(data,gpu)
%build spatial-depth based similarity matrix

%data: input data for calculating reflection spatial depth wrt to each point
%gpu: if using gpu to accelerate the calculation or not, it is only
%recommended if the data set is large.

%output:dm-similarity matrix.
%Lmatrix: norm^2 matrix(norm^2 of any vector connecting two points in data)
if  nargin<2
    gpu = false;
end

tic
[n,d]= size(data);
rn_inv = 1/(2*n-1);
dm = zeros(n, n);

tol = 1e-5;


Ematrix=zeros(n,d);
Lmatrix=zeros(n,n);
for i = 1:n
    a = data(i,:) - data;
    Lmatrix(:,i) = vecnorm(a,2,2);
    norm_a = Lmatrix(:,i);
    norm_a(norm_a < tol) = 1;
    a = a./norm_a;
    Ematrix(i,:)=sum(a,1);
end
toc
%Ematrix actually is the sum of unit vector for each obs based on basic
%spatial depth
Lmatrix_save=Lmatrix;
Lmatrix= Lmatrix.^2;
%update Lmatrix to save the norm^2 between any two points in the data set

if gpu==true
    % Check if a GPU device is available
    if gpuDeviceCount == 0
        error('No GPU device found. Please ensure you have a supported GPU device and the Parallel Computing Toolbox installed.');
    end

    data = gpuArray(data);
    rn_inv = gpuArray(rn_inv);
    Lmatrix_gpu=gpuArray(Lmatrix);
    dm=gpuArray(dm);
    two_Lmatrix=2*Lmatrix_gpu;


    for j = 1:n
        idx = ~(1:n == j);
        %matrix of reflection points
        b_temp = -2 * data(j, :) + data(idx, :);
        %matrix of length from reflection points to each point in data
        norm_bM = sqrt(two_Lmatrix(idx, j) * ones(1, n, 'gpuArray') + two_Lmatrix(j, :) - Lmatrix_gpu(idx, :));

        index = (norm_bM < tol); % determine if any 0 values
        if sum(index, "all") > 0
            norm_bM(index) = inf;
        end
        norm_bM = 1 ./ norm_bM;
        C = (norm_bM.' *b_temp);
        norm_csum = sum(norm_bM, 1).';
        C2 = norm_csum .* data;
        C = C + C2 + Ematrix;
        dm(j, :) = 1 - vecnorm(C * rn_inv, 2, 2);
    end

    dm = gather(dm);

else

    two_Lmatrix = 2 * Lmatrix;


    for j=1:n
        idx = ~(1:n==j);
        b_temp = -2*data(j,:) + data(idx,:);
        norm_bM = sqrt(two_Lmatrix(idx, j) * ones(1, n) + two_Lmatrix(j, :) - Lmatrix(idx, :));

        index=(norm_bM < tol);%determine if any 0 values
        if sum(index,"all")>0
            norm_bM(index)=inf;
        end
        norm_bM=1./norm_bM;
        C=(b_temp.'*norm_bM).';
        norm_csum = sum(norm_bM,1).';
        C2=norm_csum.*data;
        C=C+C2+Ematrix;
        dm(j,:)=1-vecnorm(C*rn_inv,2,2);
    end
end



% tic
% Bmatrix=zeros(n,d);
% c=zeros(1,d);
% for j=1:n
% med=data(j,:);
% b_temp = -2*med + data(~(1:n==j),:);
% % norm_bM = sqrt(2 * Lmatrix(:, j) * ones(1, n) + 2 * Lmatrix(j, :) - Lmatrix);
% % norm_bM(j,:)=[];
%  for i=1:n
%   norm_b=sqrt(2*Lmatrix(:,j)+2*Lmatrix(j,i)-Lmatrix(:,i));
%   norm_b(j)=[];
%   idx_zero = (norm_b < tol); % find small values
%   if any(idx_zero)
%   norm_b(idx_zero)=1;
%   c=sum(1./norm_b(~idx_zero))*data(i,:);
%   b_temp= b_temp./norm_b;
%   g=sum(b_temp(idx_zero,:),1);
%   Bmatrix(i,:)=sum(b_temp,1)+c-g;
%  else
%    c=sum(1./norm_b)*data(i,:);
%    b_temp= b_temp./norm_b;
%    Bmatrix(i,:)=sum(b_temp,1)+c;
%   end
%     b_temp=b_temp.*norm_b;
%  end
%  Bmatrix=Bmatrix+Ematrix;
%  dm2(j,:)=1-vecnorm(Bmatrix*rn_inv,2,2);
% end
% toc


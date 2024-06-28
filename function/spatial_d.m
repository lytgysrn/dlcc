function dep = spatial_d(x, data, Lmatrix)
% spatial depth (Serfling 2000 version)
% x: whose depth is to be found
% data: sample points of this distribution
% Lmatrix: n*n norm matrix of data, if x=data, with Lmatrix, the depth can
% be computed more efficiently

tol = 1e-5;
[n,D]= size(data);

if nargin < 3
N = size(x,1);
E = zeros(N, D);    
for i = 1:N
    a = x(i,:)-data;
    norm_a = vecnorm(a,2,2);
    norm_a(norm_a<tol) = 1;
    a = a./norm_a;
    E(i,:)=sum(a,1);
end
dep=1-vecnorm(E/n,2,2);
else

index=(Lmatrix < tol);
Lmatrix(index)=inf;
Lmatrix=1./Lmatrix;
C=(data.'*Lmatrix).';
norm_csum = sum(Lmatrix,1).';
C2=norm_csum.*data;
C=C2-C;
dep=1-vecnorm(C/n,2,2);

end
  
end


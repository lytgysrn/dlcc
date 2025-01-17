function sim_mat = sim_mat(N,dm0_order, s, centers,centers2)
%similarity matrix between local centers in centers and centers2, if
%missing centers2 then center2=centers
if nargin<4
     sim_mat = zeros(N);
    for i = 1:(N-1)
        for j = (i+1):N
            sim = sum(ismember(dm0_order(:, i), dm0_order(:, j))) / s;
            sim_mat(i,j) = sim;
            sim_mat(j,i) = sim;
        end
    end
    sim_mat = sim_mat + diag(ones(N,1));
elseif nargin<5
    sim_mat = zeros(N);
    for i = 1:(N-1)
        idx=dm0_order(1:s,centers(i));
        for j = (i+1):N
            sim = sum(ismember(idx, dm0_order(1:s, centers(j)))) / s;
            sim_mat(i,j) = sim;
            sim_mat(j,i) = sim;
        end
    end
    sim_mat = sim_mat + diag(ones(N,1));
else
n=length(centers2);
sim_mat=zeros(N, n);
for i=1:N
    idx=dm0_order(1:s,centers(i));
    for j=1:n
        sim_mat(i,j)=sum(ismember(idx, dm0_order(1:s, centers2(j)))) / s;
    end
end
end
end
function best_model = EM_EEV(X, Kset, t)
%EM algorithm for EEV model in GPCM family (Gaussian parsimonious clustering models)

%Input:
 % X: data of observations
 % Kset: the number of components, a range is also acceptable.
 % t: maximum number of iterations
%Output:
 % The best BIC values under different K, and the information about the
 % best model (covariance matrix, mean vector, EEV model clustering results e.t.c)
 
[N, D] = size(X);
nK = length(Kset);
Qhis = cell(1, nK);
muhis = cell(1, nK);
covhis = cell(1, nK);
alphahis = cell(1, nK);
zlist = cell(1, nK);
iter_num=zeros(0,nK);
A_save=cell(1, nK);
for l = 1:nK
    K = Kset(l);
    alpha = repmat(1/K, 1, K);
    X_kmeans = kmeans(X, K, 'MaxIter', 100, 'Replicates', 100);
    mu = splitapply(@(x) mean(x), X, X_kmeans);
    cov_mat = zeros(D, D, K);
    for k = 1:K
        idx = (X_kmeans == k);
        cov_mat(:,:,k) = cov(X(idx, :));
    end

    Q = 0;

    for h = 2:(t + 1)
        [log_lik, z] = llh_calculator(N, K, alpha, X, mu, cov_mat);
        Q(h) = log_lik;

        eps = abs(Q(h) - Q(h - 1));

        if eps > 1e-6
            sumz = sum(z, 1);
            alpha = sumz / N;
            new_D = zeros(D, D, K);
            Omega = zeros(D, D, K);

            for j = 1:K
                zx = z(:, j).*X;
                mu(j, :) = sum(zx, 1) / sumz(j);
                x_min_mu = X-mu(j,:);
                z_x_min_mu = z(:, j).*x_min_mu;
                [V, V2] = eig((z_x_min_mu' * x_min_mu));
                [eigenvalues, sort_idx] = sort(diag(V2), 'descend');
                V = V(:, sort_idx);
                V2 = diag(eigenvalues);
                new_D(:,:,j) = V;
                Omega(:,:,j) = V2;
            end

            sumOmega = sum(Omega, 3);
            aveOmega = det(sumOmega) ^ (1 / D);
            A = sumOmega / aveOmega;
            lambda = aveOmega / N;

            for k = 1:K
                cov_mat(:,:,k) = lambda * new_D(:,:,k) * A * new_D(:,:,k)';
            end
        else
            break
        end
    end

    Qhis{l} = Q(end);
    muhis{l} = mu;
    covhis{l} = cov_mat;
    alphahis{l} = alpha;
    zlist{l} = z;
    A_save{l}=diag(A);
    iter_num(l)=h-1;
end

npara = Kset - 1 + Kset * D + D + Kset .* (D * (D - 1) / 2);
BIC_save = 2 * cell2mat(Qhis) - npara * log(N);

[~, bml] = max(BIC_save);
best_model.K = Kset(bml);
best_model.mu = muhis{bml};
best_model.cov = covhis{bml};
best_model.z = zlist{bml};
[~, classes] = max(best_model.z, [], 2);
best_model.classes=classes;
best_model.loglik = Qhis{bml};
best_model.BIC = BIC_save(bml);
best_model.BIC_all = BIC_save;
best_model.iter_num=iter_num;
best_model.shape=A_save{bml};
end
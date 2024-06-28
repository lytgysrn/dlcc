function [log_lik, z] = llh_calculator(N, K, alpha, data, mu, Sig)
%log likelihood calculator for EEV model in EM algorithm

D = size(data, 2);
p = zeros(N, K);

for g = 1:K
    p(:, g) = alpha(g) * exp(log_density(data, mu(g, :), Sig(:, :, g), D));
end

row_sum_term = sum(p, 2);
z = p./row_sum_term;
row_sum_term = log(row_sum_term);
log_lik = sum(row_sum_term);
end


function dist_sq = mahalanobis_distance(X, mu, C)
diff = X - mu;
temp = diff/C;
dist_sq = sum(temp .* diff, 2);
end

function ld = log_density(x, mu, Sig, D)
log_det = log(det(Sig));
cons_value = log(2) + log(pi);
cons_value = -(D * 0.5) * cons_value- (0.5 * log_det);
rep_mu=repmat(mu, size(x, 1), 1);
mh_distance = mahalanobis_distance(x, rep_mu, Sig);
mh_term = -0.5 * mh_distance;
ld = cons_value + mh_term;
end
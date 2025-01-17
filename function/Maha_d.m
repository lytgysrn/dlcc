function depth = Maha_d(x, data, cov_method)
%Compute Mahalanobis depth 

    % x: matrix containing the points whose depth is to be computed 
    % data: matrix containing the sample points of the distribution 
    % cov_method: method for estimating the covariance matrix
   
    mu = mean(data);
    
    try
        if strcmp(cov_method, 'moment')
            cov_matrix = cov(data);

        elseif strcmp(cov_method, 'MCD')
            %  Minimum Covariance Determinant
            cov_matrix = robustcov(data, 'Method', 'olivehawkins', 'OutlierFraction', 0.25);
            
            % singuarity check
            if rcond(cov_matrix) < 1e-10
                cov_matrix = eye(size(data, 2)); % default to identity matrix
            end

        else

            error('Invalid covariance estimation method. Choose either "moment" or "MCD".');
        end

    catch
        warning('Covariance estimation failed. Using identity matrix.');
        cov_matrix = eye(size(data, 2)); 
    end

    centered_x = x - mu;
    md2 = sum(centered_x / cov_matrix .* centered_x, 2);

    % Compute the Mahalanobis depth for each point in x
    depth = 1 ./ (1 + md2);
end


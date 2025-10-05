function [centroid_x, centroid_y, axis_major_length, axis_minor_length, axis_ratio_minor_to_major, orientation_radians, orientation_degrees] = fit_ellipse_pure_moments_cov(binary_image)
    assert(islogical(binary_image), 'fit_ellipse_pure_moments_cov:TypeError', 'input image must be type of Logical');
    
    % Initialize outputs
    centroid_x = NaN;
    centroid_y = NaN;
    axis_major_length = NaN;
    axis_minor_length = NaN;
    axis_ratio_minor_to_major = NaN;
    orientation_radians = NaN;
    orientation_degrees = NaN;
    
    [rows, cols] = find(binary_image);
    
    % Fit an ellipse using the 2nd moments of the image's convex hull
    if ~isempty(rows)
        % Calculate centroid
        x_center = mean(cols);
        y_center = mean(rows);
        mu20 = sum((cols - x_center).^2);
        mu02 = sum((rows - y_center).^2);
        mu11 = sum((cols - x_center) .* (rows - y_center));
        
        % Form the covariance matrix (normalized by number of points)
        n = length(cols);
        cov_matrix = [mu20/n, mu11/n; mu11/n, mu02/n];
        
        % Get eigenvalues and eigenvectors
        [eigenvecs, eigenvals] = eig(cov_matrix);
        eigenvals = diag(eigenvals);
        
        % Sort by eigenvalue magnitude
        [eigenvals, idx] = sort(eigenvals, 'descend');
        eigenvecs = eigenvecs(:, idx);
        
        % Ellipse parameters
        major_axis = 2 * sqrt(eigenvals(1)); % major axis
        minor_axis = 2 * sqrt(eigenvals(2)); % minor axis
        orientation_radians = atan2(eigenvecs(2,1), eigenvecs(1,1));
        orientation_degrees = rad2deg(orientation_radians);  % FIX: Was deg2rad (incorrect conversion)
        axis_ratio_minor_to_major = minor_axis / major_axis;
        axis_major_length = major_axis;
        axis_minor_length = minor_axis;
        centroid_x = x_center;
        centroid_y = y_center;
    end
end
%% Second Moment Ellipse Orientation Calculation
% This script demonstrates three methods to calculate the orientation of
% the second-moment ellipse for all non-black pixels in an image

% Clear workspace
clear; clc; close all;

%% Load or create sample image
% Option A: Load your own image
% img = imread('your_image.png');

% Option B: Create a sample binary image for demonstration
img = zeros(200, 200);
img(50:150, 60:140) = 1;  % Rectangle
img(80:120, 100:180) = 1; % Another rectangle
img = img + 0.3*randn(size(img)); % Add some noise
img(img < 0) = 0; % Ensure no negative values

% Display original image
figure;
subplot(2,2,1);
imshow(img, []);
title('Original Image');

%% Method 1: Single binary mask with regionprops
fprintf('=== Method 1: Single Binary Mask ===\n');

% Create a binary mask where all non-black pixels are white
if size(img, 3) == 3  % RGB image
    binary_mask1 = any(img > 0, 3);
else  % Grayscale image
    binary_mask1 = img > 0;
end

% Fill any holes to ensure it's treated as one region
binary_mask1 = imfill(binary_mask1, 'holes');

% Use regionprops to get orientation and other properties
stats1 = regionprops(binary_mask1, 'Orientation', 'MajorAxisLength', ...
    'MinorAxisLength', 'Centroid');

if ~isempty(stats1)
    orientation1 = stats1.Orientation;
    fprintf('Orientation: %.2f degrees\n', orientation1);
    fprintf('Major Axis Length: %.2f pixels\n', stats1.MajorAxisLength);
    fprintf('Minor Axis Length: %.2f pixels\n', stats1.MinorAxisLength);
else
    fprintf('No regions detected\n');
    orientation1 = NaN;
end

% Display binary mask
subplot(2,2,2);
imshow(binary_mask1);
title('Method 1: Binary Mask');

%% Method 2: Manual second-moment calculation
fprintf('\n=== Method 2: Manual Calculation ===\n');

% Get coordinates of all non-black pixels
if size(img, 3) == 3  % RGB image
    [rows, cols] = find(any(img > 0, 3));
else  % Grayscale image
    [rows, cols] = find(img > 0);
end

if ~isempty(rows)
    % Calculate centroid
    x_center = mean(cols);
    y_center = mean(rows);
    
    % Calculate second moments
    mu20 = sum((cols - x_center).^2);
    mu02 = sum((rows - y_center).^2);
    mu11 = sum((cols - x_center) .* (rows - y_center));
    
    % Calculate orientation (in degrees)
    orientation2 = 0.5 * atan2(2*mu11, mu20 - mu02) * 180/pi;
    
    fprintf('Orientation: %.2f degrees\n', orientation2);
    fprintf('Centroid: (%.2f, %.2f)\n', x_center, y_center);
    fprintf('Second moments - mu20: %.2f, mu02: %.2f, mu11: %.2f\n', ...
        mu20, mu02, mu11);
else
    fprintf('No non-black pixels found\n');
    orientation2 = NaN;
    x_center = NaN; y_center = NaN;
end

%% Method 3: Convex hull approach
fprintf('\n=== Method 3: Convex Hull ===\n');

% Create a convex hull of all non-black pixels
if size(img, 3) == 3  % RGB image
    binary_mask3 = any(img > 0, 3);
else  % Grayscale image
    binary_mask3 = img > 0;
end

[rows3, cols3] = find(binary_mask3);
if length(rows3) >= 3  % Need at least 3 points for convex hull
    try
        k = convhull(cols3, rows3);
        hull_mask = poly2mask(cols3(k), rows3(k), size(img,1), size(img,2));
        stats3 = regionprops(hull_mask, 'Orientation', 'MajorAxisLength', ...
            'MinorAxisLength');
        
        orientation3 = stats3.Orientation;
        fprintf('Orientation: %.2f degrees\n', orientation3);
        fprintf('Major Axis Length: %.2f pixels\n', stats3.MajorAxisLength);
        fprintf('Minor Axis Length: %.2f pixels\n', stats3.MinorAxisLength);
    catch
        fprintf('Error creating convex hull\n');
        orientation3 = NaN;
        hull_mask = false(size(img,1), size(img,2));
    end
else
    fprintf('Insufficient points for convex hull\n');
    orientation3 = NaN;
    hull_mask = false(size(img,1), size(img,2));
end

% Display convex hull
subplot(2,2,3);
imshow(hull_mask);
title('Method 3: Convex Hull');

%% Visualization of ellipse orientations
subplot(2,2,4);
imshow(img, []);
hold on;

% Plot ellipse for Method 1 if valid
if ~isempty(stats1) && ~isnan(orientation1)
    % Get ellipse parameters
    cx = stats1.Centroid(1);
    cy = stats1.Centroid(2);
    a = stats1.MajorAxisLength / 2;
    b = stats1.MinorAxisLength / 2;
    theta = deg2rad(orientation1);
    
    % Create ellipse points
    t = linspace(0, 2*pi, 100);
    x_ellipse = cx + a*cos(t)*cos(theta) - b*sin(t)*sin(theta);
    y_ellipse = cy + a*cos(t)*sin(theta) + b*sin(t)*cos(theta);
    
    plot(x_ellipse, y_ellipse, 'r-', 'LineWidth', 2);
    plot(cx, cy, 'r+', 'MarkerSize', 10, 'LineWidth', 2);
end

% Plot orientation line for Method 2 if valid
if ~isnan(orientation2)
    line_length = 50;
    theta2 = deg2rad(orientation2);
    x_line = x_center + line_length * cos(theta2) * [-1, 1];
    y_line = y_center + line_length * sin(theta2) * [-1, 1];
    plot(x_line, y_line, 'b--', 'LineWidth', 2);
    plot(x_center, y_center, 'bo', 'MarkerSize', 8);
end

title('Ellipse Overlays');
legend('Method 1 Ellipse', 'Method 1 Center', 'Method 2 Orientation', ...
    'Method 2 Center', 'Location', 'best');
hold off;

%% Summary
fprintf('\n=== Summary ===\n');
fprintf('Method 1 (regionprops):     %.2f degrees\n', orientation1);
fprintf('Method 2 (manual):          %.2f degrees\n', orientation2);
fprintf('Method 3 (convex hull):     %.2f degrees\n', orientation3);

% Note: Methods 1 and 2 should give very similar results
% Method 3 might differ as it uses only the convex hull boundary
function simple_region_props(rendition_directory, out_csv_base_path, file_ext)
    % Enter the directory of renditions to process here
    RENDITION_DIRECTORY = rendition_directory;
    % Base path for output CSVs (will append method name)
    OUT_PATH_BASE = out_csv_base_path;
    file_list = dir(fullfile(RENDITION_DIRECTORY, file_ext));
    
    % Filter out Apple metadata files and other hidden files
    valid_files = {};
    valid_names = {};
    for i = 1:length(file_list)
        current_file = file_list(i);
        if ~startsWith(current_file.name, ".") && ~startsWith(current_file.name, "._")
            valid_files{end+1} = fullfile(current_file.folder, current_file.name);
            valid_names{end+1} = current_file.name;
        end
    end
    
    % Convert to proper arrays
    num_images = length(valid_files);
    image_paths = valid_files'; % Convert to column cell array
    image_names = valid_names'; % Convert to column cell array
    fprintf('Found %d valid image files\n', num_images);
    
    if num_images == 0
        error('No valid image files found in directory');
    end
    
    %% Initialize Parallel Pool
    % Check if parallel pool is running, if not, start one
    if isempty(gcp('nocreate'))
        poolobj = parpool(); % Use default number of workers
        fprintf('Started parallel pool with %d workers\n', poolobj.NumWorkers);
    else
        poolobj = gcp;
        fprintf('Using existing parallel pool with %d workers\n', poolobj.NumWorkers);
    end
    
    % Initialize storage for all three methods
    all_props_convhull = cell(num_images, 1);
    all_props_filled = cell(num_images, 1);
    all_props_moments = cell(num_images, 1);
    all_names = cell(num_images, 1);
    
    % Setup progress tracking with DataQueue
    progress_queue = parallel.pool.DataQueue;
    progress_count = 0;
    start_time = tic;
    afterEach(progress_queue, @(~) updateProgress());
    
    % Nested function for progress updates
    function updateProgress()
        progress_count = progress_count + 1;
        if mod(progress_count, max(1, floor(num_images/20))) == 0 || progress_count == num_images
            elapsed_time = toc(start_time);
            avg_time_per_image = elapsed_time / progress_count;
            remaining_images = num_images - progress_count;
            estimated_remaining_time = remaining_images * avg_time_per_image;
            
            fprintf('Processed %d/%d images (%.1f%%) - Elapsed: %.1fs, ETA: %.1fs\n', ...
                    progress_count, num_images, (progress_count/num_images)*100, ...
                    elapsed_time, estimated_remaining_time);
        end
    end
    
    fprintf('Starting processing of %d images with 3 methods...\n', num_images);
    
    parfor i = 1:num_images
        try
            img = imread(image_paths{i});
            if size(img, 3) == 3
                img = rgb2gray(img);
            end
            threshold = 0;
            img_binary = img > threshold;
            
            
            % Method 1: Convex Hull
            [cx1, cy1, maj1, min1, ratio1, rad1, deg1] = fit_ellipse_conv_hull(img_binary);
            all_props_convhull{i} = struct('centroid_x', cx1, 'centroid_y', cy1, ...
                                           'axis_major_length', maj1, 'axis_minor_length', min1, ...
                                           'axis_ratio_minor_to_major', ratio1, ...
                                           'orientation_radians', rad1, 'orientation_degrees', deg1);
            
            % Method 2: Filled Region
            [cx2, cy2, maj2, min2, ratio2, rad2, deg2] = fit_ellipse_filled_region(img_binary);
            all_props_filled{i} = struct('centroid_x', cx2, 'centroid_y', cy2, ...
                                         'axis_major_length', maj2, 'axis_minor_length', min2, ...
                                         'axis_ratio_minor_to_major', ratio2, ...
                                         'orientation_radians', rad2, 'orientation_degrees', deg2);
            
            % Method 3: Pure Moments/Covariance
            [cx3, cy3, maj3, min3, ratio3, rad3, deg3] = fit_ellipse_pure_moments_cov(img_binary);
            all_props_moments{i} = struct('centroid_x', cx3, 'centroid_y', cy3, ...
                                          'axis_major_length', maj3, 'axis_minor_length', min3, ...
                                          'axis_ratio_minor_to_major', ratio3, ...
                                          'orientation_radians', rad3, 'orientation_degrees', deg3);
            
            all_names{i} = image_names{i};
            
            % Send progress update
            send(progress_queue, i);
            
        catch ME
            fprintf('Error reading image %d (%s): %s\n', i, image_paths{i}, ME.message);
            all_props_convhull{i} = [];
            all_props_filled{i} = [];
            all_props_moments{i} = [];
            all_names{i} = image_names{i};
            send(progress_queue, i); % Still count failed attempts
            continue;
        end
    end
    
    % Final processing summary
    total_time = toc(start_time);
    fprintf('Completed processing %d images in %.2f seconds (%.3f sec/image)\n', ...
            num_images, total_time, total_time/num_images);
    
    % Create three separate tables
    methods = {'convhull', 'filled', 'moments'};
    props_data = {all_props_convhull, all_props_filled, all_props_moments};
    
    for m = 1:length(methods)
        results_table = table();
        results_table.image = all_names;

        % Pre-allocate all numeric columns
        results_table.centroid_x = nan(num_images, 1);
        results_table.centroid_y = nan(num_images, 1);
        results_table.axis_major_length = nan(num_images, 1);
        results_table.axis_minor_length = nan(num_images, 1);
        results_table.axis_ratio_minor_to_major = nan(num_images, 1);
        results_table.orientation_radians = nan(num_images, 1);
        results_table.orientation_degrees = nan(num_images, 1);
        
        for i = 1:num_images
            if ~isempty(props_data{m}{i})
                props = props_data{m}{i};
                results_table.centroid_x(i) = props.centroid_x;
                results_table.centroid_y(i) = props.centroid_y;
                results_table.axis_major_length(i) = props.axis_major_length;
                results_table.axis_minor_length(i) = props.axis_minor_length;
                results_table.axis_ratio_minor_to_major(i) = props.axis_ratio_minor_to_major;
                results_table.orientation_radians(i) = props.orientation_radians;
                results_table.orientation_degrees(i) = props.orientation_degrees;
            else
                % Handle failed processing or no regions found
                results_table.centroid_x(i) = NaN;
                results_table.centroid_y(i) = NaN;
                results_table.axis_major_length(i) = NaN;
                results_table.axis_minor_length(i) = NaN;
                results_table.axis_ratio_minor_to_major(i) = NaN;
                results_table.orientation_radians(i) = NaN;
                results_table.orientation_degrees(i) = NaN;
            end
        end
        
        % Generate output path with method name
        [path, name, ext] = fileparts(OUT_PATH_BASE);
        out_path = fullfile(path, sprintf('%s_%s%s', name, methods{m}, ext));
        
        % Save to CSV
        writetable(results_table, out_path);
        fprintf('Results for %s method saved to %s\n', methods{m}, out_path);
    end
end
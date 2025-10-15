function [centroid_x, centroid_y, axis_major_length, axis_minor_length, axis_ratio_minor_to_major, orientation_radians, orientation_degrees] = fit_ellipse_conv_hull(binary_image)
    assert(islogical(binary_image), 'fit_ellipse_conv_hull:TypeError', 'input image must be type of Logical');
    
    % Initialize outputs
    centroid_x = NaN;
    centroid_y = NaN;
    axis_major_length = NaN;
    axis_minor_length = NaN;
    axis_ratio_minor_to_major = NaN;
    orientation_radians = NaN;
    orientation_degrees = NaN;
    
    [rows, cols] = find(binary_image);
    
    if length(rows) >= 3
        try
            hull = convhull(cols, rows);
            hull_mask = poly2mask(cols(hull), rows(hull), size(binary_image,1), size(binary_image,2));
            props = regionprops(hull_mask, 'Centroid', 'Orientation', 'MajorAxisLength', 'MinorAxisLength');
            
            if ~isempty(props)
                % FIX: Handle case where props is a structure array (multiple regions)
                if length(props) > 1
                    % Multiple regions found - use the largest one
                    areas = zeros(length(props), 1);
                    for k = 1:length(props)
                        areas(k) = props(k).MajorAxisLength * props(k).MinorAxisLength;
                    end
                    [~, idx] = max(areas);
                    props = props(idx);
                end
                
                % Now props is guaranteed to be a single structure
                centroid_x = props.Centroid(1);
                centroid_y = props.Centroid(2);
                axis_major_length = props.MajorAxisLength;
                axis_minor_length = props.MinorAxisLength;  % FIX: Was AxisMinorLength (typo)
                axis_ratio_minor_to_major = axis_minor_length / axis_major_length;
                orientation_radians = deg2rad(props.Orientation);
                orientation_degrees = props.Orientation;
            end
        catch ME
            fprintf('Error creating convex hull: %s\n', ME.message);
        end
    end
end
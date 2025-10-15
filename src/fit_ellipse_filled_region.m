function [centroid_x, centroid_y, axis_major_length, axis_minor_length, axis_ratio_minor_to_major, orientation_radians, orientation_degrees] = fit_ellipse_filled_region(binary_image)
    assert(islogical(binary_image), 'fit_ellipse_filled_region:TypeError', 'input image must be type of Logical');
    
    % Initialize outputs
    centroid_x = NaN;
    centroid_y = NaN;
    axis_major_length = NaN;
    axis_minor_length = NaN;
    axis_ratio_minor_to_major = NaN;
    orientation_radians = NaN;
    orientation_degrees = NaN;
    
    % Fill the holes
    binary_mask = imfill(binary_image, "holes");
    
    % Use regionprops like before
    props = regionprops(binary_mask, 'Orientation', 'MajorAxisLength', 'MinorAxisLength', 'Centroid');
    
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
end
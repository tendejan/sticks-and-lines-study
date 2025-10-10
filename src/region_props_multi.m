ROOT_DIR = "../data/null_datasets";

% ROOT_DIR = "../data/experimental_datasets";

% Get all items in root directory
all_items = dir(ROOT_DIR);

% Filter to get only subdirectories (excluding . and ..)
subdirs = all_items([all_items.isdir] & ~ismember({all_items.name}, {'.', '..'}));

for i = 1:length(subdirs)
    subdir_name = subdirs(i).name;
    subdir_path = fullfile(ROOT_DIR, subdir_name);

    renditions_dir = fullfile(subdir_path, 'renditions');
    csv_output_path = fullfile(subdir_path, sprintf("%s_region_props.csv", subdir_name));

    if isfolder(renditions_dir)
        try
            simple_reigon_props(renditions_dir, csv_output_path, "*.jpg");
            fprintf('Processed %s\n',subdir_name);
        catch ME
            fprintf ("Error processing %s: %s\n", subdir_name, ME.message);
        end
    else
        fprintf("No rendition directory found for %s\n", subdir_name);
    end
end
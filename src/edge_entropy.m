%% Function for Calculating Entropy of Edge Orientations (Parallelized Version)

%This work is a modified MATLAB port of 'Python code to calculate edge orientation
%entropy in digital images' by Anselm Brachmann used under CC BY 4.0 
%(https://creativecommons.org/licenses/by/4.0/). 

%The original python code can be found here: https://osf.io/bd8ma/, and is
%described in the following publication: Redies, C., Brachmann, A., &
%Wagemans, J. (2017). High entropy of edge orientations characterizes
%visual artworks from diverse cultural backgrounds. Vision Research, 133,
%130-144. http://doi.org/10.1016/j.visres.2017.02.00.

%This work requires the Image Processing Toolbox and Parallel Computing Toolbox.

% Modified by Tom Endejan to read a directory and dynamically process the
% images in parallel, preserving image name as primary key for csv

function edge_entropy(image_directory, first_order_csv_out, second_order_csv_out_dir)
    
    directory = image_directory;
    file_ext = '*.jpg';
    
    % Create output directories if they don't exist
    if ~exist(first_order_csv_out, 'dir')
        mkdir(first_order_csv_out);
    end
    
    if ~exist(second_order_csv_out_dir, 'dir')
        mkdir(second_order_csv_out_dir);
    end
    
    file_list = dir(fullfile(directory, file_ext));
    file_list = file_list(~[file_list.isdir]);
    num_images=length(file_list);
    image_paths = cell(num_images, 1);
    image_names = cell(num_images, 1);
    for i= 1:num_images
        image_paths{i} = fullfile(file_list(i).folder, file_list(i).name);
        image_names{i} = file_list(i).name;
    end
    
    %init locations to store first and second order entropies
    first_order_entropy= zeros(1, num_images);
    second_order_entropy = zeros(1, num_images);
    
    all_first_order_hists = cell(num_images, 1);
    all_second_order_matrices = cell(num_images, 1);
    
    %% Set Constants
    MAX_PIXELS = 300*400; %For resizing img
    MAX_DIAGONAL = 500;
    CIRC_BINS = 48; %Angular Direction Bins
    GABOR_BINS = 24; %Orientation Bins
    BINS_VEC = linspace(0, 2*pi, GABOR_BINS+1); %Orientations expressed in Radians
    BINS_VEC = BINS_VEC(1:end-1); %0 to 345 degrees
    FILTER_SIZE = 31; %31
    
    %for generating figures
    HIRES_FILTER_SIZE = 301;
    [X_hires, Y_hires] = meshgrid(linspace((-HIRES_FILTER_SIZE/2)+1, (HIRES_FILTER_SIZE/2), HIRES_FILTER_SIZE));
    
    % Scale coordinates to match original filter scale
    scale_factor = HIRES_FILTER_SIZE / FILTER_SIZE;
    X_hires = X_hires / scale_factor;
    Y_hires = Y_hires / scale_factor;
    
    %% Create Filter Bank
    
    %Set Parameters for Gabors
    phi = pi/2; %phase offset
    lambda=8; %wavelength
    frequency=1/lambda;
    hrsf = 4; %half-response spatial frequency bandwith
    sigma = (1/(pi*frequency)) * sqrt(log(2)/2) * (2^hrsf+1)/(2^hrsf-1);% See Kruizinga & Petkov, 1999, pg. 1396
    [X,Y] = meshgrid(linspace((-FILTER_SIZE/2)+1,(FILTER_SIZE/2),FILTER_SIZE));
    
    %Create Guassian
    gaussian = exp(-(X.*X + Y.*Y)/(2*sigma*sigma));
    
    % Create high-res Gaussian
    gaussian_hires = exp(-(X_hires.*X_hires + Y_hires.*Y_hires)/(2*sigma*sigma));
    
    %Create Filters
    filter_bank=zeros(FILTER_SIZE,FILTER_SIZE,GABOR_BINS);
    filter_bank_hires = zeros(HIRES_FILTER_SIZE, HIRES_FILTER_SIZE, GABOR_BINS);
    
    for i=1:length(BINS_VEC)
        sinusoid=2*pi*frequency*(sin(BINS_VEC(i))*(X)) + 2*pi*frequency*(cos(BINS_VEC(i))*(Y));
        sinusoid_hires = 2*pi*frequency*(sin(BINS_VEC(i))*(X_hires)) + 2*pi*frequency*(cos(BINS_VEC(i))*(Y_hires));
        gabor = gaussian.*cos(sinusoid+phi);
        gabor_hires = gaussian_hires .* cos(sinusoid_hires + phi);
        filter_bank(:,:,i)=gabor;
        filter_bank_hires(:,:,i) = gabor_hires;
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
    
    %% Main Loop for Edge Orientation Entropies (Parallelized)
    parfor i=1:length(image_paths)
        
        %Read Image
        img=imread(image_paths{i});
        
        %Determine Resize Factor to Stay Within Max Pixels (120,000)
        [h,w] = size(img(:,:,1));
        if h*w > MAX_PIXELS
            a=sqrt(MAX_PIXELS / (h*w));
            h_out=floor(a*h);
            w_out=floor(a*w);
            img=imresize(img,a,'lanczos3','OutputSize',[h_out w_out]);
        end
        
         %Ensure Image is Grayscale
        if size(img,3) > 1
            img=double(rgb2gray(img));
        else
            img=double(img);
        end
        
        %Filter Image
        filt_img = zeros(size(img,1), size(img,2), size(filter_bank,3));
        for j=1:size(filter_bank,3)
            filt_img(:,:,j)=imfilter(img,filter_bank(:,:,j),'conv','symmetric');
        end
        
        %Save Orientation and Response of Strongest Filter
        %resp_bin is Orientation , resp_val is Response
        [resp_val, resp_bin]=max(filt_img,[],3);
        
        %Retain Strongest 10,000 Edges
        sorted_resp_val=sort(resp_val(:),'descend');
        cutoff=sorted_resp_val(min(10000, length(sorted_resp_val)));
        resp_val(resp_val<cutoff)=0;
        
        %Determine Locations in Image of Strongest Edges
        [ex,ey] = find(resp_val');
        
        %Create Histogram for First-Order Entropy and Normalize
        first_order_hist = zeros(1,GABOR_BINS);
        for j = 1:GABOR_BINS
            first_order_hist(j) = sum(resp_val(resp_bin==j));
        end
        first_order_hist=first_order_hist / sum(first_order_hist);
        
        %Create 3-Dimensional Histogram for Second-Order Entropy
        
        %Obtain a look-up table for distances between edges
        edge_dims = size(resp_val);
        [xx, yy] = meshgrid(linspace(-edge_dims(2),edge_dims(2),2*edge_dims(2)+1),linspace(-edge_dims(1),edge_dims(1),2*edge_dims(1)+1));
        dist = sqrt(xx.^2+yy.^2);
        
        %Obtain list of orientations for each edge
        orientations = resp_bin(sub2ind(edge_dims,ey,ex));
        
        %Pre-Allocate Array for 3D Histogram (500 X 48 X 24 Bins)
        counts=zeros([MAX_DIAGONAL, CIRC_BINS, GABOR_BINS]);
        
        %Begin Pairwise Comparisons
        for j=1:length(ex) %j is the index for the reference edge
            
            %Determine Orientation Difference
            orientations_rel = orientations - orientations(j);
            orientations_rel = mod(orientations_rel + GABOR_BINS, GABOR_BINS);
            orientations_rel=orientations_rel+1; %add 1 to avoid indexing 0 in 3D Histogram
            
            %Determine Distance
            distance_col=(ex-ex(j)+edge_dims(2))+1;
            distance_row=(ey-ey(j)+edge_dims(1))+1;
            distance_rel=uint32(round(dist(sub2ind(size(dist),distance_row,distance_col))));
            distance_rel=distance_rel+1;
            distance_rel(distance_rel>=MAX_DIAGONAL) = MAX_DIAGONAL-1;
            
            %Determine Angular Direction
            direction = round((atan2((ey-1)-(ey(j)-1),(ex-1)-(ex(j)-1))/(2*pi)*CIRC_BINS + (orientations(j)-1)/double(GABOR_BINS)*CIRC_BINS));
            direction = mod(direction+CIRC_BINS,CIRC_BINS);
            direction=direction+1;
            
            %Calculate Product of Filter Responses
            product=resp_val(sub2ind(edge_dims,ey,ex)).*resp_val(sub2ind(edge_dims,ey(j),ex(j)));
            
            %Add Products to Appropriate Bin in Histogram
            for k = 1:length(product)
                counts(distance_rel(k),direction(k),orientations_rel(k))= counts(distance_rel(k),direction(k),orientations_rel(k)) + product(k);
            end
        end
        
        %Normalize 1D Histograms of Orientation Difference Bins
        normalized_counts=zeros([MAX_DIAGONAL, CIRC_BINS, GABOR_BINS]);
        counts_sum=sum(counts,3)+ 0.00001;
        for j= 1:GABOR_BINS
            normalized_counts(:,:,j) = counts(:,:,j) ./ counts_sum;
        end
        
        %Calculate First-Order Entropy
        first_order_entropy(i)=shannon_entropy(first_order_hist);
        
        %Obtain Matrix of Second-Order Entropies at each Distance and Direction
        [d,a,~]=size(normalized_counts);
        second_order_matrix=zeros(d,a);
        second_order_matrix_nan=zeros(d,a);
        for di = 1:d %Distance
            for ai = 1:a %Angular Direction
                second_order_matrix(di,ai)=shannon_entropy(reshape(normalized_counts(di,ai,:),[1 24])); %One-Dimensional Histogram of Edge Orientation Differences
                if counts_sum(di,ai) > 1
                    second_order_matrix_nan(di,ai)=second_order_matrix(di,ai);
                else
                    second_order_matrix_nan(di,ai)=NaN;
                end
            end
        end
        
        %Calculate Average Second-Order Entropy over 20-240px Distances and
        %all Directions. Index is 21:241 because distances begin a 0.
        second_order_entropy(i) = nanmean(second_order_matrix_nan(21:241,:),'all');
    
        %save the data for this image
        local_first_order_hist = first_order_hist;
        local_second_order_matrix = second_order_matrix;
        
        %store results in temporary variables to avoid broadcast variables
        all_first_order_hists{i} = local_first_order_hist;
        all_second_order_matrices{i} = local_second_order_matrix;
        
        fprintf('Processed image %d of %d: %s\n', i, num_images, image_names{i});
    end
    
    %% Save Results to CSV Files
    results_table = table(image_names, first_order_entropy', ...
        second_order_entropy', 'VariableNames', {'image', 'first_order_entropy', ...
        'second_order_entropy'});
    
    results_path = fullfile(first_order_csv_out, 'image_entropies.csv');
    writetable(results_table, results_path);
    fprintf('Results saved to %s\n', results_path);
    
    % Now save first order histograms
    first_order_hist_table = cell2table(all_first_order_hists, 'VariableNames', {'first_order_hist'});
    first_order_hist_table.image = image_names;
    
    % Convert cell arrays of histograms to individual columns
    hist_table = first_order_hist_table(:, {'image'});
    for j = 1:GABOR_BINS
        orientation_rad = BINS_VEC(j);
        orientation_deg = orientation_rad * 180/pi;
        bin_name = sprintf('%.1fdeg', orientation_deg);
        bin_values = cellfun(@(x) x(j), all_first_order_hists);
        hist_table.(bin_name) = bin_values;
    end
    
    hist_path = fullfile(first_order_csv_out, 'first_order_histograms.csv');
    writetable(hist_table, hist_path);
    fprintf('First order histograms saved to %s\n', hist_path);
    
    % Save second order matrices to the specified directory
    for i = 1:num_images
        % Create a filename based on image name (remove extension)
        [~, img_name, ~] = fileparts(image_names{i});
        matrix_path = fullfile(second_order_csv_out_dir, sprintf('%s_second_order_matrix.csv', img_name));
        
        % Convert matrix to table
        matrix_data = all_second_order_matrices{i};
        matrix_table = array2table(matrix_data);
        
        % Add column names (angle bins)
        col_names = cell(1, CIRC_BINS);
        for j = 1:CIRC_BINS
            angle_deg = (j-1)*(360/CIRC_BINS);
            col_names{j} = sprintf('%.1fdeg', angle_deg);
        end
        matrix_table.Properties.VariableNames = col_names;
        
        % Add row names (distance bins)
        row_names = cell(MAX_DIAGONAL, 1);
        for j = 1:MAX_DIAGONAL
            dist_pixels = j-1;
            row_names{j} = sprintf('dist_%dpx', dist_pixels);
        end
        matrix_table.Properties.RowNames = row_names;
        
        % Write to CSV
        writetable(matrix_table, matrix_path, 'WriteRowNames', true);
        fprintf('Second order matrix for %s saved to %s\n', image_names{i}, matrix_path);
    end
end

%% Define the shannon_entropy function
function out = shannon_entropy(a)
    if sum(a)~=1 && sum(a)>0
        a = a / sum(a); %Normalize Histogram 
    end
    v = a>0.0; %Determine Non-Zero Bins
    out=-nansum((a.*v) .* log2(a.*v));
end

%% Helper function for saving in parfor loops
function parsave(fname, data)
    save(fname, 'data');
end
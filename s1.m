%% Step 1: Organizing Data
% Define paths
data_dir = 'F:\dcm\DCM_collaboration_data';
output_dir = 'F:\dcm\test 16';
if ~exist(output_dir, 'dir')
   mkdir(output_dir);
end
% List of ROIs
rois = {'VTA', 'SN', 'NAc', 'OFC', 'PrL','aIns','CPu','hippo'};
num_rois = length(rois);
% Number of subjects
num_subjects = 6;
% Initialize a structure to store time series data
time_series_data = struct();
% Loop over each subject and each ROI to load and concatenate the time series
for s = 1:num_subjects
   subject_id = sprintf('ID_%02d', s);
   for r = 1:num_rois
       roi_name = rois{r};
      
       % Load pre-Naloxone (EPI_1) data
       pre_file = fullfile(data_dir, sprintf('avgMap_%s_%s_post_EPI_1.mat', roi_name, subject_id));
       pre_data = load(pre_file);
       pre_fieldname = fieldnames(pre_data);
       pre_time_series = pre_data.(pre_fieldname{1});
      
       % Load post-Naloxone (EPI_2) data
       post_file = fullfile(data_dir, sprintf('avgMap_%s_%s_post_EPI_2.mat', roi_name, subject_id));
       post_data = load(post_file);
       post_fieldname = fieldnames(post_data);
       post_time_series = post_data.(post_fieldname{1});
      
       % Concatenate pre and post time series
       time_series_data.(subject_id).(roi_name) = [pre_time_series; post_time_series];
   end
end
% Save the time series data structure for further use
save(fullfile(output_dir, 'time_series_data.mat'), 'time_series_data');

%% Step2: Model Specification
% This script prepares DCM models for each subject for the specified hypotheses
% and saves the models for estimation.
% Set up the environment
addpath('F:\spm12');
spm('defaults', 'FMRI');
spm_jobman('initcfg');
% Define directories
data_dir = 'F:\dcm\test 16'; % Directory with extracted time series
output_dir = 'F:\dcm\test 16\models'; % Directory to save DCM models
if ~exist(output_dir, 'dir')
   mkdir(output_dir);
end
% Define relevant ROIs and subjects
rois = {'VTA', 'SN', 'NAc', 'OFC', 'PrL','aIns','CPu','hippo'};
subjects = 1:6;
% TR value (Repetition Time)
TR = 2; % Set your TR value here
% Load the extracted time series data
load(fullfile(data_dir, 'time_series_data.mat'), 'time_series_data');
% Number of ROIs
num_rois = length(rois);
% Define specific connections: VTA/SN to NAc/OFC/PrL
connections = [1 3; 1 4; 1 5; 2 3; 2 4; 2 5]; % [row col] pairs of allowed connections
% Total number of possible combinations for the selected B matrix connections (2^6)
num_combinations = 2^length(connections);
% Prepare the DCM models
for s = subjects
  
       subject_id = sprintf('ID_%02d', s);
      
       % Create the DCM structure
       DCM = struct();
      
       % Define matrices based on the hypothesis
       A_matrix = ones(num_rois); % Fully connected model
       B_matrix = zeros(num_rois, num_rois, num_combinations); % 3D modulation matrix
       C_matrix = zeros(num_rois, 1); % Input matrix
      
       % Generate all possible combinations for the selected connections
       for comb = 0:num_combinations-1
           comb_binary = dec2bin(comb, length(connections)); % Get binary representation of the combination
           comb_matrix = zeros(num_rois, num_rois); % Start with an empty modulation matrix
          
           % Set the selected connections based on the binary combination
           for i = 1:length(connections)
               comb_matrix(connections(i, 1), connections(i, 2)) = str2double(comb_binary(i));
           end
          
           % Store this combination in the 3D B matrix
           B_matrix(:, :, comb + 1) = comb_matrix;
       end
       % Assign time series data
       DCM.Y.y = [];
       for r = 1:num_rois
           roi_time_series = time_series_data.(subject_id).(rois{r});
           DCM.Y.y = [DCM.Y.y roi_time_series];
       end
      
       % Assign additional variables for DCM.Y
       DCM.Y.name = rois; % Names of the ROIs
       DCM.Y.dt = TR; % Sampling interval (Repetition Time)
      
       % Assign matrices
       DCM.a = A_matrix; % Connectivity matrix
       DCM.b = B_matrix; % 3D modulation matrix with all possible combinations
       DCM.c = C_matrix; % Input matrix
      
       % Define the design matrix X
       num_timepoints = size(DCM.Y.y, 1);
       X = zeros(num_timepoints, 1); % Initialize the design matrix with zeros
      
       % Define naloxone injection timing
       naloxone_timepoint = 601; % Example time point; adjust if necessary
       if naloxone_timepoint < num_timepoints
           X(naloxone_timepoint:end) = 1; % Set values to 1 from naloxone injection onwards
       end
      
       % Assign design matrix X
       DCM.X = X; % Design matrix for external input
       % Define U as zeros since it should not influence the model
       DCM.U.u = ones(size(DCM.Y.y, 1), 1); % External inputs (no influence)
       DCM.U.name = {'No Influence'}; % Naming the input
       DCM.U.dt = TR; % External input sampling interval
       % Define DCM options with updated settings
       DCM.options.nonlinear = 0; % Linear model
       DCM.options.two_state = 0; % Single state model
       DCM.options.stochastic = 0; % Deterministic model
       DCM.options.centre = 0; % Use raw (non-centered) time series data
       DCM.options.induced = 0; % Deterministic model without stochastic fluctuations
       DCM.options.TR = TR; % Time interval between scans in seconds
      
       % Save the DCM model
       save(fullfile(output_dir, sprintf('DCM_Subject_%02d.mat', s)), 'DCM');
      
       disp(['DCM model prepared and saved for Subject ' num2str(s)]);
 
end
disp('DCM models preparation completed.');

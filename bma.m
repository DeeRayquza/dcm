%% BMA

% Set up the environment
addpath('F:\spm12');
spm('defaults', 'FMRI');
spm_jobman('initcfg');

% Define directories
results_dir = 'F:\dcm\test 16\results';
output_dir = 'F:\dcm\test 16\bma'; % Directory with BMR results
BMA_results = {}; % Initialize cell array for BMA results
Ep_B = {}; % Initialize cell array to hold Ep.B results

num_subjects = 3;
num_combinations = 64;

% Initialize a cell array to hold DCM structs for each subject and combination
DCM_estimated_all_subjects = cell(num_subjects, num_combinations); 

for s = 1:num_subjects
    subject_id = sprintf('%02d', s);
    for comb = 1:num_combinations
        % Load the estimated model for each combination
        model_file = fullfile(results_dir, sprintf('DCM_estimated_Subject_%s_Hyp2_Comb_%d.mat', subject_id, comb));
        if exist(model_file, 'file')
            load(model_file, 'DCM_estimated_comb');
            
            % Check if DCM_estimated_comb is a 1x1 cell array with a struct inside it
            if iscell(DCM_estimated_comb) && ~isempty(DCM_estimated_comb{1}) && isstruct(DCM_estimated_comb{1})
                % Directly store the DCM struct inside the appropriate cell
                DCM_estimated_all_subjects{s, comb} = DCM_estimated_comb{1};  
                disp(['Loading completed for subject ' num2str(s) ' combination ' num2str(comb)]);
            else
                warning('DCM structure is missing or invalid for Subject %s, Combination %d.', subject_id, comb);
            end
        else
            warning('Model file for Subject %s, Combination %d does not exist.', subject_id, comb);
        end
    end
end
%% 

% Perform Bayesian Model Averaging (BMA) for group analysis
valid_models = reshape(DCM_estimated_all_subjects(~cellfun('isempty', DCM_estimated_all_subjects)), [], 1);

if ~isempty(valid_models)
    BMA_group = spm_dcm_bma(valid_models); % Perform BMA on valid models
    disp('BMA completed successfully.');
else
    error('No valid models found for BMA.');
end


%BMA_group = spm_dcm_bma(DCM_estimated_all_subjects);

% Save all BMA results
save(fullfile(output_dir, 'BMA.mat'), 'BMA_group');
%% 

% Extract the B matrix (modulatory connections) from BMA
B_matrix = BMA_group.Ep.B;

% Define the regions of interest (ROIs) for the plot
rois = {'NAc', 'OFC', 'PrL'}; % Only the regions connected to VTA and SN

% Create data for VTA and SN connections to NAc, OFC, PrL
VTA_connections = B_matrix(1, [3, 4, 5]); % VTA to NAc, OFC, PrL
SN_connections = B_matrix(2, [3, 4, 5]);  % SN to NAc, OFC, PrL

% Combine VTA and SN connections into one array for plotting
connections = [VTA_connections; SN_connections];

% Create a new figure for the group plot
figure;

% Bar plot of VTA and SN connections
bar(connections'); % Transpose to plot correctly

% Set x-axis labels for ROIs (NAc, OFC, PrL)
set(gca, 'XTickLabel', rois, 'XTick', 1:length(rois));

% Add legend to distinguish between VTA and SN
legend({'VTA', 'SN'}, 'Location', 'northwest');

% Uniform Y-axis range for easy comparison
ylim([-0.5, 0.5]); % Adjust this range based on the data scale

% Label each bar with the corresponding ROI values
for i = 1:length(rois)
    text(i - 0.15, VTA_connections(i), sprintf('%.2f', VTA_connections(i)), 'VerticalAlignment', 'bottom');
    text(i + 0.15, SN_connections(i), sprintf('%.2f', SN_connections(i)), 'VerticalAlignment', 'bottom');
end

% Add labels and title
xlabel('Target ROIs');
ylabel('Modulatory Effect (B matrix)');
title('BMA: VTA and SN Modulatory Connections to NAc, OFC, PrL');

% Save the group plot if needed
saveas(gcf, fullfile(output_dir, 'BMA_B_Plot.png')); % Save the group plot

disp('BMA Plotting completed.');
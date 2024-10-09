%% Group Analysis for FFX After BMR

% Set up the environment
addpath('F:\spm12');
spm('defaults', 'FMRI');
spm_jobman('initcfg');

% Define directories
input_dir = 'F:\dcm\test 16\results';
output_dir = 'F:\dcm\test 16\bma'; % Directory with BMR results
all_selected_models = {}; % Initialize cell array to hold all selected models

% Loop through subjects to load BMR results and collect selected models
for s = 1:3 % Assuming 3 subjects
    % Load selected models for the subject
    load(fullfile(input_dir, sprintf('BMR_Selected_Subject_%02d.mat', s)), 'selected_models');
    
    % Append selected models for the subject to the overall list
    all_selected_models = [all_selected_models, selected_models]; % Combine models
end

% Perform Bayesian Model Averaging (BMA) for group analysis
BMA_group = spm_dcm_bma(all_selected_models); % Perform BMA on all models together

% Save group-level BMA results
save(fullfile(output_dir, 'BMA_Group.mat'), 'BMA_group');
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
title('Fixed Effects: VTA and SN Modulatory Connections to NAc, OFC, PrL');

% Save the group plot if needed
saveas(gcf, fullfile(output_dir, 'Ep_B_Group_Plot.png')); % Save the group plot

disp('Group analysis and plotting of Ep.B completed.');
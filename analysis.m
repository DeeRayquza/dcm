% Set up the environment
addpath('F:\spm12');
spm('defaults', 'FMRI');
spm_jobman('initcfg');

% Define directories
results_dir = 'F:\dcm\test 16\results';
output_dir_rfx = 'F:\dcm\test 16\rfx'; % Directory for RFX results
output_dir_ffx = 'F:\dcm\test 16\ffx';
output_dir_plots = 'F:\dcm\test 16\plots'; % Directory to save the plots

% Create the plot directory if it doesn't exist
if ~exist(output_dir_plots, 'dir')
    mkdir(output_dir_plots);
end

% Load estimated DCM models
num_subjects = 3; % Assuming 3 subjects
num_combinations = 64; % Assuming 64 combinations

% Initialize an array to hold free energy values
free_energy_array = zeros(num_subjects, num_combinations); 

for s = 1:num_subjects
    subject_id = sprintf('%02d', s);
    for comb = 1:num_combinations
        % Load the estimated model for each combination
        model_file = fullfile(results_dir, sprintf('DCM_estimated_Subject_%s_Hyp2_Comb_%d.mat', subject_id, comb));
        if exist(model_file, 'file')
            load(model_file, 'DCM_estimated_comb');

            % Store the DCM struct inside the appropriate cell
            DCM_estimated_all_subjects{s}{comb} = DCM_estimated_comb; 

            % Extract the free energy and store it in the array
            free_energy_array(s, comb) = DCM_estimated_comb{1}.F; 
            disp(['Loading completed for subject ' num2str(s) ' combination ' num2str(comb)]);
        else
            warning('Model file for Subject %s, Combination %d does not exist.', subject_id, comb);
        end
    end
end

%% FFX

% Perform Fixed Effects BMS using spm_BMS on all models
Nsamp = 1e6;  % Number of samples for FFX
do_plot = 1;  % Set to 1 to enable plotting of results
sampling = 0;  % Analytical computation
ecp = 1;  % Compute exceedance probabilities

% Call spm_BMS to perform FFX
[ffx_alpha, ffx_exp_r, ffx_xp, ffx_pxp, ffx_bor] = spm_BMS(free_energy_array, Nsamp, do_plot, sampling, ecp);

% Optionally save results to output directory
save(fullfile(output_dir_ffx, 'FFX_BMS_Results.mat'), 'ffx_alpha', 'ffx_exp_r', 'ffx_xp', 'ffx_pxp', 'ffx_bor');

%% RFX

% Perform Random Effects BMS using spm_BMS on all models
Nsamp = 1e6;  % Number of samples for RFX
do_plot = 1;  % Set to 1 to enable plotting of results
sampling = 1;  % Analytical computation
ecp = 1;  % Compute exceedance probabilities

% Call spm_BMS to perform RFX
[rfx_alpha, rfx_exp_r, rfx_xp, rfx_pxp, rfx_bor] = spm_BMS(free_energy_array, Nsamp, do_plot, sampling, ecp);

% Optionally save results to output directory
save(fullfile(output_dir_rfx, 'RFX_BMS_Results.mat'), 'rfx_alpha', 'rfx_exp_r', 'rfx_xp', 'rfx_pxp', 'rfx_bor');

%% Visualization of FFX and RFX results

% Define model labels
model_labels = arrayfun(@(k) sprintf('Model %d', k), 1:num_combinations, 'UniformOutput', false);

% Generate and save plots for FFX and RFX
plot_ffx_rfx(ffx_alpha, rfx_alpha, rfx_xp, rfx_pxp, rfx_bor, model_labels, output_dir_plots);

% Define the function for plotting and saving the figures
function plot_ffx_rfx(alpha_ffx, alpha_rfx, xp_rfx, pxp_rfx, bor, model_labels, output_dir_plots)
    % Inputs:
    % alpha_ffx: Posterior model probabilities for FFX
    % alpha_rfx: Posterior model frequencies for RFX
    % xp_rfx: Exceedance probabilities for RFX
    % pxp_rfx: Protected exceedance probabilities for RFX
    % bor: Bayes Omnibus Risk (RFX)
    % model_labels: Cell array of strings with model labels (e.g., {'Model 1', 'Model 2', ...})
    % output_dir_plots: Directory to save plots
    
    Nk = length(alpha_ffx);  % Number of models
    if nargin < 6
        model_labels = arrayfun(@(k) sprintf('Model %d', k), 1:Nk, 'UniformOutput', false);
    end
    
    %% Plot for FFX: Posterior Model Probabilities
    figure;
    subplot(2, 1, 1);
    bar(alpha_ffx, 'FaceColor', [0.2 0.6 0.8]);
    set(gca, 'XTick', 1:Nk, 'XTickLabel', model_labels, 'XTickLabelRotation', 45);
    ylabel('Posterior Probability (FFX)');
    title('Fixed Effects (FFX) Model Selection');
    ylim([0 3]);
    
    % Save the FFX plot
    saveas(gcf, fullfile(output_dir_plots, 'FFX_Posterior_Probabilities.png'));

    %% Plot for RFX: Posterior Model Frequencies
    subplot(2, 1, 2);
    bar(alpha_rfx, 'FaceColor', [0.2 0.6 0.2]);
    set(gca, 'XTick', 1:Nk, 'XTickLabel', model_labels, 'XTickLabelRotation', 45);
    ylabel('Posterior Frequency (RFX)');
    title('Random Effects (RFX) Model Selection');
    ylim([0 3]);
    
    % Save the RFX plot
    saveas(gcf, fullfile(output_dir_plots, 'RFX_Posterior_Frequencies.png'));
    
    %% Plot for Exceedance Probabilities (RFX)
    figure;
    subplot(2, 1, 1);
    bar(xp_rfx, 'FaceColor', [0.8 0.4 0.2]);
    set(gca, 'XTick', 1:Nk, 'XTickLabel', model_labels, 'XTickLabelRotation', 45);
    ylabel('Exceedance Probability');
    title('Exceedance Probabilities (RFX)');
    ylim([0 0.1]);
    
    % Save the Exceedance Probabilities plot
    saveas(gcf, fullfile(output_dir_plots, 'RFX_Exceedance_Probabilities.png'));

    %% Plot for Protected Exceedance Probabilities (RFX)
    subplot(2, 1, 2);
    bar(pxp_rfx, 'FaceColor', [0.6 0.2 0.6]);
    set(gca, 'XTick', 1:Nk, 'XTickLabel', model_labels, 'XTickLabelRotation', 45);
    ylabel('Protected Exceedance Probability');
    title(sprintf('Protected Exceedance Probabilities (BOR = %.3f)', bor));
    ylim([0 0.1]);

    % Save the Protected Exceedance Probabilities plot
    saveas(gcf, fullfile(output_dir_plots, 'RFX_Protected_Exceedance_Probabilities.png'));
end

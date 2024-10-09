%% RFX Analysis on Variable Number of Models

% Set up the environment
addpath('F:\spm12');
spm('defaults', 'FMRI');
spm_jobman('initcfg');

% Define directories
results_dir = 'F:\dcm\test 16\results';
output_dir = 'F:\dcm\test 16\rfx'; % Directory for RFX results

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
%% 

% Perform Random Effects BMS using spm_BMS on all models
Nsamp = 1e6;  % Number of samples for RFX
do_plot = 1;  % Set to 1 to enable plotting of results
sampling = 1;  % Analytical computation
ecp = 1;  % Compute exceedance probabilities

% Call spm_BMS to perform RFX
[alpha, exp_r, xp, pxp, bor] = spm_BMS(free_energy_array, Nsamp, do_plot, sampling, ecp);

% Display the results
disp('Posterior probabilities of models (alpha):');
disp(alpha);
disp('Expected posterior probabilities (exp_r):');
disp(exp_r);
disp('Exceedance probabilities (xp):');
disp(xp);
disp('Protected exceedance probabilities (pxp):');
disp(pxp);
disp('Bayes Omnibus Risk (bor):');
disp(bor);

% Select the best model based on the highest posterior probability
[~, best_model_idx] = max(alpha);  % Find index of best model
disp('Best model index:');
disp(best_model_idx);

% Optionally save results to output directory
save(fullfile(output_dir, 'RFX_BMS_Results.mat'), 'alpha', 'exp_r', 'xp', 'pxp', 'bor', 'best_model_idx');

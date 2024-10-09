%% BMR for DCM Models

% Set up the environment
addpath('F:\spm12');
spm('defaults', 'FMRI');
spm_jobman('initcfg');

% Define directories
results_dir = 'F:\dcm\test 16\results'; % Directory with estimated DCM results
output_dir = 'F:\dcm\test 16\bmr';
subjects = 1:3; % Define subjects

% Initialize cell arrays to store models for all subjects
DCM_estimated_all_subjects = cell(length(subjects), 1);

% Load estimated DCM models
for s = subjects
    subject_id = sprintf('%02d', s);
    for comb = 1:64 % Assuming 64 combinations
        % Load the estimated model for each combination
        model_file = fullfile(results_dir, sprintf('DCM_estimated_Subject_%s_Hyp2_Comb_%d.mat', subject_id, comb));
        if exist(model_file, 'file')
            load(model_file, 'DCM_estimated_comb');
            % Store the DCM struct inside the appropriate cell
            DCM_estimated_all_subjects{s}{comb} = DCM_estimated_comb; 
            disp(['Loading completed for subject ' num2str(s) ' combination ' num2str(comb)]);
        else
            warning('Model file for Subject %s, Combination %d does not exist.', subject_id, comb);
        end
    end
end

% Perform Bayesian Model Reduction (BMR)
for s = 1:length(subjects)
    % Prepare cell array for models of the current subject
    models_for_subject = cell(1, 64); % Assuming 64 combinations
    for comb = 1:64
        models_for_subject{comb} = DCM_estimated_all_subjects{s}{comb}{1}; % Accessing the struct
    end

    % Perform BMR on the models
    [reduced_models, BMC, BMA] = spm_dcm_bmr(models_for_subject);

    % Eliminate models based on posterior probabilities
    threshold = 0.01; % Define a threshold for posterior probabilities
    selected_models = reduced_models([BMC.P] > threshold); % Keep only models with probabilities above the threshold

    % Save selected models
    save(fullfile(output_dir, sprintf('BMR_Selected_Subject_%02d.mat', subjects(s))), 'selected_models', 'BMC', 'BMA');
end

disp('BMR completed and models eliminated for all subjects.');

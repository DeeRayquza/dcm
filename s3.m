%% Step3: Estimate using parallel processing (with save outside of parfor loop)
% Set up the environment
addpath('F:\spm12');
spm('defaults', 'FMRI');
spm_jobman('initcfg');
% Define directories
model_dir = 'F:\dcm\test 16\models'; % Directory with DCM models
results_dir = 'F:\dcm\test 16\results'; % Directory to save DCM results
if ~exist(results_dir, 'dir')
   mkdir(results_dir);
end
% Define subjects
subjects = 4:6;
% Define the number of combinations from the B matrix (assuming 64)
num_combinations = 64;
% Initialize parallel pool if not already started
if isempty(gcp('nocreate'))
   parpool; % Starts a parallel pool of workers if not already started
end
% Estimate the DCM models for Hypothesis 2 for each subject
for s = subjects
   subject_id = sprintf('%02d', s);
  
   % Load the DCM model for the subject
   model_file_H2 = fullfile(model_dir, sprintf('DCM_Subject_%s.mat', subject_id));
   if exist(model_file_H2, 'file')
       load(model_file_H2, 'DCM');
      
       % Initialize a cell array to store the estimated DCMs
       DCM_estimated_all_combinations = cell(1, num_combinations);
      
       % Use parfor for fitting all combinations in parallel
       parfor comb = 1:num_combinations
           DCM_temp = DCM; % Create a copy of the DCM structure
           DCM_temp.b = DCM.b(:, :, comb); % Select the current combination
          
           try
               % Estimate the DCM model
               DCM_estimated_H2 = spm_dcm_fit(DCM_temp);
               DCM_estimated_all_combinations{comb} = DCM_estimated_H2; % Store the result in the cell array
              
               disp(['DCM model for Subject ' num2str(s) ', Combination ' num2str(comb) ' estimated.']);
           catch ME
               disp(['Error during DCM fitting for Subject ' num2str(s) ', Combination ' num2str(comb) ': ', ME.message]);
           end
       end
      
       % Save the estimated models outside of the parfor loop
       for comb = 1:num_combinations
           if ~isempty(DCM_estimated_all_combinations{comb})
               % Assign the cell element to a temporary variable
               DCM_estimated_comb = DCM_estimated_all_combinations{comb};
              
               % Save the estimated model
               save(fullfile(results_dir, sprintf('DCM_estimated_Subject_%s_Hyp2_Comb_%d.mat',subject_id, comb)), 'DCM_estimated_comb');
           end
       end
   else
       warning('Model file for Hypothesis 2, Subject %s does not exist.', subject_id);
   end
end
disp('DCM estimation completed for all subjects.');


% try to unpack some of the behavioural data in the DARPAK study.
% NOTE: 'alpha' seems to refer to 'contrast', probably because the alpha
% channel was used to modify contrast during each trial.

%load the data:
%load ('C:\Users\RyanDataPC\Documents\Stefan Honours theses\DARPAK\data\behavioural\DARPA_EEG_Data_K4.mat')

%%
% Define location of bh data:
base_data_dir = 'C:\Users\RyanDataPC\Documents\Stefan Honours theses\DARPAK\data\behavioural';
base_file_name = 'DARPA_EEG_Data_K'; % First part of the data file names

n = 22;
SubjNos = 4:25; % only subjects K4-K25 are analysed in the final regression: n=22

%%
% Set up containers to hold the data, sets of cell arrays filled with NaNs
% Remember we have subject (22) * block (7) per condition * trial (up to 25)
entropy.instructive = cell(1,n);
entropy.instructive(:) = {nan(7,25)};
entropy.monetary = cell(1,n);
entropy.monetary(:) = {nan(7,25)};

MI.instructive = cell(1,n);
MI.instructive(:) = {nan(7,25)};
MI.monetary = cell(1,n);
MI.monetary(:) = {nan(7,25)};

% container for behavioural responses:
% these are in absolute contrast difference in the chosen from the target contrast.
behav.instructive = cell(1,n);
behav.instructive(:) = {nan(7,25)};
behav.monetary = cell(1,n);
behav.monetary(:) = {nan(7,25)};

% Set up container for reward values.
% for the monetary condition, these are in dollers and are no greater than
% 0.25. For the instructive condition, these are in %contrast and are
% no greater than 15. Of course, no actual reward is given in the instructive condition,
% just information.
reward.instructive = cell(1,n);
reward.instructive(:) = {nan(7,25)};
reward.monetary = cell(1,n);
reward.monetary(:) = {nan(7,25)};

% we want to average across blocks for each subject
% and then across subjects
% so we have an indication of mean values across trials

%%
% Go through each SUBJECT, load and process their data in turn.
for S = 1:n
    
    % Load RAW bh data file for current subject:
    load (fullfile(base_data_dir, [base_file_name, num2str(SubjNos(S)), '.mat']))
    % first, work out whether first set of 7 blocks was instructive (0) or monetary (1):
    reward_block_first(S) = params.rewardFirst;
    
    % Now load the modelling results for current subject.[
    % Load block labels:
    load (fullfile(base_data_dir, 'regression_labels', ['K', num2str(SubjNos(S))], 'BlockLabels.mat'))
    % Load trial labels:
    load (fullfile(base_data_dir, 'regression_labels', ['K', num2str(SubjNos(S))], 'TrialLabels.mat'))
    % Load entropy values:
    load (fullfile(base_data_dir, 'regression_labels', ['K', num2str(SubjNos(S))], 'EntropyLabels.mat'))
    
    % Load mutual information labels:
    % NOTE: most subjects have 'NewInfoLabels', others just 'InfoLabels'.
    % We'll go with the new ones if they are present.
    if exist(fullfile(base_data_dir, 'regression_labels', ['K', num2str(SubjNos(S))], ...
            'NewInfoLabels.mat'), 'file')
        load (fullfile(base_data_dir, 'regression_labels', ['K', num2str(SubjNos(S))], ...
            'NewInfoLabels.mat'))
    else
        load (fullfile(base_data_dir, 'regression_labels', ['K', num2str(SubjNos(S))], 'InfoLabels.mat'))
    end
    
    % Now we should have the block and trial labels loaded,
    % as well as the entropy and MI.
    
    % Pull out the monetary values based on the trial and block labels.
    % Also pull out the behavioural data: absolute differences between chosen and target contrast
    % Looping across blocks and trials.
    % NOTE: rarely, sometimes 2 reward values are given for a single trial.
    % This must be a bug of some sort. To deal with this, we just take the first of
    % these values ie trialLog(blk,trial).reward(1);
    for blk = 1:7 % Blocks 1 to 7 first:
        for trial = 1:25
            idx =  find (urBlockLabels==blk & urTrialLabels == trial);
            if ~isempty(idx)
                % If blocks 1 to 7 are the monetary blocks:
                if reward_block_first(S)
                    entropy.monetary{S}(blk,trial) = urEntropyLabels(idx);
                    MI.monetary{S}(blk,trial) = urInfoLabels(idx);
                    behav.monetary{S}(blk,trial) = abs(event.rewardAlpha(blk) - trialLog(blk,trial).alpha);
                    reward.monetary{S}(blk,trial) = trialLog(blk,trial).reward(1);
                    % Otherwise if block 1 to 7 are the instructive blocks:
                else
                    entropy.instructive{S}(blk,trial) = urEntropyLabels(idx);
                    MI.instructive{S}(blk,trial) = urInfoLabels(idx);
                    behav.instructive{S}(blk,trial) = abs(event.rewardAlpha(blk) - trialLog(blk,trial).alpha);
                    reward.instructive{S}(blk,trial) = trialLog(blk,trial).reward(1);
                end
            end
        end
    end
    
    % Now blocks 8 to 14.
    % NOTE: there are 7 blocks here: numbered 8 to 14.
    % Because we are separating data out according to condition type, we therefore loop 1:7 but add 7 when indexing the blocks,
    % but not when placing the data in the 1:7 containers for each block per condition.
    for blk = 1:7
        for trial = 1:25
            idx =  find (urBlockLabels==blk+7 & urTrialLabels == trial); % add 7 to the block making them 8:14
            if ~isempty(idx)
                %If the reward block is first, these latter blocks will be the instructive ones
                if reward_block_first(S)
                    entropy.instructive{S}(blk,trial) = urEntropyLabels(idx);
                    MI.instructive{S}(blk,trial) = urInfoLabels(idx);
                    behav.instructive{S}(blk,trial) = abs(event.rewardAlpha(blk+7) - trialLog(blk+7,trial).alpha);
                    reward.instructive{S}(blk,trial) = trialLog(blk+7,trial).reward(1);
                    %Otherwise these latter blocks will be for monetary
                else
                    entropy.monetary{S}(blk,trial) = urEntropyLabels(idx);
                    MI.monetary{S}(blk,trial) = urInfoLabels(idx);
                    behav.monetary{S}(blk,trial) = abs(event.rewardAlpha(blk+7) - trialLog(blk+7,trial).alpha);
                    reward.monetary{S}(blk,trial) = trialLog(blk+7,trial).reward(1);
                end
            end
        end     % End of loop across TRIALS
    end         % End of loop across BLOCKS
end             % End of loop across SUBJECTS.

% We now should have the modelling results pulled out separately for MI and entropy,
% with a separate cell for each subject, 7 rows per condition with (max of) 25 trials.
% as well as the overall behavioural performance, and the reward values (or info) given.

%%
% Compute mean values across BLOCKS for each subject
entropy.mean_instructive = cellfun(@nanmean, entropy.instructive, 'UniformOutput', false);
entropy.mean_instructive = cat(1, entropy.mean_instructive{:}); %Concatenate across subjects into single matrix
entropy.mean_monetary = cellfun(@nanmean, entropy.monetary, 'UniformOutput', false);
entropy.mean_monetary = cat(1, entropy.mean_monetary{:});

MI.mean_instructive = cellfun(@nanmean, MI.instructive, 'UniformOutput', false);
MI.mean_instructive = cat(1, MI.mean_instructive{:}); %Concatenate across subjects into single matrix
MI.mean_monetary = cellfun(@nanmean, MI.monetary, 'UniformOutput', false);
MI.mean_monetary = cat(1, MI.mean_monetary{:});

behav.mean_instructive = cellfun(@nanmean, behav.instructive, 'UniformOutput', false);
behav.mean_instructive = cat(1, behav.mean_instructive{:}); %Concatenate across subjects into single matrix
behav.mean_monetary = cellfun(@nanmean, behav.monetary, 'UniformOutput', false);
behav.mean_monetary = cat(1, behav.mean_monetary{:});

reward.mean_instructive = cellfun(@nanmean, reward.instructive, 'UniformOutput', false);
reward.mean_instructive = cat(1, reward.mean_instructive{:}); %Concatenate across subjects into single matrix
reward.mean_monetary = cellfun(@nanmean, reward.monetary, 'UniformOutput', false);
reward.mean_monetary = cat(1, reward.mean_monetary{:});

% Save these structs of data as a mat file
save([base_data_dir, '\DARKPAK_proc_model_results.mat'], 'MI', 'entropy', 'behav', 'reward'); %, '-v7.3')

%% The below can be used to extract the behavioural performance data on the task.
%Find out how many trials in each block
% for ii = 1:length(event.rewardHistory)
%     x(ii) = size(event.rewardHistory{ii},2);
% end
%
% % Pull out chosen contrasts and compute absolute difference between
% % them and the target contrasts.
% All_alphas = nan(size(trialLog));
% Abs_contr_diff = nan(size(trialLog));
% for ii = 1:size(trialLog,1)                %loop across blocks
%     for jj = 1:size(trialLog,2)            %loop across trials
%         if ~isempty(trialLog(ii,jj).alpha) %If there's data for this trial
%             All_alphas(ii,jj) = trialLog(ii,jj).alpha;
%             % Work out abs difference between the TARGET and CHOSEN contrast:
%             Abs_contr_diff(ii,jj) = abs(event.rewardAlpha(ii) - trialLog(ii,jj).alpha) ;
%         end
%     end
% end




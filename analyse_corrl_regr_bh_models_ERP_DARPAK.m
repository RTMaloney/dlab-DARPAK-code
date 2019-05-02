% Perform various analyses on modelling, behaviour and ERP results for the DARPAK study.
% Analyses mostly follow those performed in Bennett et al 2015
% R Maloney, 22 Dec 2018

%%
n = 22; %we have modelling results for 22 subjects.

% These have been processed from matlab data in the try_to_unpack_bh_data_DARPAK.m script, in the form of a .mat file.
% Load the mat file.
% Note this also includes the behavioural performance and reward values.
model_results = load ('C:\Users\RyanDataPC\Documents\Stefan Honours theses\DARPAK\data\behavioural\DARKPAK_proc_model_results.mat');

% Load perceptual uncertainty values (sigma), instructive
instr_model = load('C:\Users\RyanDataPC\Documents\Stefan Honours theses\DARPAK\data\model\FitParameters_Instructive.mat');
monet_model = load('C:\Users\RyanDataPC\Documents\Stefan Honours theses\DARPAK\data\model\FitParameters_Reward.mat');

%% Perform correlation: behaviour and perceptual uncertainty (sigma)
% average across TRIALS for each subject.
% Reward (or feedback instruction) provided on each trial:
mean_bh_instr = nanmean(model_results.behav.mean_instructive, 2); % average in % contrast
mean_bh_monet = nanmean(model_results.behav.mean_monetary, 2);    % average in CENTS

% Perform the correlation, remembering to omit the first subject (K3), who we are missing
% data for (modelling, bh, ERPs)
% Monetary feedback:
[R,P] = corrcoef(mean_bh_monet, monet_model.x(2:end))
% Instructive feedback:
[R,P] = corrcoef(mean_bh_instr, instr_model.x(2:end))

% t-test comparing sigmas:
%[h, p] = ttest(monet_model.x(2:end), instr_model.x(2:end))

%% Now do regressions: behaviour and modelling

% First, effects of belief uncertainty (entropy), and trial number
% on behaviour. These will be done across subjects, on a trial-by-trial basis.
% We then report the mean beta, as in Bennett et al 2015.

stats_instr = [];
betas_instr = [];
stats_monet = [];
betas_monet = [];
for S = 1:n % Loop across subjects.
    
    % *** Do instructive feedback first *** :
    idx = ~isnan(model_results.behav.mean_instructive(S,:)); % grab trials where data exists
    y = model_results.behav.mean_instructive(S,idx)';
    
    x1 = ones(length(y),1); % need a column of 1s for a regression (the constant term)
    trial_no = [1:length(y)]';
    ent = model_results.entropy.mean_instructive(S,idx)'; % pull out belief uncertainty: entropy
    
    % Perform the regression.
    % the vector stats contains the R2 statistic, the F-statistic and its p-value, and an estimate of the error variance.
    [b,bint,~,~,stats_instr(S,:)] = regress(y,[x1, trial_no, ent]);
    
    betas_instr(S,:) = b'; % set aside betas: each row is a subject, each column for the 3 predictors
    
    % *** Now do monetary feedback in the same way *** :
    idx = ~isnan(model_results.behav.mean_monetary(S,:)); % grab trials where data exists
    y = model_results.behav.mean_monetary(S,idx)';
    
    x1 = ones(length(y),1); % need a column of 1s for a regression (the constant term)
    trial_no = [1:length(y)]';
    ent = model_results.entropy.mean_monetary(S,idx)'; % pull out belief uncertainty: entropy
    
    % Perform the regression.
    % the vector stats contains the R2 statistic, the F-statistic and its p-value, and an estimate of the error variance.
    [b,bint,~,~,stats_monet(S,:)] = regress(y,[x1, trial_no, ent]);
    
    betas_monet(S,:) = b'; % set aside betas: each row is a subject, each column for the 3 predictors
    
end

% DF for the F test of the regression are n-1 for the total, n-p for the error term,
% and p-1 for the model, where p is the number of predictors (including the constant).
% We have 3 predictors, so the model df are 3-1 = 2, and n=22, so 22-3=19 for the error.
% ie F(2,19).

%%
% Now we have all the betas and regression values,
% we can do 1-sample t-tests across the betas to see if they differ from zero
% remember, the degrees of freedom are n-1 (but these are given in the stats struct).
[h,p,ci,stats] = ttest(betas_instr(:,2))

%% Do regressions on the ERP components now.
% We can't do the trial-by-trial analyses (or at least, Karen did them but they didn't work).
% So let's look at single values per subject: a regular single regression across subjects.
% REMEMBER: we need to omit 3 subjects, taking us to 19 (not 22 or 23).
% For the modelling, we have subjects 4:25 (n=22).
% For the ERPs, we drop 10, 15 and 20, giving us n=19.

% NOTE: running the regressions with or without the mean reward variable does not change the
% strength of the R2 values or their significance, but it does (obviously) change the Beta values

%% P3as: regression
% with reward magnitude (averaged across trials per subject),
% belief update size (mutual information) and P3 peak amplitude.
% We will do one regression for each electrode.

% Extract peak P3 amplitudes across subjects:
% largest positive peak in 250-550 ms post-feedback
% channels AFz, Fz, FCz, Cz, CPz
% Only interested in bins 5 and 6 (overall instructive vs overall monetary).
ALLERP = pop_geterpvalues( 'C:\Users\RyanDataPC\Documents\Stefan Honours theses\DARPAK\data\erp\current\DARPAK_all_subjects_list_consistent_bins_RM.txt', [ 250 550],...
    5:6, [ 32 37 38 47 48] , 'Baseline', 'none', 'Binlabel', 'on', 'FileFormat', 'wide', 'Filename',...
    'C:\Users\RyanDataPC\Documents\Stefan Honours theses\DARPAK\data\erp\current\mean_P3s_all_subjects_250-550ms_epoch.txt', 'Fracreplace', 'NaN', 'InterpFactor',  1, 'Measure',...
    'peakampbl', 'Neighborhood',  3, 'PeakOnset',  1, 'Peakpolarity', 'positive', 'Peakreplace', 'absolute', 'Resolution',  3, 'SendtoWorkspace',...
    'on' );

N = size(ERP_MEASURES,3);           % should be 19 subjects
num_bins = size(ERP_MEASURES,1);
num_electrodes = size(ERP_MEASURES,2);

% Set up subject ID vector:
IDs = [4:9, 11:14, 16:19, 21:25];

%% Separate regression for each electrode: P3a
% bin1 = instructive, bin2 = monetary. Do bins separately
bn = 2;
x1 = ones(N,1); % need a column of 1s for a regression (the constant term)
chan_names = {'CPz', 'AFz', 'Fz', 'FCz', 'Cz'};
b =[];
stats=[];
%%
for chan = 1:num_electrodes
    
    d=[]; % pull out ERP data for this bin and electrode/channel
    for s=1:N
        d_tmp = squeeze(ERP_MEASURES(bn,chan,s));
        d = [d, d_tmp];
    end
    % Pull out reward values and belief update sizes,
    % but remember, omit subjects 10, 15 and 20
    if bn == 1
        mean_reward = zeros(N,1); % for the instructive feedback, reward is all 0
        mean_update_size = nanmean(model_results.MI.mean_instructive, 2);
        mean_update_size = mean_update_size(ismember(4:25, IDs));
    elseif bn == 2
        mean_reward = nanmean(model_results.reward.mean_monetary,2); % values in dollars
        % Here, get rid of the subjects we are lacking in ERP data for.
        mean_reward = mean_reward(ismember(4:25,IDs));
        mean_update_size = nanmean(model_results.MI.mean_monetary, 2);
        mean_update_size = mean_update_size(ismember(4:25, IDs));
    end
    
    % Perform the regression for this electrode:
    % Remember, the vector stats contains the R2 statistic, the F-statistic and its p-value, and an estimate of the error variance.
    fprintf(chan_names{chan})
    [b(:,chan), ~, ~, ~, stats(chan,:)] = regress(d',[x1, mean_reward, mean_update_size]);
    
    % Omit the reward predictor from the model:
    %[b(:,chan), ~, ~, ~, stats(chan,:)] = regress(d',[x1, mean_update_size]);
    
    % *** another way to do the regression (get p-values for each regressor):
    lm = fitlm([x1, mean_reward, mean_update_size], d, 'linear') %> give t-stats for each regressor
    %anova(lm,'summary') %> gives the f-test of the whole regression
    
end % end of loop across electrodes
b=b';

%% Now do LPPs amplitudes: regression with belief update size

ALLERP = pop_geterpvalues( 'C:\Users\RyanDataPC\Documents\Stefan Honours theses\DARPAK\data\erp\current\DARPAK_all_subjects_list_consistent_bins_RM.txt', [ 550 900],...
    5:6, [ 31 32 48] , 'Baseline', 'none', 'Binlabel', 'on', 'FileFormat', 'wide', 'Filename',...
    'C:\Users\RyanDataPC\Documents\Stefan Honours theses\DARPAK\data\erp\current\mean_LPP_amplitude_all_subjects_550-900ms_epoch.txt', 'Fracreplace', 'NaN', 'InterpFactor',  1, 'Measure',...
    'meanbl', 'PeakOnset',  1, 'Resolution',  3, 'SendtoWorkspace', 'on' );

N = size(ERP_MEASURES,3);           % should be 19 subjects
num_bins = size(ERP_MEASURES,1);
num_electrodes = size(ERP_MEASURES,2);

% Set up subject ID vector:
IDs = [4:9, 11:14, 16:19, 21:25];

%%
% bin1 = instructive, bin2 = monetary. Do bins separately
bn = 2;
x1 = ones(N,1); % need a column of 1s for a regression (the constant term)
chan_names = {'Pz', 'CPz', 'Cz'};
b =[];
stats=[];
%%
for chan = 1:num_electrodes
    
    d=[]; % pull out ERP data for this bin and electrode/channel
    for s=1:N
        d_tmp = squeeze(ERP_MEASURES(bn,chan,s));
        d = [d, d_tmp];
    end
    % Pull out reward values and belief update sizes,
    % but remember, omit subjects 10, 15 and 20
    if bn == 1
        mean_reward = zeros(N,1); % for the instructive feedback, reward is all 0
        mean_update_size = nanmean(model_results.MI.mean_instructive, 2);
        mean_update_size = mean_update_size(ismember(4:25, IDs));
    elseif bn == 2
        mean_reward = nanmean(model_results.reward.mean_monetary,2); % values in dollars
        % Here, get rid of the subjects we are lacking in ERP data for.
        mean_reward = mean_reward(ismember(4:25,IDs));
        mean_update_size = nanmean(model_results.MI.mean_monetary, 2);
        mean_update_size = mean_update_size(ismember(4:25, IDs));
    end
    
    % Perform the regression for this electrode:
    % Remember, the vector stats contains the R2 statistic, the F-statistic and its p-value, and an estimate of the error variance.
    fprintf(chan_names{chan})
    [b(:,chan), ~, ~, ~, stats(chan,:)] = regress(d',[x1, mean_reward, mean_update_size]);
    
    % *** another way to do the regression (get p-values for each regressor):
    lm = fitlm([x1, mean_reward, mean_update_size], d, 'linear') %> give t-stats for each regressor
    %anova(lm,'summary') %> gives the f-test of the whole regression
    
end % end of loop across electrodes
b=b';

%% Regression of FRNs with * belief uncertainty *
% Here we load the FRNs computed using the (modified) version
% of Daniel Bennett's code, getFRN.m
% Remember, the structure of the data is in the form of a 19 row * 30 column matrix
% such that subjects are each in a different row, then electrodes in order 1:5, with each of the six BINS within an electrode
% ie [electrode1([bin 1:6]), electode2([bin 1:6]), ...] forming the 6*5 =30 columns
% *** So this is different to the above format of the other ERP components ***.
load 'C:\Users\RyanDataPC\Documents\Stefan Honours theses\DARPAK\data\erp\current\FRN_amplitudes_for_analysis.mat'

num_bins = 2;
num_electrodes = 5;
N = size(frnAmplitude,1);

% Set up subject ID vector:
IDs = [4:9, 11:14, 16:19, 21:25];

%% Separate regression for each electrode: FRNs
% bin1 = instructive, bin2 = monetary. Do bins separately
bn = 2;
x1 = ones(N,1); % need a column of 1s for a regression (the constant term)
chan_names = {'CPz', 'AFz', 'Fz', 'FCz', 'Cz'};
b =[];
stats=[];
%%
for chan = 1:num_electrodes
    
    % pull out ERP data for this bin and electrode/channel
    % Pull out reward values and belief update sizes,
    % but remember, omit subjects 10, 15 and 20
    if bn == 1
        mean_reward = zeros(N,1); % for the instructive feedback, reward is all 0
        mean_uncert = nanmean(model_results.entropy.mean_instructive, 2);
        mean_uncert = mean_uncert(ismember(4:25, IDs));
        
        d = frnAmplitude(:,5:6:end); % pull out instructive (bin 5) data only, for each electrode
        
    elseif bn == 2
        mean_reward = nanmean(model_results.reward.mean_monetary,2); % values in dollars
        % Here, get rid of the subjects we are lacking in ERP data for.
        mean_reward = mean_reward(ismember(4:25,IDs));
        mean_uncert = nanmean(model_results.entropy.mean_monetary, 2);
        mean_uncert = mean_uncert(ismember(4:25, IDs));
        
        d = frnAmplitude(:,6:6:end); % pull out monetary (bin 6) data only, for each electrode
    end
    
    % Perform the regression for this electrode:
    % Remember, the vector stats contains the R2 statistic, the F-statistic and its p-value, and an estimate of the error variance.
    fprintf(chan_names{chan})
    [b(:,chan), ~, ~, ~, stats(chan,:)] = regress(d(:,chan),[x1, mean_reward, mean_uncert]);
    
    % *** another way to do the regression (get p-values for each regressor):
    lm = fitlm([x1, mean_reward, mean_uncert], d(:,chan), 'linear') %> give t-stats for each regressor
    %anova(lm,'summary') %> gives the f-test of the whole regression
    
end % end of loop across electrodes
b=b';
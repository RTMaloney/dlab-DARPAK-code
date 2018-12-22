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
% remember, the degrees of freedom are n-1

[h,p,ci,stats] = ttest(betas_instr(:,2))

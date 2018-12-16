

n = 22; %we have modelling results for 22 subjects.

%%
% # import the behavioural data excel sheet 
bh_data = importdata('C:\Users\RyanDataPC\Dropbox\DARPAK\material\Behavioural_Data.xlsx');

%%
% # Pull out just the averages for plotting. They are in row 24
% # There are missing data for some subjects after trial 15, so we'll only plot those
rewrd_trials = bh_data.data.TrialByTrial(24,1:15);   %# The 15 reward trials first
instr_trials = bh_data.data.TrialByTrial(24,26:40);  %# The 15 instructive trials

% # Also pull out the standard errors; they are row 25
rewrd_trials_se = bh_data.data.TrialByTrial(25,1:15);
instr_trials_se = bh_data.data.TrialByTrial(25,26:40); 

%%
% Now make the line plots in a new figure
figure(1)
mon_line = plot(rewrd_trials*100, 'r');
hold on
ins_line = plot(instr_trials*100, 'k');

% Convert performance data into a percentage contrast when plotting
errorbar(rewrd_trials*100, rewrd_trials_se*100, 'r', 'LineWidth', 2, 'CapSize', 0)
errorbar(instr_trials*100, instr_trials_se*100, 'k', 'LineWidth', 2,'CapSize', 0)

%Now we can insert this legend
[~, hobj, ~, ~] = legend([mon_line ins_line],{'Monetary feedback', 'Instructive feedback'}, 'FontSize',14)
% Make the linewidth bigger on the legend
hl = findobj(hobj,'type','line');
set(hl,'LineWidth',2)
legend('boxoff')

% Adjust the axis limits
xlim([0 16])
ylim([0 30])
xlabel('Trial Number')
ylabel('Choice error (% contrast)')

% Final formatting:
set(gca, 'box', 'off', 'TickDir', 'out', ...
    'FontSize', 14)

%%  ***** Now load and plot the initial modelling results. *****

% These have been processed from matlab data in the try_to_unpack_bh_data_DARPAK.m script, in the form of a .mat file.

% Load the mat file
model_results = load ('C:\Users\RyanDataPC\Documents\Stefan Honours theses\DARPAK\data\behavioural\DARKPAK_proc_model_results.mat');

% Compute mean MI across subjects for instructive condition
MI_mean_instr = nanmean(model_results.MI.mean_instructive);
MI_mean_monet = nanmean(model_results.MI.mean_monetary);

% And now do the same for entropy
entropy_mean_instr = nanmean(model_results.entropy.mean_instructive);
entropy_mean_monet = nanmean(model_results.entropy.mean_monetary);

% Now compute the standard errors of these values, for plotting.
MI_SE_instr = nanstd(model_results.MI.mean_instructive) / sqrt(n);
MI_SE_monet = nanstd(model_results.MI.mean_monetary) / sqrt(n);

entropy_SE_instr = nanstd(model_results.entropy.mean_instructive) / sqrt(n);
entropy_SE_monet = nanstd(model_results.entropy.mean_monetary) / sqrt(n);

%% Plot the modelling results.
% Only plot the first 15 trials, because they are missing for some subjects for trials >15

% Plot entropy (belief uncertainty).
figure(2)
errorbar(entropy_mean_monet(1:15), entropy_SE_monet(1:15), 'r', 'CapSize', 0, 'LineWidth', 2)
hold on
errorbar(entropy_mean_instr(1:15), entropy_SE_instr(1:15),'k', 'CapSize', 0, 'LineWidth', 2)

% Adjust the axis limits
xlim([0 16])
ylim([4 8])
xlabel('Trial Number')
ylabel('Belief uncertainty (a.u.)')

% Final formatting:
set(gca, 'box', 'off', 'TickDir', 'out', ...
    'FontSize', 14)

%%
% Now plot mutual information (belief update magnitude)
figure(3)
errorbar(MI_mean_monet(1:15), MI_SE_monet(1:15), 'r', 'CapSize', 0, 'LineWidth', 2)
hold on
errorbar(MI_mean_instr(1:15), MI_SE_instr(1:15),'k', 'CapSize', 0, 'LineWidth', 2)

% Adjust the axis limits
xlim([0 16])
ylim([0 0.5])
xlabel('Trial Number')
ylabel('Belief update magnitude (a.u.)')

% Final formatting:
set(gca, 'box', 'off', 'TickDir', 'out', ...
    'FontSize', 14)


% Extract and perform anovas on ERP components across subjects, from the DARPAK study.
% Some simple plots are also produced.
% R Maloney, Dec 2018

%  ---------------------------------------------- %
%% -------------------- P3as -------------------- %%
%  ---------------------------------------------- %

% Extract peak P3 amplitudes across subjects:
% largest positive peak in 250-550 ms post-feedback
% channels AFz, Fz, FCz, Cz, CPz
ALLERP = pop_geterpvalues( 'C:\Users\RyanDataPC\Documents\Stefan Honours theses\DARPAK\data\erp\current\DARPAK_all_subjects_list_consistent_bins_RM.txt', [ 250 550],...
  1:8, [ 32 37 38 47 48] , 'Baseline', 'none', 'Binlabel', 'on', 'FileFormat', 'wide', 'Filename',...
 'C:\Users\RyanDataPC\Documents\Stefan Honours theses\DARPAK\data\erp\current\mean_P3s_all_subjects_250-550ms_epoch.txt', 'Fracreplace', 'NaN', 'InterpFactor',  1, 'Measure',...
 'peakampbl', 'Neighborhood',  3, 'PeakOnset',  1, 'Peakpolarity', 'positive', 'Peakreplace', 'absolute', 'Resolution',  3, 'SendtoWorkspace',...
 'on' );

%%
 %From Daniel B's email:
%I definitely think it’s something that might be worth looking at further
%this time around. I’m not sure how the results will come out, but it could
%be interesting to look at a 2x2 factorial analysis of the P3 amplitude,
%where the factors are condition (rewarding, instructive) and belief update
%size (small, large). Eyeballing the grand average plot in your email,
%there could be an odd (and potentailly interesting) sort of interaction
%going on between these factors, such that the large/small belief updates
%are more separated in the instructive condition than the rewarding condition.
% Presumably, also electrode loc

% Determine number of subjects, bins, electrodes
N = size(ERP_MEASURES,3);           % should be 19 subjects
num_bins = size(ERP_MEASURES,1);
num_electrodes = size(ERP_MEASURES,2);

% Set up subject ID vector:
IDs = [4:9, 11:14, 16:19, 21:25];

%%  First up, lets plot the mean P3a amplitudes across SUBJECTS.
% Pull out the means and standard errors.
for ii=1:num_bins       % loop across BINS
    for jj=1:num_electrodes % loop across ELECTRODES        
        P3a_ampl_means(ii,jj) = mean(squeeze(ERP_MEASURES(ii,jj,:))); % compute mean
        P3a_ampl_SEs(ii,jj) = std(squeeze(ERP_MEASURES(ii,jj,:)))/ sqrt(N); % compute SEM   
    end
end

% Just plot overall instructive/monetary for the 5 electrodes (bins 5 and 6):
figure
plot_mean_ampl = [P3a_ampl_means(5,:); P3a_ampl_means(6,:)]';
plot_SEs_ampl = [P3a_ampl_SEs(5,:); P3a_ampl_SEs(6,:)]';

b1 = bar(plot_mean_ampl,'BarWidth',0.95, 'LineWidth', 1.5);

% Put errorbars on the appropriate locations on the bars. This is not straightforward.
% Find out the number of clusters along the x-axis
ngroups = size(plot_mean_ampl, 1);
% And the number of bars per cluster
nbars = size(plot_mean_ampl, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5)); % I think 0.8 is the default bar width
ln_cols = {'k', 'k'}; % Specify what the error bar colours should be.
% Set the position of each error bar in the centre of the main bar
hold on
for ii = 1:nbars
    % Calculate center of each bar before placing error bars
    jj = [(1:ngroups) - groupwidth/2 + (2*ii-1) * groupwidth / (2*nbars)];
    errorbar(jj, plot_mean_ampl(:,ii), plot_SEs_ampl(:,ii), ln_cols{ii}, 'linestyle', 'none', 'LineWidth', 1.5);
end

% Set color for each bar face 
b1(1).FaceColor = 'k';
b1(2).FaceColor = 'r'; 

% Some formatting for the bar plot. Also put in labels on the X-axis values.
set(gca, 'box', 'off', 'TickDir', 'out', 'YGrid', 'on', 'FontSize', 14, ...
    'XTickLabels', {'CPz','AFz','Fz','FCz','Cz'})
ylabel ('Amplitude (uV)');
xlabel('Electrode site')
title('Mean P3a amplitudes, 250-550 ms post-feedback')
legend('Instructive', 'Monetary') % put in legend for just the data

%% Plot publication-quality bar graph of mean P3a amplitudes averaged across electrodes:
% We will plot these to the left of the data split according to condition/electrode
% (on the same axes). Be sure to use the same figure

% Reassign values to plot:
m = mean(plot_mean_ampl);
se = std(plot_mean_ampl)/ sqrt(num_electrodes); % SE across electrodes
figure
% Set the points on the x-axis as -1 and 0, so the bars are to the left
% of the collapsed data (drawn below).
b2(1) = bar([-1 0],[m(1) nan], 'BarWidth',0.5, 'LineWidth', 2, 'FaceColor', 'k');
hold on
b2(2) = bar([-1 0],[nan m(2)], 'BarWidth',0.5, 'LineWidth', 2, 'FaceColor', 'r');

% place errorbars 
errorbar([-1 0], m, se, 'k', 'linestyle', 'none', 'LineWidth', 2, 'CapSize', 0);
%xlim([0.5 2.5])
ylim([0 20])

% Some formatting for the bar plot. 
set(gca, 'box', 'off', 'TickDir', 'out', 'YGrid', 'on', 'FontSize', 14, ...
    'XTickLabels', {'',''})

%% Now plot P3a peak amplitudes split according to belief update size
figure
% reassign values to plot:
plot_mean_ampl = [P3a_ampl_means(1,:); P3a_ampl_means(2,:); P3a_ampl_means(3,:); P3a_ampl_means(4,:)]';
plot_SEs_ampl = [P3a_ampl_SEs(1,:); P3a_ampl_SEs(2,:); P3a_ampl_SEs(3,:); P3a_ampl_SEs(4,:)]';

b1 = bar(plot_mean_ampl,1, 'Grouped', 'BarWidth', 0.75, 'LineWidth', 2);

% Put errorbars on the appropriate locations on the bars. This is not straightforward.
% Find out the number of clusters along the x-axis
ngroups = size(plot_mean_ampl, 1);
% And the number of bars per cluster
nbars = size(plot_mean_ampl, 2);
% Calculating the width for each bar group
groupwidth = min(1, nbars/(nbars + 1.5)); % I think 0.8 is the default bar width
ln_cols = {'k'}; % Specify what the error bar colours should be.
% Set the position of each error bar in the centre of the main bar
hold on
for ii = 1:nbars
    % Calculate center of each bar before placing error bars
    jj = [(1:ngroups) - groupwidth/2 + (2*ii-1) * groupwidth / (2*nbars)];
    errorbar(jj, plot_mean_ampl(:,ii), plot_SEs_ampl(:,ii), ln_cols{1}, 'linestyle', 'none', 'LineWidth', 2, 'CapSize',0);
end

% Set color for each bar face 
b1(1).FaceColor = 'w';           % Instr, small
b1(2).FaceColor = [0.5 0.5 0.5]; % Instr, large
b1(3).FaceColor = [1 0 0.5];     % Monetary, small
b1(4).FaceColor = [1 0.5 0];     % Monetary, large

% publication-quality formatting:
% set(gca, 'box', 'off', 'TickDir', 'out', 'YGrid', 'on', 'FontSize', 14, 'XTick', [-1:5], 'XTickLabels', {'','','','','','',''})

% Some formatting for the bar plot. Also put in labels on the X-axis values.
set(gca, 'box', 'off', 'TickDir', 'out', 'YGrid', 'on', 'FontSize', 14, ...
    'XTickLabels', {'CPz','AFz','Fz','FCz','Cz'})
ylabel ('Amplitude (uV)');
xlabel('Electrode site')
title('Mean peak P3a amplitude, 250-550 ms post-feedback')
legend('Instructive, small update', 'Instructive, large update', ...
    'Monetary, small update', 'Monetary, large update') % put in legend for just the data

%% Set up variables for the 3-way anova for the P3a.
% These are data from the first 4 bins. 
% Test small/large belief update (2) * feedback condition (2) * electrode (5)
% * N subjects 2*2*5*19 = 380
% So we want one long vector for each item
electrodes = reshape(repmat({'CPz','AFz','Fz','FCz','Cz'},1,2*2*N)', 1, 2*2*num_electrodes*N)';
belief_size = repmat([repmat({'small'},1,num_electrodes), repmat({'large'},1,num_electrodes)], 1, N*2)';
feedback_cond = repmat([repmat({'instructive'},1,num_electrodes*2), repmat({'monetary'},1,num_electrodes*2)], 1, N)';

% Set up subject vector: N=19
% Use the subject codes included in this analysis. Remember, S1-4, 10, and 15 are omitted
subject = reshape( repmat(IDs',1,num_bins/2*num_electrodes)', ... 
    1, num_bins/2*num_electrodes*N)'; 

% Loop across all subjects, pulling out data from bins 1:4,
% and rearranging it so bins*electrodes are stacked in COLUMNS,
% one column for each subject.
d=[];
for s=1:N
    d_tmp = squeeze(ERP_MEASURES(1:num_bins/2,:,s));
    d = [d, reshape(d_tmp', 1, num_bins/2*num_electrodes)'];
end
d = d(:); % rearrange the data, so subjects are stacked on top in one long vector: 1,1,1,....1*20, 2,2,2,...n

%% Perform the 3-way ANOVA on the peak P3a amplitudes
% Now we should have all we need for the anova.
% Use subject as a random factor (stronger mixed-effects model):
[p,tbl,stats] = anovan(d, {subject feedback_cond belief_size electrodes}, ... 
    'model',3,'random', 1, 'varnames',{'ID' 'feedback'  'belief size'  'electrode site'}); %Include 3-way interaction

% If none of the 2-way interactions are significant, then it follows that there is no point doing any 2-way anovas.

% Omit the random factor (weaker model), to allow for multiple comparison tests 
[p,tbl,stats] = anovan(d, {feedback_cond belief_size electrodes}, 'model','full', 'varnames',{'feedback'  'belief size'  'electrode site'}); 
% Do multiple comparisons on electrode site (ie, the 3rd grouping variable).
multcompare(stats, 'Dimension', 3, 'CType', 'bonferroni')

% compute the effect size: partial eta-squared for each factor in the table
%(see here for explanation: https://psychohawks.wordpress.com/2010/10/31/effect-size-for-analysis-of-variables-anova/)
eta_sq = [cell2mat({tbl{2:end-1,2}})./ tbl{end,2}]'; % divide effect SS by total SS


%  ---------------------------------------------- %
%% -------------------- LPPs -------------------- %%
%  ---------------------------------------------- %


% Now, extract the LPP data: mean amplitude
% Extract the MEAN amplitude (ie voltage) in 550-900 ms post feedback for the LPP (late positive potential)
% at electrodes Cz(48), CPz(32) and Pz(31)
clear 
close all
ALLERP = pop_geterpvalues( 'C:\Users\RyanDataPC\Documents\Stefan Honours theses\DARPAK\data\erp\current\DARPAK_all_subjects_list_consistent_bins_RM.txt', [ 550 900],...
  1:8, [ 31 32 48] , 'Baseline', 'none', 'Binlabel', 'on', 'FileFormat', 'wide', 'Filename',...
 'C:\Users\RyanDataPC\Documents\Stefan Honours theses\DARPAK\data\erp\current\mean_LPP_amplitude_all_subjects_550-900ms_epoch.txt', 'Fracreplace', 'NaN', 'InterpFactor',  1, 'Measure',...
 'meanbl', 'PeakOnset',  1, 'Resolution',  3, 'SendtoWorkspace', 'on' );

N = size(ERP_MEASURES,3);           % should be 19 subjects
num_bins = size(ERP_MEASURES,1);
num_electrodes = size(ERP_MEASURES,2);

% Set up subject ID vector:
IDs = [4:9, 11:14, 16:19, 21:25];
%%  plot the mean LPP amplitudes across SUBJECTS.
% Pull out ALL the means and standard errors.
for ii=1:num_bins           % loop across BINS
    for jj=1:num_electrodes % loop across ELECTRODES        
        LPP_ampl_means(ii,jj) = mean(squeeze(ERP_MEASURES(ii,jj,:))); % compute mean
        LPP_ampl_SEs(ii,jj) = std(squeeze(ERP_MEASURES(ii,jj,:)))/ sqrt(N); % compute SEM   
    end
end

% Just plot overall instructive/monetary for the 5 electrodes (bins 5 and 6):
figure(3)
plot_mean_ampl = [LPP_ampl_means(5,:); LPP_ampl_means(6,:)]';
plot_SEs_ampl = [LPP_ampl_SEs(5,:); LPP_ampl_SEs(6,:)]';

b2 = bar(plot_mean_ampl,'BarWidth',0.95, 'LineWidth', 1.5);

% Put errorbars on the appropriate locations on the bars. This is not straightforward.
% Find out the number of clusters along the x-axis
ngroups = size(plot_mean_ampl, 1);
% And the number of bars per cluster
nbars = size(plot_mean_ampl, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5)); % I think 0.8 is the default bar width
ln_cols = {'k'}; % Specify what the error bar colours should be.
% Set the position of each error bar in the centre of the main bar
hold on
for ii = 1:nbars
    % Calculate center of each bar before placing error bars
    jj = [(1:ngroups) - groupwidth/2 + (2*ii-1) * groupwidth / (2*nbars)];
    errorbar(jj, plot_mean_ampl(:,ii), plot_SEs_ampl(:,ii), ln_cols{1}, 'linestyle', 'none', 'LineWidth', 1.5);
end

% Set color for each bar face 
b2(1).FaceColor = 'k';
b2(2).FaceColor = 'r'; 

% Some formatting for the bar plot. Also put in labels on the X-axis values.
set(gca, 'box', 'off', 'TickDir', 'out', 'YGrid', 'on', 'FontSize', 14, ...
    'XTickLabels', {'Pz','CPz','Cz'})
ylabel ('Amplitude (uV)');
xlabel('Electrode site')
title('Mean LPP amplitude in 550-900 ms post-feedback')
legend('Instructive', 'Monetary') % put in legend for just the data

%% Plot publication-quality bar graph of mean LPP amplitudes averaged across electrodes:

m = mean(plot_mean_ampl);
se = std(plot_mean_ampl)/ sqrt(num_electrodes); % SE across electrodes
figure
b2(1) = bar([m(1) nan], 'BarWidth',0.5, 'LineWidth', 1.5, 'FaceColor', 'k');
hold on
b2(2) = bar([nan m(2)], 'BarWidth',0.5, 'LineWidth', 1.5, 'FaceColor', 'r');

% place errorbars 
errorbar([1 2], m, se, 'k', 'linestyle', 'none', 'LineWidth', 1.5);
xlim([0.5 2.5])
ylim([0 20])

% Some formatting for the bar plot. 
set(gca, 'box', 'off', 'TickDir', 'out', 'YGrid', 'on', 'FontSize', 14, ...
    'XTickLabels', {'',''})

%% Now plot LPP mean amplitudes split according to belief update size
figure(4)
plot_mean_ampl = [LPP_ampl_means(1,:); LPP_ampl_means(2,:); LPP_ampl_means(3,:); LPP_ampl_means(4,:)]';
plot_SEs_ampl = [LPP_ampl_SEs(1,:); LPP_ampl_SEs(2,:); LPP_ampl_SEs(3,:); LPP_ampl_SEs(4,:)]';

b3 = bar(plot_mean_ampl,'BarWidth',0.95, 'LineWidth', 1.5);

% Put errorbars on the appropriate locations on the bars. This is not straightforward.
% Find out the number of clusters along the x-axis
ngroups = size(plot_mean_ampl, 1);
% And the number of bars per cluster
nbars = size(plot_mean_ampl, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5)); % I think 0.8 is the default bar width
ln_cols = {'k'}; % Specify what the error bar colours should be.
% Set the position of each error bar in the centre of the main bar
hold on
for ii = 1:nbars
    % Calculate center of each bar before placing error bars
    jj = [(1:ngroups) - groupwidth/2 + (2*ii-1) * groupwidth / (2*nbars)];
    errorbar(jj, plot_mean_ampl(:,ii), plot_SEs_ampl(:,ii), ln_cols{1}, 'linestyle', 'none', 'LineWidth', 1.5);
end

% Set color for each bar face 
b3(1).FaceColor = [0.5 0.5 0.5]; % Instr, small
b3(2).FaceColor = 'k';           % Instr, large
b3(3).FaceColor = [1 0.5 0.5];   % Monetary, small
b3(4).FaceColor = 'r';           % Monetary, large

% for publication-quality:
%set(gca, 'box', 'off', 'TickDir', 'out', 'YGrid', 'on', 'FontSize', 14, 'XTickLabels', {'','',''})
ylim([0 20]) % to match P3a amplitudes

% Some formatting for the bar plot. Also put in labels on the X-axis values.
set(gca, 'box', 'off', 'TickDir', 'out', 'YGrid', 'on', 'FontSize', 14, ...
    'XTickLabels', {'Pz','CPz','Cz'})
ylabel ('Amplitude (uV)');
xlabel('Electrode site')
title('Mean LPP amplitude in 550-900 ms post-feedback')
legend('Instructive, small update', 'Instructive, large update', ...
    'Monetary, small update', 'Monetary, large update') % put in legend for just the data


%% Set up variables for the 3-way anova for LPP
% These are data from the first 4 bins. 
% Test small/large belief update (2) * feedback condition (2) * electrode (3)
% * N subjects 2*2*3*19 = 228
% So we want one long vector for each item
electrodes = reshape(repmat({'Pz','CPz','Cz'},1,2*2*N)', 1, 2*2*num_electrodes*N)';
belief_size = repmat([repmat({'small'},1,num_electrodes), repmat({'large'},1,num_electrodes)], 1, N*2)';
feedback_cond = repmat([repmat({'instructive'},1,num_electrodes*2), repmat({'monetary'},1,num_electrodes*2)], 1, N)';

% Set up subject vector: N=19
% Use the subject codes included in this analysis. Remember, S1-4, 10, and 15 are omitted
subject = reshape( repmat(IDs',1,num_bins/2*num_electrodes)', ... 
    1, num_bins/2*num_electrodes*N)'; 

% Loop across all subjects, pulling out data from bins 1:4,
% and rearranging it so bins*electrodes are stacked in COLUMNS,
% one column for each subject.
d=[];
for s=1:N
    d_tmp = squeeze(ERP_MEASURES(1:num_bins/2,:,s));
    d = [d, reshape(d_tmp', 1, num_bins/2*num_electrodes)'];
end
d = d(:); % rearrange the data, so subjects are stacked on top in one long vector: 1,1,1,....1*20, 2,2,2,...n

%% Perform the 3-way ANOVA for the LPP
% Now we should have all we need for the anova.
% Use subject as a random factor (stronger mixed-effects model):
[p,tbl,stats] = anovan(d, {subject feedback_cond belief_size electrodes}, ... 
    'model',3,'random', 1, 'varnames',{'ID' 'feedback'  'belief size'  'electrode site'}); %Include 3-way interaction

% Perform 3-way anova without the random factor
[p,tbl,stats] = anovan(d, {feedback_cond belief_size electrodes}, ... 
    'model','full', 'varnames',{'feedback'  'belief size'  'electrode site'}); 

% If none of the 2-way interactions are significant, then it follows that there is no point doing any 2-way anovas.
[p,tbl,stats] = anovan(d, {belief_size electrodes}, ... 
    'model','full', 'varnames',{'belief size'  'electrode site'}); 

% compute the effect size: partial eta-squared for each factor in the table
%(see here for explanation: https://psychohawks.wordpress.com/2010/10/31/effect-size-for-analysis-of-variables-anova/)
eta_sq = [cell2mat({tbl{2:end-1,2}})./ tbl{end,2}]'; % divide effect SS by total SS


%  ---------------------------------------------- %
%% -------------------- FRNs -------------------- %%
%  ---------------------------------------------- %

% Now, analyse the FRN amplitudes
% Here we load the FRNs computed using the (modified) version 
% of Daniel Bennett's code, getFRN.m 
% Remember, the structure of the data is in the form of a 19 row * 30 column matrix
% such that subjects are each in a different row, then electrodes in order 1:5, with each of the six BINS within an electrode
% ie [electrode1([bin 1:6]), electode2([bin 1:6]), ...] forming the 6*5 =30 columns
% So this is different to the above format of the other ERP components.
clear
close all
load 'C:\Users\RyanDataPC\Documents\Stefan Honours theses\DARPAK\data\erp\current\FRN_amplitudes_for_analysis.mat'

num_bins = 6;
num_electrodes = 5;
N = size(frnAmplitude,1);

% Set up subject ID vector:
IDs = [4:9, 11:14, 16:19, 21:25];

%% Plot the overall FRNs for bins 5 and 6: not split according to belief update.
% Just plot overall instructive/monetary for the 5 electrodes (bins 5 and 6):
figure(5)
plot_mean_ampl = [mean(frnAmplitude(:,5:6:end)); mean(frnAmplitude(:,6:6:end))]';
plot_SEs_ampl = [std(frnAmplitude(:,5:6:end)) / sqrt(N); ...
                 std(frnAmplitude(:,6:6:end)) / sqrt(N)]';

b3 = bar(plot_mean_ampl,'BarWidth',0.95, 'LineWidth', 1.5);

% Put errorbars on the appropriate locations on the bars. This is not straightforward.
% Find out the number of clusters along the x-axis
ngroups = size(plot_mean_ampl, 1);
% And the number of bars per cluster
nbars = size(plot_mean_ampl, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5)); % I think 0.8 is the default bar width
ln_cols = {'k'}; % Specify what the error bar colours should be. Need to do this because
% black is not the default.
% Set the position of each error bar in the centre of the main bar
hold on
for ii = 1:nbars
    % Calculate center of each bar before placing error bars
    jj = [(1:ngroups) - groupwidth/2 + (2*ii-1) * groupwidth / (2*nbars)];
    errorbar(jj, plot_mean_ampl(:,ii), plot_SEs_ampl(:,ii), ln_cols{1}, 'linestyle', 'none', 'LineWidth', 1.5);
end

% Set color for each bar face
b3(1).FaceColor = 'k';
b3(2).FaceColor = 'r'; 

% Some formatting for the bar plot. Also put in labels on the X-axis values.
set(gca, 'box', 'off', 'TickDir', 'out', 'YGrid', 'on', 'FontSize', 14, ...
    'XTickLabels', {'CPz','AFz','Fz','FCz','Cz'})
ylabel ('Amplitude (uV)');
xlabel('Electrode site')
title('Mean FRN amplitude in 200-550 ms post-feedback')
legend('Instructive', 'Monetary') % put in legend for just the data

%% Plot publication-quality bar graph of mean FRN amplitudes averaged across electrodes:

m = mean(plot_mean_ampl);
se = std(plot_mean_ampl)/ sqrt(num_electrodes); % SE across electrodes
figure
b2(1) = bar([m(1) nan], 'BarWidth',0.5, 'LineWidth', 1.5, 'FaceColor', 'k');
hold on
b2(2) = bar([nan m(2)], 'BarWidth',0.5, 'LineWidth', 1.5, 'FaceColor', 'r');

% place errorbars 
errorbar([1 2], m, se, 'k', 'linestyle', 'none', 'LineWidth', 1.5);
xlim([0.5 2.5])
ylim([0 20])

% Some formatting for the bar plot. 
set(gca, 'box', 'off', 'TickDir', 'out', 'YGrid', 'on', 'FontSize', 14, ...
    'XTickLabels', {'',''})

%% Plot FRNs split according to large/small belief update
figure(6)
plot_mean_ampl = [mean(frnAmplitude(:,1:6:end)); mean(frnAmplitude(:,2:6:end)); ... % bins 1 and 2: Instructive, small & Large
                  mean(frnAmplitude(:,3:6:end)); mean(frnAmplitude(:,4:6:end))]'; % bins 3 and 4: Monetary, small and large
plot_SEs_ampl = [std(frnAmplitude(:,1:6:end)) / sqrt(N); ...
                 std(frnAmplitude(:,2:6:end)) / sqrt(N); ...
                 std(frnAmplitude(:,3:6:end)) / sqrt(N); ...
                 std(frnAmplitude(:,4:6:end)) / sqrt(N)]';

b4 = bar(plot_mean_ampl,'BarWidth',0.95, 'LineWidth', 1.5);

% Put errorbars on the appropriate locations on the bars. This is not straightforward.
% Find out the number of clusters along the x-axis
ngroups = size(plot_mean_ampl, 1);
% And the number of bars per cluster
nbars = size(plot_mean_ampl, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5)); % I think 0.8 is the default bar width
ln_cols = {'k'}; % Specify what the error bar colours should be. Need to do this because
% black is not the default.
% Set the position of each error bar in the centre of the main bar
hold on
for ii = 1:nbars
    % Calculate center of each bar before placing error bars
    jj = [(1:ngroups) - groupwidth/2 + (2*ii-1) * groupwidth / (2*nbars)];
    errorbar(jj, plot_mean_ampl(:,ii), plot_SEs_ampl(:,ii), ln_cols{1}, 'linestyle', 'none', 'LineWidth', 1.5);
end

% Set color for each bar face 
b4(1).FaceColor = [0.5 0.5 0.5]; % Instr, small
b4(2).FaceColor = 'k';           % Instr, large
b4(3).FaceColor = [1 0.5 0.5];   % Monetary, small
b4(4).FaceColor = 'r';           % Monetary, large

% publication-quality formatting:
set(gca, 'box', 'off', 'TickDir', 'out', 'YGrid', 'on', 'FontSize', 14, 'XTickLabels', {'','','','',''})
ylim([0 20])

% Some formatting for the bar plot. Also put in labels on the X-axis values.
set(gca, 'box', 'off', 'TickDir', 'out', 'YGrid', 'on', 'FontSize', 14, ...
    'XTickLabels', {'CPz','AFz','Fz','FCz','Cz'})
ylabel ('Amplitude (uV)');
xlabel('Electrode site')
title('Mean FRN amplitude in 200-550 ms post-feedback')
legend('Instructive, small update', 'Instructive, large update', ...
    'Monetary, small update', 'Monetary, large update')

%% Set up variables for the 3-way anova for FRN
% These are data from the first 4 bins. 
% Test small/large belief update (2) * feedback condition (2) * electrode (5)
% * N subjects 2*2*5*19 = 380
% So we want one long vector for each item
% ***NOTE*** because of the way the FRNs are computed, remember that the structure of the FRN 
% data matrix is different to the P3 and LPP analyses above (ie FRN varies bins within electrodes, rather than electrodes within bins)

electrodes = repmat([repmat({'CPz'},1,num_bins-2), repmat({'AFz'},1,num_bins-2), repmat({'Fz'},1,num_bins-2), ...
         repmat({'FCz'},1,num_bins-2), repmat({'Cz'},1,num_bins-2)], 1, N);
belief_size = repmat([{'small'}, {'large'}],1, 2*N*num_electrodes);
feedback_cond = repmat([repmat({'instructive'},1,2), repmat({'monetary'},1,2)], 1, N*num_electrodes);

% Set up subject vector: N=19
% Use the subject codes included in this analysis. Remember, S1-4, 10, and 15 are omitted
subject = reshape( repmat(IDs',1,4*num_electrodes)', ... 
    1, 4*num_electrodes*N);

% Loop across all subjects, pulling out data from bins 1:4,
d=[];
% We want to pull out bins 1:4 for each consecutive electrode:
for s=1:N
    for eth = 1:num_electrodes
        idx = 5*eth-(5-eth);
        d_tmp = frnAmplitude(s, idx:idx+3); % pull out data for bins1:4, eth electrode.
        d = [d, d_tmp];                     % Set aside data into one long vector
    end
end

%% Perform the 3-way ANOVA for the FRN
% Now we should have all we need for the anova.
% Use subject as a random factor (stronger mixed-effects model):
[p,tbl,stats] = anovan(d, {subject feedback_cond belief_size electrodes}, ... 
    'model','full','random', 1, 'varnames',{'ID' 'feedback'  'belief size'  'electrode site'}); %Include 3-way interaction

% Omit the random factor (weaker model), to allow for multiple comparison tests 
[p,tbl,stats] = anovan(d, {feedback_cond belief_size electrodes}, 'model', 3, 'varnames',{'feedback'  'belief size'  'electrode site'}); 

% compute the effect size: partial eta-squared for each factor in the table
%(see here for explanation: https://psychohawks.wordpress.com/2010/10/31/effect-size-for-analysis-of-variables-anova/)
eta_sq = [cell2mat({tbl{2:end-1,2}})./ tbl{end,2}]'; % divide effect SS by total SS

%% Do 2-way anova for FRN, over simple feedback averages
%   *** Not significant ***
% These are data from the last 2 bins of 6. 
% Test feedback condition (2) * electrode (5)
% * N subjects 2*5*19 = 190
% So we want one long vector for each item
% ***NOTE*** because of the way the FRNs are computed, remember that the structure of the FRN 
% data matrix is different to the P3 and LPP analyses above (ie FRN varies bins within electrodes, rather than electrodes within bins)

electrodes = repmat([repmat({'CPz'},1,2), repmat({'AFz'},1,2), repmat({'Fz'},1,2), ...
             repmat({'FCz'},1,2), repmat({'Cz'},1,2)], 1, N);
feedback_cond = repmat([{'instructive'},{'monetary'}], 1, N*num_electrodes);

% Set up subject vector: N=19
% Use the subject codes included in this analysis. Remember, S1-3, 10, 15 and 20 are omitted
subject = reshape( repmat(IDs',1,2*num_electrodes)', ... 
    1, 2*num_electrodes*N);

% Loop across all subjects, pulling out data from bins 1:4,
d=[];
% We want to pull out bins 1:4 for each consecutive electrode:
for s=1:N
    for eth = 1:num_electrodes
        idx = eth*5 + (eth-1);
        d_tmp = frnAmplitude(s, idx:idx+1); % pull out data for bins5:6, eth electrode.
        d = [d, d_tmp];                     % Set aside data into one long vector
    end
end

% Now we have all we need for the 2-way ANOVA. 
[p,tbl,stats] = anovan(d, {feedback_cond electrodes}, ... 
    'model','full', 'varnames',{'feedback' 'electrode site'}); %Include 2-way interaction

% compute the effect size: partial eta-squared for each factor in the table
eta_sq = [cell2mat({tbl{2:end-1,2}})./ tbl{end,2}]'; % divide effect SS by total SS

%% 2-way anova over simple mean amplitudes for the 2 feedback conditions, and electrode site
% *** not using *** 
% Test feedback condition (2) * electrode (5)
% * N subjects 2*5*19 = 190
% Need to re-define all the grouping variables here:
% electrodes = reshape(repmat({'CPz','AFz','Fz','FCz','Cz'},1,2*N)', 1, 2*num_electrodes*N)';
% feedback_cond = repmat([repmat({'instructive'},1,num_electrodes), repmat({'monetary'},1,num_electrodes)], 1, N)';
% 
% Loop across all subjects, pulling out data from bins 5:6,
% and rearranging it so bins*electrodes are stacked in COLUMNS,
% one column for each subject.
% d=[];
% for s=1:N
%     d_tmp = squeeze(ERP_MEASURES(5:6,:,s));
%     d = [d, reshape(d_tmp', 1, 2*num_electrodes)'];
% end
% d = d(:); % rearrange the data, so subjects are stacked on top in one long vector: 1,1,1,....1*20, 2,2,2,...n
% 
% Set up subject vector: N=19
% Use the subject codes included in this analysis. Remember, S1-4, 10, and 15 are omitted
% subject = reshape( repmat([5:9, 11:14, 16:25]',1,2*num_electrodes)', ... 
%     1, 2*num_electrodes*N)'; 
% 
% Now we have all we need for the 2-way ANOVA. 
% [p,tbl,stats] = anovan(d, {subject feedback_cond electrodes}, ... 
%     'model',2,'random', 1, 'varnames',{'ID' 'feedback' 'electrode site'}); %Include 2-way interaction
% Main effect of electrode site only.


%% *** Alternatively ***
% *** not using *** : the data were wrong (2 subjects swapped).
% Fit 3-way repeated measures model using 'fitrm' and 'ranova'
% This allows us to do multiple-comparisons tests of the main effects
% (slightly more confusing way). See here:
% https://au.mathworks.com/matlabcentral/answers/140799-3-way-repeated-measures-anova-pairwise-comparisons-using-multcompare

% we already have our table containing the responses in the saved ERPlab output:
%importdata 'C:\Users\RyanDataPC\Documents\Stefan Honours theses\DARPAK\data\erp\current\mean_P3s_all_subjects_250-550ms_epoch.txt' 

% meanP3sallsubjects250550msepoch = meanP3sallsubjects250550msepoch(:,1:20); % just look at first 4 bins
% meanP3sallsubjects250550msepoch.Properties.VariableNames; %where the variable names can be found
% 
% Create a table reflecting the within subject factors 'feedback condition', 'belief size', and 'electrode site' and their levels
% factorNames = {'feedback',  'belief_size',  'electrode_site'};
% within = table([repmat({'instructive'},1,num_electrodes*2), repmat({'monetary'},1,num_electrodes*2)]', ... % feedback
%     repmat([repmat({'small'},1,num_electrodes), repmat({'large'},1,num_electrodes)], 1, 2)', ...           % belief size
%     repmat({'CPz','AFz','Fz','FCz','Cz'},1,num_bins/2)', ...                                               % electrode site
%     'VariableNames',factorNames);
% 
% fit the repeated measures model
% rm = fitrm(meanP3sallsubjects250550msepoch, ...
%     [meanP3sallsubjects250550msepoch.Properties.VariableNames{1}, '-', ... %Specify all variables
%     meanP3sallsubjects250550msepoch.Properties.VariableNames{20}, '~1'], ...
%     'WithinDesign',within);
% 
% run my repeated measures anova here
% [ranovatbl,A,C,D] = ranova(rm,'WithinModel', 'feedback*belief_size*electrode_site');
% 
% multiple comparisons tests for the main effects
% multcompare(rm, 'electrode_site', 'ComparisonType','Bonferroni') %
% multcompare(rm,'feedback')
% multcompare(rm,'belief_size')

% make pairwise comparisons for the two-way interactions
%multcompare(rm,'feedback','By','electrode_site')




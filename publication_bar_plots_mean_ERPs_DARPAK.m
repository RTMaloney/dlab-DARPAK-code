% This code loads and plots mean amplitudes for ERP components in 
% the DARPAK study: these are publication-quality and indented for import into 
% Adobe Illustrator.
% R Maloney, Jan 2019

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


% Determine number of subjects, bins, electrodes
N = size(ERP_MEASURES,3);           % should be 19 subjects
num_bins = size(ERP_MEASURES,1);
num_electrodes = size(ERP_MEASURES,2);

% Set up subject ID vector:
IDs = [4:9, 11:14, 16:19, 21:25];

%%
% Pull out the means and standard errors.
for ii=1:num_bins       % loop across BINS
    for jj=1:num_electrodes % loop across ELECTRODES
        P3a_ampl_means(ii,jj) = mean(squeeze(ERP_MEASURES(ii,jj,:))); % compute mean
        P3a_ampl_SEs(ii,jj) = std(squeeze(ERP_MEASURES(ii,jj,:)))/ sqrt(N); % compute SEM
    end
end

% Pull out means for bins 5 and 6 only:
plot_mean_ampl = [P3a_ampl_means(5,:); P3a_ampl_means(6,:)]';

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

%% Now plot the mean amplitudes collapsed across update size and electrode:
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
set(gca, 'box', 'off', 'TickDir', 'out', 'YGrid', 'on', 'FontSize', 14, 'XTick', [-1:5], 'XTickLabels', {'','','','','','',''})
xlim([-1.5 5.5])

%  ---------------------------------------------- %
%% -------------------- LPPs -------------------- %%
%  ---------------------------------------------- %

% Now, extract the LPP data: mean amplitude
% Extract the MEAN amplitude (ie voltage) in 550-900 ms post feedback for the LPP (late positive potential)
% at electrodes Cz(48), CPz(32) and Pz(31)
clear
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
plot_mean_ampl = [LPP_ampl_means(5,:); LPP_ampl_means(6,:)]';
plot_SEs_ampl = [LPP_ampl_SEs(5,:); LPP_ampl_SEs(6,:)]';

%% Plot publication-quality bar graph of mean LPP amplitudes averaged across electrodes:

m = mean(plot_mean_ampl);
se = std(plot_mean_ampl)/ sqrt(num_electrodes); % SE across electrodes
figure
b2(1) = bar([-1 0],[m(1) nan], 'BarWidth',0.5, 'LineWidth', 2, 'FaceColor', 'k');
hold on
b2(2) = bar([-1 0],[nan m(2)], 'BarWidth',0.5, 'LineWidth', 2, 'FaceColor', 'r');

% place errorbars
errorbar([-1 0], m, se, 'k', 'linestyle', 'none', 'LineWidth', 2,'CapSize',0);
ylim([0 20])

%xlim([-1.5 5.5])

%% Now plot LPP mean amplitudes split according to belief update size
% reassign values:
plot_mean_ampl = [LPP_ampl_means(1,:); LPP_ampl_means(2,:); LPP_ampl_means(3,:); LPP_ampl_means(4,:)]';
plot_SEs_ampl = [LPP_ampl_SEs(1,:); LPP_ampl_SEs(2,:); LPP_ampl_SEs(3,:); LPP_ampl_SEs(4,:)]';

b3 = bar(plot_mean_ampl,1, 'Grouped', 'BarWidth',0.75, 'LineWidth', 2);

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
b3(1).FaceColor = 'w';           % Instr, small
b3(2).FaceColor = [0.5 0.5 0.5]; % Instr, large
b3(3).FaceColor = [1 0 0.5];     % Monetary, small
b3(4).FaceColor = [1 0.5 0];     % Monetary, large

% publication-quality formatting:
set(gca, 'box', 'off', 'TickDir', 'out', 'YGrid', 'on', 'FontSize', 14, 'XTick', [-1:3], 'XTickLabels', {'','','','',''})
xlim([-1.5 3.5])

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
load 'C:\Users\RyanDataPC\Documents\Stefan Honours theses\DARPAK\data\erp\current\FRN_amplitudes_for_analysis.mat'

num_bins = 6;
num_electrodes = 5;
N = size(frnAmplitude,1);

% Set up subject ID vector:
IDs = [4:9, 11:14, 16:19, 21:25];

% Extract data for plotting:
plot_mean_ampl = [mean(frnAmplitude(:,5:6:end)); mean(frnAmplitude(:,6:6:end))]';
plot_SEs_ampl = [std(frnAmplitude(:,5:6:end)) / sqrt(N); ...
    std(frnAmplitude(:,6:6:end)) / sqrt(N)]';

%% Plot publication-quality bar graph of mean FRN amplitudes averaged across electrodes:

m = mean(plot_mean_ampl);
se = std(plot_mean_ampl)/ sqrt(num_electrodes); % SE across electrodes
figure
b2(1) = bar([-1, 0], [m(1) nan], 'BarWidth',0.5, 'LineWidth', 2, 'FaceColor', 'k');
hold on
b2(2) = bar([-1, 0], [nan m(2)], 'BarWidth',0.5, 'LineWidth', 2, 'FaceColor', 'r');

% place errorbars
errorbar([-1 0], m, se, 'k', 'linestyle', 'none', 'LineWidth', 2, 'CapSize', 0);
ylim([0 20])

%%
% reassign values to plot:
plot_mean_ampl = [mean(frnAmplitude(:,1:6:end)); mean(frnAmplitude(:,2:6:end)); ... % bins 1 and 2: Instructive, small & Large
    mean(frnAmplitude(:,3:6:end)); mean(frnAmplitude(:,4:6:end))]'; % bins 3 and 4: Monetary, small and large
plot_SEs_ampl = [std(frnAmplitude(:,1:6:end)) / sqrt(N); ...
    std(frnAmplitude(:,2:6:end)) / sqrt(N); ...
    std(frnAmplitude(:,3:6:end)) / sqrt(N); ...
    std(frnAmplitude(:,4:6:end)) / sqrt(N)]';

b4 = bar(plot_mean_ampl,1,'Grouped','BarWidth',0.75, 'LineWidth', 2);

% Put errorbars on the appropriate locations on the bars. This is not straightforward.
% Find out the number of clusters along the x-axis
ngroups = size(plot_mean_ampl, 1);
% And the number of bars per cluster
nbars = size(plot_mean_ampl, 2);
% Calculating the width for each bar group
groupwidth = min(1, nbars/(nbars + 1.5)); % I think 0.8 is the default bar width
ln_cols = {'k'}; % Specify what the error bar colours should be. Need to do this because
% black is not the default.
% Set the position of each error bar in the centre of the main bar
hold on
for ii = 1:nbars
    % Calculate center of each bar before placing error bars
    jj = [(1:ngroups) - groupwidth/2 + (2*ii-1) * groupwidth / (2*nbars)];
    errorbar(jj, plot_mean_ampl(:,ii), plot_SEs_ampl(:,ii), ln_cols{1}, 'linestyle', 'none', 'LineWidth', 2, 'CapSize',0);
end

% Set color for each bar face
b4(1).FaceColor = 'w';           % Instr, small
b4(2).FaceColor = [0.5 0.5 0.5]; % Instr, large
b4(3).FaceColor = [1 0 0.5];     % Monetary, small
b4(4).FaceColor = [1 0.5 0];     % Monetary, large

% publication-quality formatting:
set(gca, 'box', 'off', 'TickDir', 'out', 'YGrid', 'on', 'FontSize', 14, 'XTick', [-1:5], 'XTickLabels', {'','','','','','',''})
xlim([-1.5 5.5])

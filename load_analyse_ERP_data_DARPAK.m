%%

% Extract peak P3 amplitudes across subjects:
% largest positive peak in 250-550 ms post-feedback
% channels AFz, Fz, FCz, Cz, CPz
ALLERP = pop_geterpvalues( 'C:\Users\RyanDataPC\Documents\Stefan Honours theses\DARPAK\data\erp\current\P3s_all_subjects_list.txt', [ 250 550],...
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


%% 
% First up, lets plot the mean P3a amplitudes across SUBJECTS.
% Pull out the means and standard errors.
for ii=1:num_bins       % loop across BINS
    for jj=1:num_electrodes % loop across ELECTRODES        
        P3a_ampl_means(ii,jj) = mean(squeeze(ERP_MEASURES(ii,jj,:))); % compute mean
        P3a_ampl_SEs(ii,jj) = std(squeeze(ERP_MEASURES(ii,jj,:)))/ sqrt(N); % compute SEM   
    end
end

% Just plot overall instructive/monetary for the 5 electrodes:
figure
b1 = bar([P3a_ampl_means(5,:); nan(1,5)]', 1, 'FaceColor', 'k'); % plot instuctive
hold on
b2 = bar([nan(1,5); P3a_ampl_means(6,:)]', 1, 'FaceColor', 'r'); % plot monetary
legend([b1(1), b2(end)], 'Instructive', 'Monetary') % put in legend for just the data
errorbar([P3a_ampl_means(5,:); nan(1,5)]', [P3a_ampl_SEs(5,:); nan(1,5)]', 'k', 'LineStyle', 'none')


set(gca, 'TickDir', 'out', ...
'XTickLabel', {'CPz','AFz','Fz','FCz','Cz'})


%%
% Set up variables for the 3-way anova.
% These are data from the first 4 bins. 
% Test small/large belief update (2) * feedback condition (2) * electrode (5)
% * N subjects 2*2*5*19 = 380
% So we want one long vector for each item
electrodes = reshape(repmat({'CPz','AFz','Fz','FCz','Cz'},1,2*2*N)', 1, 2*2*num_electrodes*N)';
belief_size = repmat([repmat({'small'},1,num_electrodes), repmat({'large'},1,num_electrodes)], 1, N*2)';
feedback_cond = repmat([repmat({'instructive'},1,num_electrodes*2), repmat({'monetary'},1,num_electrodes*2)], 1, N)';

% Set up subject vector: N=19
% Use the subject codes included in this analysis. Remember, S1-4, 10, and 15 are omitted
subject = reshape( repmat([5:9, 11:14, 16:25]',1,num_bins/2*num_electrodes)', ... 
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

%% 
% Now we should have all we need for the anova.
% Use subject as a random factor (stronger mixed-effects model):
[p,tbl,stats] = anovan(d, {subject feedback_cond belief_size electrodes}, ... 
    'model',3,'random', 1, 'varnames',{'ID' 'feedback'  'belief size'  'electrode site'}); %Include 3-way interaction

% If none of the 2-way interactions are significant, then it follows that there is no point doing any 2-way anovas.

% Do some t-tests to look at main effects in more depth.


% Omit the random factor (weaker model):
%[p,tbl,stats] = anovan(d, {feedback_cond belief_size electrodes}, 'model',3, 'varnames',{'feedback'  'belief size'  'electrode site'}); 

% compute the effect size: partial eta-squared for each factor in the table
%(see here for explanation: https://psychohawks.wordpress.com/2010/10/31/effect-size-for-analysis-of-variables-anova/)
eta_sq = [cell2mat({tbl{2:end-1,2}})./ tbl{end,2}]'; % divide effect SS by total SS


%% *** Alternatively ***
% Fit 3-way repeated measures model using 'fitrm' and 'ranova'
% (slightly more confusing way). See here:
% https://au.mathworks.com/matlabcentral/answers/140799-3-way-repeated-measures-anova-pairwise-comparisons-using-multcompare

% we already have our table containing the responses in the saved ERPlab output:
%importdata 'C:\Users\RyanDataPC\Documents\Stefan Honours theses\DARPAK\data\erp\current\mean_P3s_all_subjects_250-550ms_epoch.txt' 
meanP3sallsubjects250550msepoch = meanP3sallsubjects250550msepoch(:,1:20); % just look at first 4 bins
meanP3sallsubjects250550msepoch.Properties.VariableNames; %where the variable names can be found

% Create a table reflecting the within subject factors 'feedback condition', 'belief size', and 'electrode site' and their levels
factorNames = {'feedback',  'belief_size',  'electrode_site'};
within = table([repmat({'instructive'},1,num_electrodes*2), repmat({'monetary'},1,num_electrodes*2)]', ... % feedback
    repmat([repmat({'small'},1,num_electrodes), repmat({'large'},1,num_electrodes)], 1, 2)', ...           % belief size
    repmat({'CPz','AFz','Fz','FCz','Cz'},1,num_bins/2)', ...                                               % electrode site
    'VariableNames',factorNames);

% fit the repeated measures model
rm = fitrm(meanP3sallsubjects250550msepoch, ...
    [meanP3sallsubjects250550msepoch.Properties.VariableNames{1}, '-', ... %Specify all variables
    meanP3sallsubjects250550msepoch.Properties.VariableNames{20}, '~1'], ...
    'WithinDesign',within);

% run my repeated measures anova here
[ranovatbl,A,C,D] = ranova(rm,'WithinModel', 'feedback*belief_size*electrode_site')

% make pairwise comparisons for the two-way interactions


%%
% Extract peak POSITIVE amplitudes across subjects:
% This is to calculate the height of the FRN in terms of the peak-to-peak distance between
% the most negative peak in 200-550 ms epoch post-feedback and the 
% immediately preceding positive peak (for which I will set a 0-550 ms epoch). 
% channels AFz, Fz, FCz, CPz 
ALLERP = pop_geterpvalues( 'C:\Users\RyanDataPC\Documents\Stefan Honours theses\DARPAK\data\erp\current\P3s_all_subjects_list.txt', [ 0 550],...
  1:8, [ 32 37 38 47 48] , 'Baseline', 'none', 'Binlabel', 'on', 'FileFormat', 'wide', 'Filename',...
 'C:\Users\RyanDataPC\Documents\Stefan Honours theses\DARPAK\data\erp\current\all_subjects_0-550ms_peak_positive_ampl_for_FRN.xls', 'Fracreplace', 'NaN', 'InterpFactor',  1, 'Measure',...
 'peakampbl', 'Neighborhood',  3, 'PeakOnset',  1, 'Peakpolarity', 'positive', 'Peakreplace', 'absolute', 'Resolution',  3, 'SendtoWorkspace',...
 'on' );

% Now, for the above calculation, do the same but extract the most NEGATIVE peak
% in the window 200-550 ms 
% channels AFz, Fz, FCz, CPz
ALLERP = pop_geterpvalues( 'C:\Users\RyanDataPC\Documents\Stefan Honours theses\DARPAK\data\erp\current\P3s_all_subjects_list.txt', [ 200 550],...
  1:8, [ 32 37 38 47 48] , 'Baseline', 'none', 'Binlabel', 'on', 'FileFormat', 'wide', 'Filename',...
 'C:\Users\RyanDataPC\Documents\Stefan Honours theses\DARPAK\data\erp\current\all_subjects_200-550ms_peak_negative_ampl_for_FRN.xls', 'Fracreplace', 'NaN', 'InterpFactor',  1,...
 'Measure', 'peakampbl', 'Neighborhood',  3, 'PeakOnset',  1, 'Peakpolarity', 'negative', 'Peakreplace', 'absolute', 'Resolution',  3, 'SendtoWorkspace',...
 'on' );

% Extract the MEAN amplitude (ie voltage) in 550-900 ms post feedback for the LPP (late positive potential)
% at electrodes Cz, CPz and Pz
ALLERP = pop_geterpvalues( 'C:\Users\RyanDataPC\Documents\Stefan Honours theses\DARPAK\data\erp\current\P3s_all_subjects_list.txt', [ 550 900],...
  1:8, [ 31 32 48] , 'Baseline', 'none', 'Binlabel', 'on', 'FileFormat', 'wide', 'Filename',...
 'C:\Users\RyanDataPC\Documents\Stefan Honours theses\DARPAK\data\erp\current\mean_LPP_amplitude_all_subjects_550-900ms_epoch.xls', 'Fracreplace', 'NaN', 'InterpFactor',  1, 'Measure',...
 'meanbl', 'PeakOnset',  1, 'Resolution',  3, 'SendtoWorkspace', 'on' );
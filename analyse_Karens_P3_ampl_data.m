
% Load and analyse the P3 amplitude data, from an Excel file,
% and probably what was put in Karen Sasmita's thesis.
% R Maloney, 


% get the data
P3ampl = importdata ('C:\Users\RyanDataPC\Documents\Stefan Honours theses\DARPAK\data\erp\current\karens_P3_amp_data.xlsx')

N = size(P3ampl,1);           % should be 19 subjects
num_bins = 2;
num_electrodes = 4;


% Set up variables for the 2-way anova 
% Test  feedback condition (2) * electrode (4)
% * N subjects 2*4*19 = 152
% So we want one long vector for each item

electrodes = repmat([repmat({'Pz'},1,2), repmat({'CPz'},1,2), ...
             repmat({'FCz'},1,2), repmat({'Cz'},1,2)], 1, N);
feedback_cond = repmat([{'monetary'}, {'instructive'}], 1, N*num_electrodes);
subject = reshape( repmat([4:9, 11:14, 16:19, 21:25]',1,2*num_electrodes)', ... 
    1, 2*num_electrodes*N);

P3ampl = reshape(P3ampl', 1, num_bins*num_electrodes*N);

%% Perform the 2-way ANOVA on the peak P3a amplitudes
% Now we should have all we need for the anova.
% Use subject as a random factor (stronger mixed-effects model):
[p,tbl,stats] = anovan(P3ampl, {subject feedback_cond electrodes}, ... 
    'model',2,'random', 1, 'varnames',{'ID' 'feedback' 'electrode site'}); %Include 2-way interaction

% Omit the random factor (weaker model), to allow for multiple comparison tests 
[p,tbl,stats] = anovan(P3ampl, {feedback_cond electrodes}, 'model','full', 'varnames',{'feedback' 'electrode site'}); 

% Do multiple comparisons on electrode site (ie, the 3rd grouping variable).
multcompare(stats, 'Dimension', 1, 'CType', 'bonferroni')



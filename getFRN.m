% Script to extract Feedback-related negativity measure.
% Written by D Bennett
% Modified by R Maloney, Dec 2018


%% specify analysis
analysisChannels = [32 37 38 47 48]; %[32 33 37 38 47 48];
analysisBins = [1:6];
baseline = [-500 0];
erpSets = 1:19;
frnWindow = [200 550];

%% Load the data
[ERP ALLERP] = pop_loaderp( 'filename', {'K4_precsd_RM.erp', 'K5_precsd.erp', 'K6_precsd.erp', 'K7_precsd.erp', 'K8_precsd.erp', 'K9_precsd.erp', 'K11_precsd.erp',...
    'K12_precsd.erp', 'K13_precsd.erp', 'K14_precsd.erp', 'K16_precsd.erp', 'K17_precsd.erp', 'K18_precsd.erp', 'K19_precsd.erp',...
    'K21_precsd.erp', 'K22_precsd.erp', 'K23_precsd.erp', 'K24_precsd.erp', 'K25_precsd.erp'}, 'filepath',...
    'C:\Users\RyanDataPC\Documents\Stefan Honours theses\DARPAK\data\erp\current\' );

%% load FRN amplitude
% Note that the 'Filename' variable is what the results are saved as (not loaded from).
ALLERP = pop_geterpvalues( ALLERP, frnWindow, analysisBins, analysisChannels , 'Baseline', baseline, 'Binlabel', 'on', 'Erpsets',  erpSets,...
    'FileFormat', 'long', 'Filename', 'C:\Users\RyanDataPC\Documents\Stefan Honours theses\DARPAK\data\erp\current\all_subjects_200-550ms_peak_negative_ampl_for_FRN.txt', 'Fracreplace', 'NaN', 'InterpFactor',  1, 'Measure', 'peakampbl',...
    'Neighborhood',  3, 'Peakpolarity', 'negative', 'Peakreplace', 'absolute', 'Resolution',  3, 'SendtoWorkspace', 'on' );

% The resulting data is a 3D array, ERP_MEASURES, which is 6 rows (bins) * 5 columns (electrodes) * 19 (subjects).

% rearrange the amplitude data,
% so that we have each subject in a different row, then we have electrodes in order 1:5, with each of the six BINS within an electrode
% ie [electrode1([bin 1:6]), electode2([bin 1:6]), ...]
% So bins increment within electrodes, and then electrodes increment
negativePeakAmplitude = reshape(ERP_MEASURES,[numel(analysisChannels) * numel(analysisBins), numel(erpSets)])';


%% get FRN latency
ALLERP = pop_geterpvalues( ALLERP, frnWindow, analysisBins, analysisChannels , 'Baseline', baseline, 'Binlabel', 'on', 'Erpsets',  erpSets,...
    'FileFormat', 'long', 'Filename', 'C:\Users\RyanDataPC\Desktop\x1.txt', 'Fracreplace', 'NaN', 'InterpFactor',  1, 'Measure', 'peaklatbl',...
    'Neighborhood',  3, 'Peakpolarity', 'negative', 'Peakreplace', 'absolute', 'Resolution',  3, 'SendtoWorkspace', 'on' );

% As with the amplitudes above,
% the resulting data is a 3D array, ERP_MEASURES, which is 6 rows (bins) * 5 columns (electrodes) * 19 (subjects).

% rearrange the latency data, as with the amplitude data above.
negativePeakLatency = reshape(ERP_MEASURES,[numel(analysisChannels) * numel(analysisBins), numel(erpSets)])';

%%
% Now compute the preceding positive peak, based on the negative peak
% latencies computed above.
% This places the data in the same format as the other variables:
% namely, subjects in a different row, then electrodes in order 1:5, with each of the six BINS within an electrode
% ie [electrode1([bin 1:6]), electode2([bin 1:6]), ...]
precedingPositivePeak = nan(size(negativePeakAmplitude));
for set = 1:numel(erpSets) % loop across all subjects
    
    chanCounter = 1; % chanCounter indexes the columns, and increments with each BIN (for a given ELECTRODE)  
    for chan = 1:numel(analysisChannels) % loop across electrodes      
        
        for bn = 1:numel(analysisBins)   % loop across bins
            
            ALLERP = pop_geterpvalues( ALLERP, [0 negativePeakLatency(set,chanCounter)], analysisBins(bn), analysisChannels(chan) , 'Baseline', baseline, 'Binlabel', 'on', 'Erpsets',  erpSets(set),...
                'FileFormat', 'long', 'Filename', 'C:\Users\RyanDataPC\Desktop\x1.txt', 'Fracreplace', 'NaN', 'InterpFactor',  1, 'Measure', 'peakampbl',...
                'Neighborhood',  3, 'Peakpolarity', 'positive', 'Peakreplace', 'absolute', 'Resolution',  3, 'SendtoWorkspace', 'on' );
            precedingPositivePeak(set,chanCounter) = ERP_MEASURES; % set aside the amplitude
            chanCounter = chanCounter + 1;
        end
    end
end

%% Compute the FRN amplitude! A difference in the 2 peaks
% Ideally, the resulting differences will all be positive
frnAmplitude = precedingPositivePeak - negativePeakAmplitude;

% Only a tiny minority are negative, so just make them all positive:
frnAmplitude = abs(frnAmplitude);

% Save this resulting matrix
save('C:\Users\RyanDataPC\Documents\Stefan Honours theses\DARPAK\data\erp\current\FRN_amplitudes_for_analysis.mat', 'frnAmplitude')



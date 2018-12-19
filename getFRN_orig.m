% Daniel Bennett's original code for extracting the FRN amplitudes in the DARPAK study.


% specify analysis 
analysisChannels = [32 33 37 38 47 48];
analysisBins = [5 6];
baseline = [-500 0];
erpSets = 1:19;
frnWindow = [200 550];

% get FRN amplitude
ALLERP = pop_geterpvalues( ALLERP, frnWindow, analysisBins, analysisChannels , 'Baseline', baseline, 'Binlabel', 'on', 'Erpsets',  erpSets,...
 'FileFormat', 'long', 'Filename', 'C:\Users\dbennett1\Desktop\x1.xls', 'Fracreplace', 'NaN', 'InterpFactor',  1, 'Measure', 'peakampbl',...
 'Neighborhood',  3, 'Peakpolarity', 'negative', 'Peakreplace', 'absolute', 'Resolution',  3, 'SendtoWorkspace', 'on' );
negativePeakAmplitude = reshape(ERP_MEASURES,[numel(analysisChannels) * 2, numel(erpSets)])';

% get FRN latency
ALLERP = pop_geterpvalues( ALLERP, frnWindow, analysisBins, analysisChannels , 'Baseline', baseline, 'Binlabel', 'on', 'Erpsets',  erpSets,...
 'FileFormat', 'long', 'Filename', 'C:\Users\dbennett1\Desktop\x1.xls', 'Fracreplace', 'NaN', 'InterpFactor',  1, 'Measure', 'peaklatbl',...
 'Neighborhood',  3, 'Peakpolarity', 'negative', 'Peakreplace', 'absolute', 'Resolution',  3, 'SendtoWorkspace', 'on' );
negativePeakLatency = reshape(ERP_MEASURES,[numel(analysisChannels) * 2, numel(erpSets)])';

precedingPositivePeak = nan(size(negativePeakAmplitude));

for set = 1:numel(erpSets)
    chanCounter = 0;
    for chan = 1:numel(analysisChannels)
        
        chanCounter = chanCounter + 1;
        
        ALLERP = pop_geterpvalues( ALLERP, [0 negativePeakLatency(set,chanCounter)], analysisBins(1), analysisChannels(chan) , 'Baseline', baseline, 'Binlabel', 'on', 'Erpsets',  erpSets(set),...
            'FileFormat', 'long', 'Filename', 'C:\Users\dbennett1\Desktop\x1.xls', 'Fracreplace', 'NaN', 'InterpFactor',  1, 'Measure', 'peakampbl',...
            'Neighborhood',  3, 'Peakpolarity', 'positive', 'Peakreplace', 'absolute', 'Resolution',  3, 'SendtoWorkspace', 'on' );
        precedingPositivePeak(set,chanCounter) = ERP_MEASURES;
        
        chanCounter = chanCounter + 1;
        
        ALLERP = pop_geterpvalues( ALLERP, [0 negativePeakLatency(set,chanCounter)], analysisBins(2), analysisChannels(chan) , 'Baseline', baseline, 'Binlabel', 'on', 'Erpsets',  erpSets(set),...
            'FileFormat', 'long', 'Filename', 'C:\Users\dbennett1\Desktop\x1.xls', 'Fracreplace', 'NaN', 'InterpFactor',  1, 'Measure', 'peakampbl',...
            'Neighborhood',  3, 'Peakpolarity', 'positive', 'Peakreplace', 'absolute', 'Resolution',  3, 'SendtoWorkspace', 'on' );
        precedingPositivePeak(set,chanCounter) = ERP_MEASURES;
    end
end

frnAmplitude = precedingPositivePeak - negativePeakAmplitude;




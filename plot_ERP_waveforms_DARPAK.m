

clear
close all

%%
%Load the grand average data:
ERP = pop_loaderp( 'filename', 'grandaverage_precsd.erp', 'filepath',...
    'C:\Users\RyanDataPC\Documents\Stefan Honours theses\DARPAK\data\erp\current\' );

% ERP.bindata is a channel (68) * time (1536) * bin (9) 3D array
% We're primarily interested in bins 1:6

% Plot P3a amplitudes for the 5 electrodes/channels we analysed:
% CPz (32), AFz (37), Fz (38), FCz (47) and Cz (48)

electrode = [32 37 38 47 48];
elc_labels = {'CPz (32)', 'AFz (37)', 'Fz (38)', 'FCz (47)', 'Cz (48)'};

%%
for chan = 1:length(electrode)
    
    % First, just overall waveforms for feedback conditions:
    figure
    subplot(1,2,1)
    plot(ERP.times, ERP.bindata(electrode(chan), :, 6), 'r-', 'LineWidth', 2) %plot rewarding condition
    hold on
    plot(ERP.times, ERP.bindata(electrode(chan), :, 5), 'k-', 'LineWidth', 2) %plot instructive condition
    xlim([-100 1000]) % Same x-limits as Bennett manuscript
    xlabel('Time from feedback onset (ms)')
    ylabel('Amplitude (uV)')
    legend({'Monetary feedback', 'Instructive feedback'})
    set(gca, 'TickDir', 'out', 'box', 'off', 'FontSize', 14, 'YGrid', 'on')
    
    % Now plot data split according to belief update size
    subplot(1,2,2)
    plot(ERP.times, ERP.bindata(electrode(chan), :, 4), 'r-', 'LineWidth', 2) %plot rewarding, large
    hold on
    %plot(ERP.times, ERP.bindata(electrode(chan), :, 3), 'r:', 'LineWidth', 3) %plot rewarding, small
    plot(ERP.times, ERP.bindata(electrode(chan), :, 3), 'Color', [1 0.5 0.5], 'LineStyle', ':', 'LineWidth', 3) %plot rewarding, small
    plot(ERP.times, ERP.bindata(electrode(chan), :, 2), 'k-', 'LineWidth', 2)  %plot instructive, large
    %plot(ERP.times, ERP.bindata(electrode(chan), :, 1), 'k:', 'LineWidth', 3) %plot instructive, large
    plot(ERP.times, ERP.bindata(electrode(chan), :, 1), 'Color', [0.5 0.5 0.5], 'LineStyle', ':', 'LineWidth', 3) %plot instructive, small
    xlim([-100 1000]) % Same x-limits as Bennett manuscript
    legend({'Monetary: large update', 'Monetary: small update', ...
        'Instructive: large update', 'Instructive: small update'})
    set(gca, 'TickDir', 'out', 'box', 'off', 'FontSize', 14, 'YGrid', 'on')
    xlabel('Time from feedback onset (ms)')
    
    suptitle (['P3a window, electrode ' elc_labels{chan}]);
    
end


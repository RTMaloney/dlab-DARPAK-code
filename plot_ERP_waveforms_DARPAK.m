
% Load and plot the grand-averaged ERP waveforms, both as simple line plots
% and as a scalp map.

clear
close all

% Load custom colormap: made using colormapeditor for orientation study (see Maloney & Clifford, 2015, Neuroimage).
% this is a modification of standard jet colormap with less saturated hot (red) end
% and is used later for the scalp plots
load ('C:\Users\RyanDataPC\Documents\Stefan Honours theses\DARPAK\dlab-DARPAK-code\OrnColMap.mat')

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

%% Just plot overall feedback conditions: site Fz (ch 38).
% Publication quality (for import to Illustrator)
% We use Fz because that is what Daniel chose for the manuscript;
% but as you can see in the above figures, all electrode sites show the effect pretty well

figure
% First, plot shaded bar showing where the P3a time window was defined. 
fill([250:550, 550:-1:250], ... % x values to fill along
    [ones(1,301)*15, ones(1,301)*-2], ... % y values to fill along: the y axis limits
    [0.85 0.85 0.85],'LineStyle', 'none')
hold on
plot(ERP.times, ERP.bindata(38, :, 6), 'r-', 'LineWidth', 2) %plot rewarding condition
plot(ERP.times, ERP.bindata(38, :, 5), 'k-', 'LineWidth', 2) %plot instructive condition
xlim([-100 1000]) % Same x-limits as Bennett manuscript
ylim([-2 15])

set(gca, 'TickDir', 'out', 'box', 'off', 'FontSize', 14) %, 'YGrid', 'on')


%% try to make scalp maps (P3a):

% We are looking at bins 5 (instructive), 6 (monetary) and 9 (monetary-instructive) only

% NOTE: we need to make map limits consistent for both feedback conditions
% for nice scalp plots.
% ** ALSO: need to make sure we don't include the external electrodes: 65:68 in the plot. 

latencies = [250 550]; % latencies to measure mean within (P3a analysis window)

datap = [];
% Pull out data for this bin:
for c=1:64
    datap(c,:) = ERP.bindata(c,:,9);
end

% Find the indices of the values in which we're interested in analysing:
latIndx = (ERP.times >= latencies(1) & ERP.times <= latencies(2)); 

% Just plot mean
data2plot = mean(datap(:,latIndx), 2);
% limit of color map:
maplimit = [0 max(data2plot)];
% full range:
maplimit = [min(data2plot) max(data2plot)];

maplimit = [0 12]; % for instructive and monetary plots
maplimit = [-2 2]; % for the difference between the two


figure
%topoplot( data2plot, ERP.chanlocs)
topoplot( data2plot, ERP.chanlocs,'maplimits', maplimit)
set(gcf, 'Colormap', OrnColMap) %set the custom colormap 
colorbar('FontSize', 14, 'TickDirection', 'out')

% other topoplot options (from pop_scalplot):
% ,...
%     'style', smapstyle, 'plotrad',mplotrad, 'headrad', mheadrad,'emarker', {'.','k',[],1},...
%     'numcontour', mapnumcontour, 'maplimits', maplimit, 'colormap', clrmap,'electrodes', elestyle, 'nosedir', mapview);


% ERP = pop_scalplot( ERP,  5, [ 250 550] , 'Blc', 'none', 'Colorbar', 'on', 'Colormap', 'hsv', 'Electrodes', 'on', 'FontName', 'Courier New',...
%  'FontSize',  10, 'Legend', 'bn-la', 'Maplimit', 'absmax', 'Mapstyle', 'both', 'Maptype', '2D', 'Mapview', '+X', 'Plotrad',  0.55, 'Position',...
%  [ 224.333 50.6667 791.333 514], 'Value', 'mean' );

%% Look at LPP window 

latencies = [550 900]; % latencies to measure mean within (LPP analysis window)

datap=[];
% Pull out data for this bin:
for c=1:64
    datap(c,:) = ERP.bindata(c,:,9);
end

% Find the indices of the values in which we're interested in analysing:
latIndx = (ERP.times >= latencies(1) & ERP.times <= latencies(2)); 

% Just plot mean
data2plot = mean(datap(:,latIndx), 2);
% limit of color map:
maplimit = [0 max(data2plot)];
% full range:
maplimit = [min(data2plot) max(data2plot)];

maplimit = [0 10] % map limit for mean instructive and monetary condtions
%maplimit = [-2 2]; % map limit for difference map

figure
%topoplot( data2plot, ERP.chanlocs)
topoplot( data2plot, ERP.chanlocs,'maplimits', maplimit)
set(gcf, 'Colormap', OrnColMap) %set the custom colormap 
colorbar('FontSize', 14, 'TickDirection', 'out')

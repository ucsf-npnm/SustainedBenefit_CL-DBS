close all
clear all
clc
%%
% Import Redcap data (from csv)

patientID = 'PR01';

[histogramHourly,redcapDataFile,ecogCatalog,studyVisitDates,stimChangeDates,...
    stimTurnOn,stimTurnOff,detectorChangeDates,chNames,enrollment,stage2date_start,...
    stage2date_end,stage3date_start,ltfudate_start,patientColor] = PresidioPatientData_PR01(patientID);

endDateForBiomarkerAnalysis = datetime('now');
crossover1Dates = datetime('2021-06-16'):caldays(1):datetime('2021-08-13');
crossover2Dates = datetime('2022-04-25'):caldays(1):datetime('2022-05-06');
excludeDates_crossoverPeriod = [crossover1Dates crossover2Dates];

savedWavelet = 'PR01_SelectData';

%%
% Channel names
chNames = {'VCVS1 - VCVS2','VCVS3 - VCVS4','Amyg1 - Amyg2','Amyg3 - Amyg4'};
numChans = length(chNames);

%%
% Wavelet
Fs = 250;
FArray = 1:125;
rangeCycles = [20 20]; % Larger number of cycles provides better frequency precision, but poorer temporal precision

load(savedWavelet)

% Apply endDateForBiomarkerAnalysis to loaded selectData
selectData = selectData(selectData.completion_pt_timestamp < endDateForBiomarkerAnalysis,:);

% Exclude crossover dates
clear indicesToRemove
for iDate = 1:size(selectData,1)
    toTest_dates = datetime(selectData.completion_pt_timestamp(iDate),'Format','dd-MMM-uuuu');
    selectDataDates_string(iDate) = string(toTest_dates);
end
indicesToRemove = find(ismember(selectDataDates_string,string(excludeDates_crossoverPeriod)));
selectData(indicesToRemove,:) = [];
numScores = size(selectData,1);

%%
% Average power within canonical frequency bands
delta = [1 4];
theta = [4 9] ;
alpha = [9 13];
beta = [13 30];
lowGamma = [31 70];
highGamma = [71 125];

allBands = {theta, alpha, beta, lowGamma, highGamma};
bandNames = {'Theta','Alpha','Beta','LowGamma','HighGamma'};

for iBand = 1:length(allBands)
    freqIndices{iBand} = find( (FArray > allBands{iBand}(1)) & (FArray <= allBands{iBand}(2)) );
end

% Add column to selectData with canonical freq band average power (iChan x iBand)
for iData = 1:size(selectData,1)
    if ~isempty(selectData.power{iData})
        currentData = selectData.power{iData};
        for iBand = 1:length(allBands)
            bandPower(:,iBand) = nanmean(currentData(:,freqIndices{iBand}),2);
        end
        selectData.bandPower(iData) = {bandPower};
    end
end

%%
% Gather all power data to z-score across recordings
clear power
for iData = 1:size(selectData,1)
    power(iData,:,:) = selectData.power{iData};
end

% Calculate z-scored power (across recordings, for each channel and frequency)
for iChan = 1:numChans
    zPower(:,iChan,:) = zscore(squeeze(power(:,iChan,:)),[],1);
end

%%
% Sliding window correlation of symptom vs power

select_VASD = selectData.vas_depression;
select_HAMD = selectData.hamd_total;

allMetrics = {select_VASD; select_HAMD};
metricNames = {'VAS-Depression','HAMD'};

allTimestamps = selectData.completion_pt_timestamp;
allTimestamps_relative = days(allTimestamps - enrollment);

firstDate_relative = allTimestamps_relative(1);

correlationType = 'Spearman'; % Spearman or Pearson
%%
% Calculate sliding window correlation of symptom vs band-averaged RAW power

calculationWindow = 180; % days
slidingWindow = 30; % days

startDate = firstDate_relative;
endDate = allTimestamps_relative(end);

clear allR_bands allp_bands slidingIndices
dateCounter = 1;
for iDate = floor(startDate):slidingWindow:floor(endDate)

    currentDateRange = [floor(startDate) floor(startDate + calculationWindow)];
    % Find indices corresponding to currentDateRange
    currentIndices = intersect(find(allTimestamps_relative > currentDateRange(1)),...
        find(allTimestamps_relative < currentDateRange(2)));
    slidingIndices{dateCounter} = currentIndices;

    for iMetric = 1:length(metricNames)
        disp(['Calculating correlations for ' metricNames{iMetric} ' (by frequency band)'])
        for iChan = 1:numChans
            for iFreq = 1:length(freqIndices)
                currentFreqs = freqIndices{iFreq};
                if ~isempty(currentIndices)
                    powerToTest = squeeze(nanmean(zPower(currentIndices,iChan,currentFreqs),3)); % Average z-power across freq indices for current band
                    [R,p] = corr((allMetrics{iMetric}(currentIndices)),powerToTest,'type',correlationType);
                    allR_bands(dateCounter,iMetric,iChan,iFreq) = R;
                    allp_bands(dateCounter,iMetric,iChan,iFreq) = p;
                else
                    allR_bands(dateCounter,iMetric,iChan,iFreq) = NaN;
                    allp_bands(dateCounter,iMetric,iChan,iFreq) = NaN;
                end
            end
        end
    end
    dateCounter = dateCounter + 1;

    % Iterate startDate by slidingWindow
    startDate = startDate + slidingWindow;
end

%%
% Plot sliding window correlation of symptom vs band-averaged RAW power
turboWithWhite = [[1 1 1]; turbo];

xTickLabels = firstDate_relative:slidingWindow:floor(endDate);
iMetric = 2; % 1 = VAS-D, 2 = HAMD-6
selectChan = 4; % Amyg 3 - Amyg 4

clear toPlotR toPlotP
figure
clf
set(gcf,'Position',[250.6000  615.4000  977.6000  178.4000])
toPlotR = squeeze(allR_bands(:,iMetric,selectChan,:))';
toPlotP = squeeze(allp_bands(:,iMetric,selectChan,:))';

% Create mask of above using p-values (only show if p < 0.05)
sig_pValues = ones(size(toPlotP));
nonSigIndices = find(toPlotP >=0.05);
sig_pValues(nonSigIndices) = NaN;

imagesc(toPlotR);
hold on

% External function to create outline of significant correlations
ClusterMap = ones(size(toPlotP));
ClusterMap(toPlotP >= 0.05) = 0;
[ClusterVertices] = ClusterBorder(ClusterMap);
Perimeters = ClusterVertices(:,3);
Nperimeters = max(Perimeters);
for iPerimeter = 1:Nperimeters()
    patch('Faces',1:sum(Perimeters == iPerimeter),'Vertices',...
        ClusterVertices(Perimeters == iPerimeter,1:2),...
        'FaceColor','none', 'EdgeColor','k','LineWidth',1)
end

set(gca,'Ydir','normal')
set(gca,'xtick',[])
yticks(1:length(bandNames))
yticklabels(bandNames)
currentYLim = get(gca,'Ylim');
line([days(firstDate_relative + stage2date_start - enrollment) days(firstDate_relative + stage2date_start - enrollment)],currentYLim,'Color','r')
line([days(firstDate_relative + ltfudate_start - enrollment) days(firstDate_relative + ltfudate_start - enrollment)],currentYLim,'Color','r')
a = colorbar;
ylabel(a,'Correlation coefficient (R)')

title([metricNames{iMetric} ': ' num2str(calculationWindow) ' days averaged, '...
    num2str(slidingWindow) ' day sliding window; RAW POWER'])

colormap(turboWithWhite)

%%
% Same as above, but with RELATIVE power
% Calcaulte sliding window correlation of symptom vs band-averaged RELATIVE power
% (Relative power calculation first, then z-score across recordings)

% Pull out average bandPower values
for iData = 1:numScores
    allBandPower(iData,:,:) = selectData.bandPower{iData};
end

% Relative power calculation
summedBandPower = nansum(allBandPower,3);
relativeBandPower = allBandPower ./ summedBandPower;

% Calculate z-score relative power across recordings (for each channel and band)
for iChan = 1:numChans
    zRelativeBandPower(:,iChan,:) = zscore(squeeze(relativeBandPower(:,iChan,:)),[],1);
end

calculationWindow = 180; % days
slidingWindow = 7; % days

startDate = firstDate_relative;
endDate = allTimestamps_relative(end);

clear allR_bands allp_bands slidingIndices
dateCounter = 1;
for iDate = floor(startDate):slidingWindow:floor(endDate)

    currentDateRange = [floor(startDate) floor(startDate + calculationWindow)];
    % Find indices corresponding to currentDateRange
    currentIndices = intersect(find(allTimestamps_relative > currentDateRange(1)),...
        find(allTimestamps_relative < currentDateRange(2)));
    slidingIndices{dateCounter} = currentIndices;

    for iMetric = 1:length(metricNames)
        disp(['Calculating correlations for ' metricNames{iMetric} ' (by frequency band)'])
        for iChan = 1:numChans
            for iFreq = 1:length(freqIndices)
                if ~isempty(currentIndices)
                    powerToTest = squeeze(zRelativeBandPower(currentIndices,iChan,iFreq));
                    [R,p] = corr((allMetrics{iMetric}(currentIndices)),powerToTest,'type',correlationType);
                    allR_bands(dateCounter,iMetric,iChan,iFreq) = R;
                    allp_bands(dateCounter,iMetric,iChan,iFreq) = p;
                else
                    allR_bands(dateCounter,iMetric,iChan,iFreq) = NaN;
                    allp_bands(dateCounter,iMetric,iChan,iFreq) = NaN;
                end
            end
        end
    end
    dateCounter = dateCounter + 1;

    % Iterate startDate by slidingWindow
    startDate = startDate + slidingWindow;
end


%%
% Plot sliding window correlation of symptom vs band-averaged RELATIVE power

xTickLabels = firstDate_relative:slidingWindow:floor(endDate);
iMetric = 2; % 1 = VAS-D, 2 = HAMD-6
selectChan = 4; % Amyg 3 - Amyg 4

clear toPlotR toPlotP
figure
clf
set(gcf,'Position',[250.6000  615.4000  977.6000  178.4000])
toPlotR = squeeze(allR_bands(:,iMetric,selectChan,:))';
toPlotP = squeeze(allp_bands(:,iMetric,selectChan,:))';

% Create mask of above using p-values (only show if p < 0.05)
sig_pValues = ones(size(toPlotP));
nonSigIndices = find(toPlotP >=0.05);
sig_pValues(nonSigIndices) = NaN;

imagesc(toPlotR);
hold on

% External function to create outline of significant correlations
ClusterMap = ones(size(toPlotP));
ClusterMap(toPlotP >= 0.05) = 0;
[ClusterVertices] = ClusterBorder(ClusterMap);
Perimeters = ClusterVertices(:,3);
Nperimeters = max(Perimeters);
for iPerimeter = 1:Nperimeters()
    patch('Faces',1:sum(Perimeters == iPerimeter),'Vertices',...
        ClusterVertices(Perimeters == iPerimeter,1:2),...
        'FaceColor','none', 'EdgeColor','k','LineWidth',1)
end
set(gca,'Ydir','normal')
set(gca,'xtick',[])
yticks(1:length(bandNames))
yticklabels(bandNames)
currentYLim = get(gca,'Ylim');
line([days(firstDate_relative + stage2date_start - enrollment) days(firstDate_relative + stage2date_start - enrollment)],currentYLim,'Color','r')
line([days(firstDate_relative + ltfudate_start - enrollment) days(firstDate_relative + ltfudate_start - enrollment)],currentYLim,'Color','r')
a = colorbar;
ylabel(a,'Correlation coefficient (R)')

title([metricNames{iMetric} ': ' num2str(calculationWindow) ' days averaged, '...
    num2str(slidingWindow) ' day sliding window; RELATIVE POWER'])
colormap(turboWithWhite)

%%
% Averaging within trial periods

stage2date_start_relative = days(stage2date_start - enrollment);
stage2date_end_relative = days(stage2date_end - enrollment);
stage3date_start_relative = days(stage3date_start - enrollment);
ltfudate_start_relative = days(ltfudate_start - enrollment);

dateRangeStart = stage2date_start_relative;
dateRangeEnd = stage2date_end_relative;

% dateRangeStart = stage2date_end_relative;
% dateRangeEnd = ltfudate_start_relative;

% dateRangeStart = ltfudate_start_relative;
% dateRangeEnd = days(datetime('now') - enrollment);

% dateRangeStart = stage2date_start_relative;
% dateRangeEnd = days(datetime('now') - enrollment);


startDate = max(firstDate_relative,dateRangeStart);
endDate = min(allTimestamps_relative(end),dateRangeEnd);

clear allR_bands allp_bands
currentDateRange = [floor(startDate) endDate];
% Find indices corresponding to currentDateRange
currentIndices = intersect(find(allTimestamps_relative > currentDateRange(1)),...
    find(allTimestamps_relative < currentDateRange(2)));

for iMetric = 1:length(metricNames)
    disp(['Calculating correlations for ' metricNames{iMetric} ' (by frequency band)'])
    for iChan = 1:numChans
        for iFreq = 1:length(freqIndices)
            if ~isempty(currentIndices)
                powerToTest = squeeze(zRelativeBandPower(currentIndices,iChan,iFreq));
                [R,p] = corr((allMetrics{iMetric}(currentIndices)),powerToTest,'type',correlationType);
                allR_bands(iMetric,iChan,iFreq) = R;
                allp_bands(iMetric,iChan,iFreq) = p;
            else
                allR_bands(iMetric,iChan,iFreq) = NaN;
                allp_bands(iMetric,iChan,iFreq) = NaN;
            end
        end
    end
end

iMetric = 2;
iChan = 4;

% Create mask of above using p-values (only show if p < 0.05)
toPlotR = squeeze(allR_bands(iMetric,iChan,:));
toPlotP = squeeze(allp_bands(iMetric,iChan,:));

sig_pValues = ones(size(toPlotP));
nonSigIndices = toPlotP >= 0.05;
sig_pValues(nonSigIndices) = NaN;

mask = repmat(isnan((sig_pValues)),1,1,3);

figure
clf
set(gcf,'Position',[ 208.2000  329.8000  204.8000  385.6000])
imagesc((toPlotR),[-1 1]);
set(gca,'ydir','normal')
hold on
h = imagesc(mask);
h.AlphaData = isnan(sig_pValues);
a = colorbar;
ylabel(a,'Correlation coefficient (R)')
yticks(1:length(bandNames))
yticklabels((bandNames))
colormap turbo

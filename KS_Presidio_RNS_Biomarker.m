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
if strcmp(patientID,'PR01')
    crossover1Dates = datetime('2021-06-16'):caldays(1):datetime('2021-08-13');
    crossover2Dates = datetime('2022-04-25'):caldays(1):datetime('2022-05-06');
    excludeDates_crossoverPeriod = [crossover1Dates crossover2Dates];
    daysToExclude = 60; % Number of days of neural data to exclude following Stage 2 implant
end
savedWavelet = 'PR01_SelectData.mat';
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

% If applicable, remove days at the beginning of Stage 2
selectData = selectData(selectData.completion_pt_timestamp > stage2date_start + daysToExclude,:);

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
        allBandPower(iData,:,:) = bandPower;
    end
end

% Zscore allBandPower
for iChan = 1:numChans
    zAllBandPower(:,iChan,:) = zscore(squeeze(allBandPower(:,iChan,:)),[],1);
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
                    powerToTest = squeeze(zAllBandPower(currentIndices,iChan,iFreq));
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

% Using date in the middle of calculation window for x-axis
xTickLabels = firstDate_relative + calculationWindow/2:slidingWindow:floor(endDate) + calculationWindow/2;

iMetric = 2; % 1 = VAS-D, 2 = HAMD-6
selectChan = 4; % 4 = Amyg 3 - Amyg 4

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

imagesc(toPlotR,[-1 1]);
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
a = colorbar;
ylabel(a,'Correlation coefficient (R)')

title([metricNames{iMetric} ': ' num2str(calculationWindow) ' days averaged, '...
    num2str(slidingWindow) ' day sliding window; RAW POWER'])

colormap(turboWithWhite)


%%
% Shuffle HAMD-6SR data 
% Calculate sliding window correlation of symptom vs band-averaged RAW power

verum_toPlotR = toPlotR;
verum_toPlotP = toPlotP;

iMetric = 2; % 1 = VAS-D, 2 = HAMD-6
selectChan = 4; % 4 = Amyg 3 - Amyg 4

calculationWindow = 180; % days
slidingWindow = 7; % days

numShuffle = 10000;

startDate = firstDate_relative;
endDate = allTimestamps_relative(end);
numDates = length(floor(startDate):slidingWindow:floor(endDate));

shuffledR_bands = NaN(numShuffle,numDates,length(freqIndices));

for iShuffle = 1:numShuffle
    disp(['Shuffle: ' num2str(iShuffle) ' of ' num2str(numShuffle)])
    % Shuffle scores
    currentMetric = allMetrics{iMetric};
    shuffledMetric = shuffle_preserve_acorr(currentMetric);

    startDate = firstDate_relative;
    endDate = allTimestamps_relative(end);

    dateCounter = 1;
    for iDate = floor(startDate):slidingWindow:floor(endDate)
        currentDateRange = [floor(startDate) floor(startDate + calculationWindow)];
        % Find indices corresponding to currentDateRange
        currentIndices = intersect(find(allTimestamps_relative > currentDateRange(1)),...
            find(allTimestamps_relative < currentDateRange(2)));
    
        iChan = selectChan;
        for iFreq = 1:length(freqIndices)
            currentFreqs = freqIndices{iFreq};
            if ~isempty(currentIndices)
                clear powerToTest R p
                powerToTest = squeeze(zAllBandPower(currentIndices,iChan,iFreq));
                [R,p] = corr((shuffledMetric(currentIndices)),powerToTest,'type',correlationType);
                shuffledR_bands(iShuffle,dateCounter,iFreq) = R;
            end
        end
        dateCounter = dateCounter + 1;



        % Iterate startDate by slidingWindow
        startDate = startDate + slidingWindow;
    end
end

permutationPvalues = NaN(length(freqIndices),numDates);

for iDate = 1:numDates
    for iFreq = 1:length(freqIndices)
        currentVerumR_bands = verum_toPlotR(iFreq,iDate);
        currentShuffledR_bands = shuffledR_bands(:,iDate,iFreq);
        permutationPvalues(iFreq,iDate) = (sum(abs(currentShuffledR_bands) >= abs(currentVerumR_bands))) / (numShuffle); 
    end
end

selectSigLevel = 0.05;

% Non-significant pValues == NaN
permutation_sig_pValues = permutationPvalues;
permutation_nonSigIndices = find(permutationPvalues >= selectSigLevel);
permutation_sig_pValues(permutation_nonSigIndices) = NaN;

figure
clf
set(gcf,'Position',[250.6000  615.4000  977.6000  178.4000])
imagesc((permutation_sig_pValues),[0 selectSigLevel])
hold on
set(gca,'Ydir','normal')
set(gca,'xtick',[])
yticks(1:length(bandNames))
yticklabels(bandNames)
currentYLim = get(gca,'Ylim');
a = colorbar;
ylabel(a,'P-value')

copperWithWhite = [[1 1 1]; copper];
colormap(copperWithWhite)


figure
clf
set(gcf,'Position',[250.6000  615.4000  977.6000  178.4000])
imagesc(verum_toPlotR,[-1 1])
% imagesc(xTickLabels,[],verum_toPlotR,[-1 1])
hold on

% External function to create outline of significant correlations
ClusterMap = ones(size(permutationPvalues));
ClusterMap(permutationPvalues >= selectSigLevel) = 0;
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
a = colorbar;
ylabel(a,'Correlation coefficient (R)')
title([metricNames{iMetric} ': ' num2str(calculationWindow) ' days averaged, '...
    num2str(slidingWindow) ' day sliding window; RAW POWER'])

turboWithWhite = [[1 1 1]; turbo];
colormap(turboWithWhite)

%%
% Quantify score ranges for each significant biomarker period

startDate = firstDate_relative;
endDate = allTimestamps_relative(end);
allDates = floor(startDate):slidingWindow:floor(endDate);

clear symptomsToTest
% For each timepoint in the imagesc, slidingIndices has the indices (referring back to selectData) 
for iFreq = 1:length(freqIndices)
    sigIndices = find(permutation_sig_pValues(iFreq,:) > 0);
    currentIndices = unique(vertcat(slidingIndices{sigIndices}));
    symptomsToTest{iFreq} = allMetrics{iMetric}(currentIndices);
end

% Create grouping vector to use for boxplot
scores = vertcat(symptomsToTest{:});
n = cellfun(@numel,symptomsToTest);
groupID = repelem(1:numel(symptomsToTest),n);

figure
clf
set(gcf,'Position',[1     1   246   423]);
boxplot(scores,groupID,'MedianStyle','line','Color',[0 0.4470 0.7410],'Widths',0.5,...
    'Labels',bandNames)
ylabel(metricNames{iMetric})

[p,table,stats] = kruskalwallis(scores,groupID);
[results,~,~,gNames] = multcompare(stats);

%%
% Save shuffle output
% save('PR01_ShuffledRNSBiomarker_HAMD.mat')

%%
% Averaging within trial periods

stage2date_start_relative = days(stage2date_start - enrollment);
stage2date_end_relative = days(stage2date_end - enrollment);
stage3date_start_relative = days(stage3date_start - enrollment);
ltfudate_start_relative = days(ltfudate_start - enrollment);

dateRangeStart = stage2date_start_relative;
dateRangeEnd = days(datetime('now') - enrollment);

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
                powerToTest = squeeze(zAllBandPower(currentIndices,iChan,iFreq));
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
iBand = 4;

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

% Scatter plot of zScored band neural power vs symptom score
figure
clf
set(gcf,'Position',[568   370   509   289])
powerToPlot = squeeze(zAllBandPower(:,iChan,iBand));
scoresToPlot = allMetrics{iMetric};
scatter(powerToPlot,scoresToPlot,30,'filled','k')
h = lsline;
h.Color = 'r';
xlabel(['Amyg3-Amyg4 ' bandNames{iBand} ' Power'])
ylabel(metricNames{iMetric})
axis square
currentYLim = get(gca,'ylim');
ylim([0 currentYLim(2)])
title(['R = ' num2str(allR_bands(iMetric,iChan,iBand)) '; p = ' num2str(allp_bands(iMetric,iChan,iBand))])

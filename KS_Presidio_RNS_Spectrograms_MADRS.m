clear all
close all
clc
%%
patientID = 'PR01';
[histogramHourly,redcapDataFile,ecogCatalog,studyVisitDates,stimChangeDates,...
    stimTurnOn,stimTurnOff,detectorChangeDates,chNames,enrollment,stage2date_start,...
    stage2date_end,stage3date_start,ltfudate_start,patientColor] = PresidioPatientData(patientID);

load('PR01_ComprehensiveMADRS_07032025.mat')

clinicianScales = PR01_MADRS;

stimChangeDates_relative = days(stimChangeDates - enrollment);
stimTurnOn_relative = days(stimTurnOn - enrollment);
stimTurnOff_relative = days(stimTurnOff - enrollment);

load([pwd filesep patientID '_SelectMagnetData_07162025.mat'])
longitudinalSpectrogramFile = [patientID '_LongitudinalSpectrograms_07162025.mat'];

startAnalysisDate = enrollment;
endAnalysisDate = datetime('now');

crossover1Dates = []; % Removed
crossover2Dates = []; % Removed
excludeDates_crossoverPeriod = [crossover1Dates crossover2Dates];

%%
% Load Redcap data
redcapData = readtable(redcapDataFile);

%%
% Remove visit days from redcap data and plot scores
% Indices corresponding to start and end dates for analysis
indices_surveys = intersect(find(redcapData.completion_pt_timestamp >= startAnalysisDate),...
    find(redcapData.completion_pt_timestamp <= endAnalysisDate));

% Indices of study visit dates
redcap_allDates = datetime(redcapData.completion_pt_timestamp,...
    'Format','dd-MMM-uuuu');
redcapAllDates_string = string(redcap_allDates);
indicesToKeep = find(~ismember(redcapAllDates_string,string(studyVisitDates)));

% Indices to plot
indicesToPlot = intersect(indices_surveys,indicesToKeep);

selectTime = redcapData.completion_pt_timestamp(indicesToPlot);
selectData = redcapData(indicesToPlot,:);

%%
% Wavelet
Fs = 250;
FArray = 1:125;
rangeCycles = [20 20]; % Larger number of cycles provides better frequency precision, but poorer temporal precision

numECOG = size(catalog_selectMagnetData,1);
load(longitudinalSpectrogramFile)

%%
% All magnet recordings = Plot spectrograms: Across all recordings and then smoothed based on days -- one sample per day

% Vector of days
elapsedDays_allMagnetRecordings = days(catalog_selectMagnetData.Timestamp - startAnalysisDate);

allDays_allRecordings = floor(elapsedDays_allMagnetRecordings(1)):floor(elapsedDays_allMagnetRecordings(end));
firstDay_allRecordings = floor(elapsedDays_allMagnetRecordings(1));

% Temporal smoothing of spectrograms (mean across previous numDaysToSmooth)
numDaysToSmooth = 1; % Using the past x days for smoothing
zLongitudinalSpectrum_byDays = NaN(4,length(allDays_allRecordings),length(FArray));
longitudinalSpectrum_byDays = NaN(4,length(allDays_allRecordings),length(FArray));
for iDay = 1:length(allDays_allRecordings)
    clear indicesForAverage
    currentDay = allDays_allRecordings(iDay);
    % Only if there is enough historical data for smoothing
    if iDay >= numDaysToSmooth
        % Find samples which correspond to the past numDaysToSmooth
        indicesForAverage = find(elapsedDays_allMagnetRecordings >= currentDay-numDaysToSmooth & elapsedDays_allMagnetRecordings <= currentDay);
        zLongitudinalSpectrum_byDays(:,iDay,:) = squeeze(mean(zLongitudinalSpectrum(:,indicesForAverage,:),2));
        longitudinalSpectrum_byDays(:,iDay,:) = squeeze(mean(longitudinalSpectrum(indicesForAverage,:,:),1)); % Order of variables for longitudinalSpectrum is different than zLongitudinalSpectrum
    end
end

% Create mask to make NaN values white
maskValues = isnan(zLongitudinalSpectrum_byDays);

figure
clf
set(gcf,'Position',[217.8000  460.2000  968.8000  330.4000])
iChan = 4;
subplot(2,1,1)
imagesc([], FArray,squeeze(zLongitudinalSpectrum(iChan,:,:))',[-2 2]);
set(gca,'YDir','normal')
title(chNames{iChan})
xlabel('Recording Number')
ylabel('Frequency [Hz]')
colorbar

subplot(2,1,2)
imagesc(allDays_allRecordings, FArray,squeeze(zLongitudinalSpectrum_byDays(iChan,:,:))',[-2 2]);
hold on
currentMask = squeeze(maskValues(iChan,:,:))';
toPlotMask = repmat(currentMask,1,1,3);
h = imagesc(allDays_allRecordings,[],toPlotMask);
h.AlphaData = currentMask;
set(gca,'YDir','normal')
title([chNames{iChan} ' Smoothed (' num2str(numDaysToSmooth) ' day)'])
xlabel('Elapsed Days')
ylabel('Frequency [Hz]')
colorbar

sgtitle('All Magnet Recordings')
%%
% Plot indication of what elapsed day the magnet recordings come from
figure
clf
set(gcf,'Position',[488 438 937 85.4000])
stem(elapsedDays_allMagnetRecordings,ones(length(elapsedDays_allMagnetRecordings),1),'Marker','|','Color','black')
xlim([elapsedDays_allMagnetRecordings(1) elapsedDays_allMagnetRecordings(end)])
xlabel('Days of Trial Participation')

%%
% Band power / relative power for all magnet recordings

% Average power within canonical frequency bands
delta = [1 4];
theta = [4 9] ;
alpha = [9 13];
beta = [13 30];
lowGamma = [31 70];
highGamma = [71 125];
allGamma = [31 125];

allBands = {theta, alpha, beta, lowGamma, highGamma};
bandNames = {'Theta','Alpha','Beta','LowGamma','HighGamma'};

for iBand = 1:length(allBands)
    freqIndices{iBand} = find( (FArray > allBands{iBand}(1)) & (FArray <= allBands{iBand}(2)) );
end

% Band power calculation for longitudinalSpectrum_byDays (all recordings
% averaged per day)
clear bandPower
for iData = 1:size(longitudinalSpectrum_byDays,2)
    clear currentData
    currentData = squeeze(longitudinalSpectrum_byDays(:,iData,:));
    for iBand = 1:length(allBands)
        bandPower(:,iBand) = nanmean(currentData(:,freqIndices{iBand}),2);
    end
    longitudinalSpectrum_byDays_bandPower(iData,:,:) = bandPower;
end

% Relative power calculation
clear summedBandPower longitudinalSpectrum_byDays_bandPower_relativeBandPower
summedBandPower = nansum(longitudinalSpectrum_byDays_bandPower,3);
longitudinalSpectrum_byDays_bandPower_relativeBandPower = longitudinalSpectrum_byDays_bandPower ./ summedBandPower;
%%
% Plot high gamma power and relative high gamma power across days of trial
% participation (excluding first 60 days after implant)
daysToExclude = 60;
selectChan = 4; % Amyg 3- Amyg 4
selectBand = 5; % High gamma

figure
clf
set(gcf,'Position',[253 431.4000 935.2000 420.0000])
subplot(2,1,1)
toPlot = longitudinalSpectrum_byDays_bandPower(:,selectChan,selectBand);
scatter(allDays_allRecordings(daysToExclude:end),toPlot(daysToExclude:end),10,'filled','k')
h = lsline;
h.Color = 'r';
xlabel('Days of Trial Participation')
title('Amyg3-Amyg4: High Gamma Power')
ylabel('High Gamma Power')

subplot(2,1,2)
toPlot = longitudinalSpectrum_byDays_bandPower_relativeBandPower(:,selectChan,selectBand);
scatter(allDays_allRecordings(daysToExclude:end),toPlot(daysToExclude:end),10,'filled','k')
h = lsline;
h.Color = 'r';
xlabel('Days of Trial Participation')
ylabel('Relative High Gamma Power')
title('Amyg3-Amyg4: Relative High Gamma Power')

%%
% MADRS indices
MADRS_timestamps = clinicianScales.Date;

for iDate = 1:length(MADRS_timestamps)
    toTest_MADRSdates = datetime(MADRS_timestamps(iDate),'Format','dd-MMM-uuuu');
    allMADRSDates_string(iDate) = string(toTest_MADRSdates);
end

indicesToRemove = find(ismember(allMADRSDates_string,string(excludeDates_crossoverPeriod)));
indicesToKeep = setdiff(1:size(MADRS_timestamps,1),indicesToRemove);
MADRS_timestamps = MADRS_timestamps(indicesToKeep,:);
MADRS_total = clinicianScales.MADRS_Total(indicesToKeep,:);

MADRS_relativeTimestamps = days(MADRS_timestamps - enrollment);

%%
% Pull out segment of spectral activity around each MADRS score and look at
% correlation

% Subselect MADRS scores from Stage 2 and beyond (when there is
% corresponding neural data)
selectIndices = find(MADRS_relativeTimestamps > days(stage2date_start - enrollment));
selectMADRS_relativeTimestamps = MADRS_relativeTimestamps(selectIndices);
selectMADRS_total = MADRS_total(selectIndices);

%%
% Longitudinal biomarker (MADRS vs relative power in across full FArray) 

stage2date_start_relative = days(stage2date_start - enrollment);
stage2date_end_relative = days(stage2date_end - enrollment);
stage3date_start_relative = days(stage3date_start - enrollment);
ltfudate_start_relative = days(ltfudate_start - enrollment);

selectDaysToAverage = [1 7 14 21 28];

dateRangeStart = stage2date_start_relative;
dateRangeEnd = days(datetime('now') - enrollment);

clear neuralDataToPlot
neuralDataToPlot = longitudinalSpectrum_byDays;

selectChan = 4;
clear pValues correlationValues
for iDays = 1:length(selectDaysToAverage)
    currentDaysToAverage = selectDaysToAverage(iDays);

    % Pulling out data with paired MADRS and spectral data
    counter = 1;
    indicesToTest = [];
    clear selectSpectrum selectMADRS_total_withSpectrum
    for iScore = 1:length(selectMADRS_total)
        scoreDate = selectMADRS_relativeTimestamps(iScore);
        if (scoreDate > (firstDay_allRecordings + currentDaysToAverage)) && (scoreDate >= dateRangeStart)
            if (floor(scoreDate)-firstDay_allRecordings) < size(neuralDataToPlot,2) && (scoreDate <= dateRangeEnd)
                clear currentCalc
                currentCalc = squeeze(nanmean(neuralDataToPlot(selectChan,floor(scoreDate)-firstDay_allRecordings-currentDaysToAverage...
                    :floor(scoreDate)-firstDay_allRecordings,:),2));
                if sum(sum(~isnan(currentCalc)))
                    selectSpectrum(counter,:) = currentCalc;
                    counter = counter + 1;
                    indicesToTest = [indicesToTest iScore];
                end
            end
        end
    end
    selectMADRS_total_withSpectrum = MADRS_total(indicesToTest);
    length(indicesToTest)

    % Keeping full FArray without binning
    for iFreq = 1:length(FArray)
        clear currentFreqs dataToTest toRemove
        dataToTest = squeeze(selectSpectrum(:,iFreq));
        [R,P] = corrcoef(selectMADRS_total_withSpectrum,dataToTest);
        pValues(iDays,iFreq) = P(2);
        correlationValues(iDays,iFreq) = R(2);
    end
end

% External function to create outline of significant correlations
ClusterMap = ones(size(pValues));
ClusterMap(pValues >= 0.05) = 0;
[ClusterVertices] = ClusterBorder(ClusterMap);

figure
clf
set(gcf,'Position',[280.2000 597 876 201.6000])
imagesc(correlationValues,[-0.5 0.5]);
hold on
Perimeters = ClusterVertices(:,3);
Nperimeters = max(Perimeters);
for iPerimeter = 1:Nperimeters()
    patch('Faces',1:sum(Perimeters == iPerimeter),'Vertices',...
        ClusterVertices(Perimeters == iPerimeter,1:2),...
        'FaceColor','none', 'EdgeColor','k','LineWidth',1)
end
a = colorbar;
ylabel(a,'Correlation coefficient (R)')
xlabel('Frequency [Hz]')
ylabel('Number of days to average for power')
yticks(1:length(selectDaysToAverage))
yticklabels(selectDaysToAverage)
title('Correlation between MADRS and spectral power')
colormap jet

%%
% Longitudinal biomarker (MADRS vs relative power averaged in canonical frequency bands) - averaged per trial period

selectDaysToAverage = [1 7 14 21 28];

% dateRangeStart = stage2date_start_relative;
% dateRangeEnd = stage2date_end_relative;

% dateRangeStart = stage2date_end_relative;
% dateRangeEnd = ltfudate_start_relative;

% dateRangeStart = ltfudate_start_relative;
% dateRangeEnd = days(datetime('now') - enrollment);

dateRangeStart = stage2date_start_relative;
dateRangeEnd = days(datetime('now') - enrollment);

clear neuralDataToPlot
neuralDataToPlot = longitudinalSpectrum_byDays_bandPower_relativeBandPower;

selectChan = 4;
clear pValues correlationValues
for iDays = 1:length(selectDaysToAverage)
    currentDaysToAverage = selectDaysToAverage(iDays);

    % Pulling out data with paired MADRS and spectral data
    counter = 1;
    indicesToTest = [];
    clear selectSpectrum selectMADRS_total_withSpectrum
    for iScore = 1:length(selectMADRS_total)
        scoreDate = selectMADRS_relativeTimestamps(iScore);
        if (scoreDate > (firstDay_allRecordings + currentDaysToAverage)) && (scoreDate >= dateRangeStart)
            if (floor(scoreDate)-firstDay_allRecordings) < size(neuralDataToPlot,1) && (scoreDate <= dateRangeEnd)
                clear currentCalc
                currentCalc = squeeze(nanmean(neuralDataToPlot(floor(scoreDate)-firstDay_allRecordings-currentDaysToAverage...
                    :floor(scoreDate)-firstDay_allRecordings,selectChan,:),1));
                if sum(sum(~isnan(currentCalc)))
                    selectSpectrum(counter,:) = currentCalc;
                    counter = counter + 1;
                    indicesToTest = [indicesToTest iScore];
                end
            end
        end
    end
    selectMADRS_total_withSpectrum = MADRS_total(indicesToTest);
    length(indicesToTest)

    for iBand = 1:length(allBands)
        clear currentFreqs dataToTest toRemove
        dataToTest = squeeze(selectSpectrum(:,iBand));
        [R,P] = corrcoef(selectMADRS_total_withSpectrum,dataToTest);
        pValues(iDays,iBand) = P(2);
        correlationValues(iDays,iBand) = R(2);
    end

    % Scatter of MADRS vs relative band power
    selectBand = 5;

    figure
    clf
    scatter(selectMADRS_total_withSpectrum,selectSpectrum(:,selectBand),'filled','k')
    hold on
    h = lsline;
    h.Color = 'r';
    title({[chNames{selectChan} ': ' num2str(currentDaysToAverage) ' days neural activity'];...
        ['R = ' num2str(R(2)) ' p = ' num2str(P(2))]})
    xlabel('MADRS')
    ylabel(['Relative Power: ' bandNames{selectBand}])
end

% Create mask of above using p-values (only show if p < 0.05)
sig_pValues = ones(size(pValues));
nonSigIndices = pValues >= 0.05;
sig_pValues(nonSigIndices) = NaN;

mask = repmat(isnan(sig_pValues),1,1,3);

figure
clf
set(gcf,'Position',[208.2000  513.8000  381.6000  201.6000])
imagesc(correlationValues,[-0.5 0.5]);
hold on
h = imagesc(mask);
h.AlphaData = isnan(sig_pValues);
a = colorbar;
ylabel(a,'Correlation coefficient (R)')
xticks(1:length(bandNames))
xticklabels(bandNames)
ylabel('Number of days to average for power')
yticks(1:length(selectDaysToAverage))
yticklabels(selectDaysToAverage)
title('Correlation between MADRS and spectral power')
colormap jet
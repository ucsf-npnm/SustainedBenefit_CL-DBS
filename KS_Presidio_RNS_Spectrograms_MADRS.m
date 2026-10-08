clear all
close all
clc
%%
patientID = 'PR01';
[histogramHourly,redcapDataFile,ecogCatalog,studyVisitDates,stimChangeDates,...
    stimTurnOn,stimTurnOff,detectorChangeDates,chNames,enrollment,stage2date_start,...
    stage2date_end,stage3date_start,ltfudate_start,patientColor] = PresidioPatientData_PR01(patientID);

load('PR01_ComprehensiveMADRS.mat')
clinicianScales = PR01_MADRS;

stimChangeDates_relative = days(stimChangeDates - enrollment);
stimTurnOn_relative = days(stimTurnOn - enrollment);
stimTurnOff_relative = days(stimTurnOff - enrollment);

load('PR01_SelectMagnetData.mat')
longitudinalSpectrogramFile = 'PR01_LongitudinalSpectrograms.mat';

startAnalysisDate = enrollment;

if strcmp(patientID,'PR01')
    crossover1Dates = datetime('2021-06-16'):caldays(1):datetime('2021-08-13');
    crossover2Dates = datetime('2022-04-25'):caldays(1):datetime('2022-05-06');
    excludeDates_crossoverPeriod = [crossover1Dates crossover2Dates];
    daysToExclude = 60; % Number of days of neural data to exclude following Stage 2 implant
end
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
% participation (excluding first x days after implant)
selectChan = 4; % Amyg 3- Amyg 4
selectBand = 5; % High gamma

figure
clf
set(gcf,'Position',[253 431.4000 935.2000 420.0000])
subplot(2,1,1)
toPlot = longitudinalSpectrum_byDays_bandPower(:,selectChan,selectBand);
scatter(allDays_allRecordings(1+daysToExclude:end),toPlot(1+daysToExclude:end),10,'filled','k')
h = lsline;
h.Color = 'r';
xlabel('Days of Trial Participation')
title('Amyg3-Amyg4: High Gamma Power')
ylabel('High Gamma Power')

subplot(2,1,2)
toPlot = longitudinalSpectrum_byDays_bandPower_relativeBandPower(:,selectChan,selectBand);
scatter(allDays_allRecordings(1+daysToExclude:end),toPlot(1+daysToExclude:end),10,'filled','k')
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

% Starting to look at biomarker 60 days after implant (to
% allow for electrode stabilization)

selectIndices = find(MADRS_relativeTimestamps > days(stage2date_start + daysToExclude - enrollment));
selectMADRS_relativeTimestamps = MADRS_relativeTimestamps(selectIndices);
selectMADRS_total = MADRS_total(selectIndices);

%%
% Longitudinal biomarker (MADRS vs relative power in across full FArray) 

stage2date_start_relative = days(stage2date_start - enrollment);
stage2date_end_relative = days(stage2date_end - enrollment);
stage3date_start_relative = days(stage3date_start - enrollment);
ltfudate_start_relative = days(ltfudate_start - enrollment);

%%
% Scatter plot of correlation between MADRS and previous 7-days zScored band neural power
currentDaysToAverage = 7;

dateRangeStart = stage2date_start_relative;
dateRangeEnd = days(datetime('now') - enrollment);

% Zscore bandPower across recordings (e.g. across days)
clear neuralDataToPlot
for iChan = 1:4
    neuralDataToPlot(iChan,:,:) = longitudinalSpectrum_byDays_bandPower(:,iChan,:);

    % Calculate zScore
    clear temp;
    temp = squeeze(longitudinalSpectrum_byDays_bandPower(:,iChan,:));
    neuralDataToPlot_zScore(iChan,:,:) = normalize(temp,1);

end

iChan = 4;
iBand = 4;

% Pulling out data with paired MADRS and spectral data
counter = 1;
indicesToTest = [];
clear selectSpectrum selectMADRS_total_withSpectrum

for iScore = 1:length(selectMADRS_total)
    scoreDate = selectMADRS_relativeTimestamps(iScore);
    if (scoreDate > (firstDay_allRecordings + currentDaysToAverage)) && (scoreDate >= dateRangeStart)
        if (floor(scoreDate)-firstDay_allRecordings) < size(neuralDataToPlot_zScore,2) && (scoreDate <= dateRangeEnd)
            clear currentCalc
            currentCalc = squeeze(nanmean(neuralDataToPlot_zScore(selectChan,floor(scoreDate)-firstDay_allRecordings-currentDaysToAverage...
                :floor(scoreDate)-firstDay_allRecordings-1,:),2));
            if sum(sum(~isnan(currentCalc)))
                selectSpectrum(counter,:) = currentCalc;
                counter = counter + 1;
                indicesToTest = [indicesToTest iScore];
            end
        end
    end
end
selectMADRS_total_withSpectrum = selectMADRS_total(indicesToTest);

% Calculate correlation
[R,p] = corr(selectSpectrum(:,iBand),selectMADRS_total_withSpectrum);

figure
clf
set(gcf,'Position',[568   370   509   289])
scatter(selectSpectrum(:,iBand),selectMADRS_total_withSpectrum,30,'filled','k')
h = lsline;
h.Color = 'r';
xlabel(['Amyg3-Amyg4 ' bandNames{iBand} ' Power'])
ylabel('MADRS')
axis square
currentYLim = get(gca,'ylim');
ylim([0 currentYLim(2)])
title(['R = ' num2str(R) '; p = ' num2str(p)])




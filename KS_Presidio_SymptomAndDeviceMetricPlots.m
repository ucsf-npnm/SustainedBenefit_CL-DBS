clear all
close all
clc
%%
patientID = 'PR01';

[histogramHourly,redcapDataFile,ecogCatalog,studyVisitDates,stimChangeDates,...
    stimTurnOn,stimTurnOff,detectorChangeDates,chNames,enrollment,stage2date_start,...
    stage2date_end,stage3date_start,ltfudate_start,patientColor] = PresidioPatientData(patientID);

load('PR01_ComprehensiveMADRS_07032025.mat')
MADRS = PR01_MADRS;

redcapFile = 'PR01_AmbulatoryRedcap_07032025.csv';
histogramHourly = 'Histogram_Hourly_07032025.csv';

resetTime = 6; % In PST [should start getting stim within 7am bin]

crossover1Dates = []; % Removed
crossover2Dates = []; % Removed
excludeDates_crossoverPeriod = [crossover1Dates crossover2Dates];
closedLoopEnabledDate = []; % Removed
excludeDatesBeforeStimOn = stage2date_start:caldays(1):closedLoopEnabledDate; % Exclude dates before stim was turned on
plotDetectors = [1,2];
stage1date_start = []; % Removed
stage1date_start_relative = days(stage1date_start - enrollment);

%%
% Load redcap file
surveyData = readtable(redcapFile);
%%
% Plotting self-reported VAS / HAMD

% Remove surveys from days of study visits
allDates = surveyData.completion_pt_timestamp;
for iDate = 1:length(allDates)
    toTest = datetime(allDates(iDate),'Format','dd-MMM-uuuu');
    allDates_string(iDate) = string(toTest);
end

indicesToKeep = find(~ismember(allDates_string,[string(studyVisitDates) string(excludeDates_crossoverPeriod)]));
selectTime = allDates(indicesToKeep);
selectSurveyData = surveyData(indicesToKeep,:);

%%
% Collecting behavior about detector
detector = readtable(histogramHourly);
selectDetectorIndices = find(detector.RegionStartTime >= stage2date_start + 1);

% Clean table such that only full days are included, starting at resetTime
tempStartCounterIndices = find(hour(detector.RegionStartTime) == resetTime);
tempIndex = find(ismember(tempStartCounterIndices,selectDetectorIndices),1) ;
indexToStart = tempStartCounterIndices(tempIndex); % Day after implantation, index of detector corresponding to resetTime

tempDetector = detector(indexToStart:end,:);

clear indicesToKeep indicesToRemove
% Remove the times for calendar days that are part of the removal list --
% the method below will clean up the time before/after
allDatesToExclude = [excludeDates_crossoverPeriod,...
    studyVisitDates studyVisitDates+1 excludeDatesBeforeStimOn];
indicesToRemove = find(ismember(tempDetector.RegionStartTime,allDatesToExclude));
indicesToKeep = setdiff(1:size(tempDetector,1),indicesToRemove);
selectDetector = tempDetector(indicesToKeep,:);
clear indicesToKeep indicesToRemove

tempElapsed = datetime('now') - selectDetector.RegionStartTime;
tempDiff = diff(tempElapsed);
problemIndices = find(tempDiff ~= duration(-1,0,0));

% Find all indices where difference in time is not 1 hour (daylight savings change OR insufficient interrogation)

% For each problematic index, need to determine all indices for that 'day'
% to remove (starting at resetTime); Can't do this mathematically, b/c
% don't know how many hours are missing
allIndicesToRemove = [];

for iIndex = 1:length(problemIndices)
    currentIndex = problemIndices(iIndex);
    currentTime = selectDetector.RegionStartTime(currentIndex);
    currentHour = hour(currentTime);
    continueLoopA = 1;
    continueLoopB = 1;
    currentIndicesToRemove = [];

    % Loop through indices going backwards until you hit the most recent preceding reset time
    priorIndices = 1;
    while continueLoopA == 1
        if ~isequal(hour(selectDetector.RegionStartTime(currentIndex - priorIndices)),resetTime)
            priorIndices = priorIndices + 1;
        else
            continueLoopA = 0;
        end
    end
    % Loop through indices going forwards until you hit the next reset time
    nextIndices = 1;
    while continueLoopB == 1
        if ~isequal(hour(selectDetector.RegionStartTime(currentIndex + nextIndices)),resetTime)
            nextIndices = nextIndices + 1;
        else
            continueLoopB = 0;
        end
    end
    currentIndicesToRemove = currentIndex - priorIndices: currentIndex + nextIndices - 1;
    allIndicesToRemove = [allIndicesToRemove currentIndicesToRemove];
end
indicesToKeep = setdiff(1:size(selectDetector,1),allIndicesToRemove);

% CleanedDetector should be used for all subsequent analyses
cleanedDetector = selectDetector(indicesToKeep,:);
clear allIndicesToRemove indicesToKeep selectDetector

% Find all indices for resetTime
startCounterIndices = find(hour(cleanedDetector.RegionStartTime) == resetTime);

% Sum metrics by 24-hour period (starting at resetTime) and during waking
% hours
numWakingHours = 15; % 6am to 9pm, by including 15 bins starting at 6am
numSleepingHours = 9; % 9pm to 6am
for iIndex = 1:length(startCounterIndices)-1
    currentStartIndex = startCounterIndices(iIndex);
    allDates_startTime(iIndex) = datetime(cleanedDetector.RegionStartTime(currentStartIndex),...
        'Format','dd-MMM-uuuu');

    episodeStartsWithRX(iIndex) = sum(cleanedDetector.EpisodeStartsWithRX(currentStartIndex:startCounterIndices(iIndex+1)-1));
    episodeStartsWithRX_wakingHours(iIndex) = sum(cleanedDetector.EpisodeStartsWithRX(currentStartIndex:currentStartIndex+numWakingHours-1));
    episodeStartsWithRX_sleepingHours(iIndex) = sum(cleanedDetector.EpisodeStartsWithRX(currentStartIndex+numWakingHours:currentStartIndex+numWakingHours+numSleepingHours-1));
end

% Dates relative to enrollment
allDates_relative = days(allDates_startTime - enrollment);

%%
% MADRS timestamps
MADRS_timestamps = datetime(MADRS.Date);
for iDate = 1:length(MADRS_timestamps)
    toTest_MADRSdates = datetime(MADRS_timestamps(iDate),'Format','dd-MMM-uuuu');
    allMADRSDates_string(iDate) = string(toTest_MADRSdates);
end

% If patient has undergone crossover, exclude MADRS scores during the
% crossover
clear indicesToRemove indicesToKeep
if strcmp(patientID,'PR01')
    indicesToRemove = find(ismember(allMADRSDates_string,string(excludeDates_crossoverPeriod)));
    indicesToKeep = setdiff(1:size(MADRS_timestamps,1),indicesToRemove);
    MADRS_timestamps = MADRS_timestamps(indicesToKeep,:);
    MADRS_total = MADRS.MADRS_Total(indicesToKeep,:);

    crossoverPeriod_relative = days(excludeDates_crossoverPeriod - enrollment);
    crossoverPeriod1_relative = 698:756;
    crossoverPeriod2_relative = 1011:1022;
else
    MADRS_total = MADRS.MADRS_Total;
end

% Days elapsed relative to enrollment
MADRS_relativeTimestamps = days(MADRS_timestamps - enrollment);
selectTime_relative = days(selectTime - enrollment); % survey time
stage2date_start_relative = days(stage2date_start - enrollment);
stage2date_end_relative = days(stage2date_end - enrollment);
stage3date_start_relative = days(stage3date_start - enrollment);
ltfudate_start_relative = days(ltfudate_start - enrollment);

%%
% MADRS figure
figure
clf
% set(gcf,'Position',[262 403 1272 316])
set(gcf,'Position',[118.6000 462.6000 1272 399.2000])
scatter(MADRS_relativeTimestamps,MADRS_total,[],[0 0.4470 0.7410],'filled')
hold on
currentXLim = get(gca,'xlim');
currentYLim = get(gca,'ylim');
line(currentXLim,[6.5 6.5],'Color','black','LineWidth',1,'LineStyle',':')
line(currentXLim,[19.5 19.5],'Color','black','LineWidth',1,'LineStyle',':')
line(currentXLim,[34.5 34.5],'Color','black','LineWidth',1,'LineStyle',':')

% Indication of crossover periods
rectangle('Position',[crossoverPeriod1_relative(1) 0 crossoverPeriod1_relative(end)-crossoverPeriod1_relative(1) currentYLim(2)],...
    'FaceColor',[0.7 0.7 0.7], 'EdgeColor','none')
rectangle('Position',[crossoverPeriod2_relative(1) 0 crossoverPeriod2_relative(end)-crossoverPeriod2_relative(1) currentYLim(2)
    ],'FaceColor',[0.7 0.7 0.7], 'EdgeColor','none')

% Indication of different stage start times
line([stage1date_start_relative stage1date_start_relative],currentYLim,'Color','red','LineWidth',1,'LineStyle','--')
line([stage2date_start_relative stage2date_start_relative],currentYLim,'Color','red','LineWidth',1,'LineStyle','--')
line([stage3date_start_relative stage3date_start_relative],currentYLim,'Color','red','LineWidth',1,'LineStyle','--')
line([ltfudate_start_relative ltfudate_start_relative],currentYLim,'Color','red','LineWidth',1,'LineStyle','--')

ylabel('MADRS Score')
xlim([0 MADRS_relativeTimestamps(end)])
xlabel('Days of Trial Participation')
title([patientID ' MADRS'])

%%
% Comprehensive time range
endLim = max([selectTime_relative(end),MADRS_relativeTimestamps(end),allDates_relative(end)]);

%%
% Scatter plots: VAS-D on top; Number daily stims (24-hours) on bottom
figure
clf
set(gcf,'Position',[118.6000 332.2000 1272 529.6000])
subplot(2,1,1)
scatter(selectTime_relative,selectSurveyData.vas_depression,10,'k','filled')
hold on
currentYLim = get(gca,'ylim');
% line([stage1date_start_relative stage1date_start_relative],currentYLim,'Color','red','LineWidth',1,'LineStyle','--')
line([stage2date_start_relative stage2date_start_relative],currentYLim,'Color','red','LineWidth',1,'LineStyle','--')
line([stage3date_start_relative stage3date_start_relative],currentYLim,'Color','red','LineWidth',1,'LineStyle','--')
line([ltfudate_start_relative ltfudate_start_relative],currentYLim,'Color','red','LineWidth',1,'LineStyle','--')
xlim([0 endLim])
ylim([0 100])
ylabel('VAS Depression')

subplot(2,1,2)
scatter(allDates_relative,episodeStartsWithRX,10,'k','filled')
hold on
currentYLim = get(gca,'ylim');
% line([stage1date_start_relative stage1date_start_relative],currentYLim,'Color','red','LineWidth',1,'LineStyle','--')
line([stage2date_start_relative stage2date_start_relative],currentYLim,'Color','red','LineWidth',1,'LineStyle','--')
line([stage3date_start_relative stage3date_start_relative],currentYLim,'Color','red','LineWidth',1,'LineStyle','--')
line([ltfudate_start_relative ltfudate_start_relative],currentYLim,'Color','red','LineWidth',1,'LineStyle','--')
xlim([0 endLim])
ylabel('Stims/day')
xlabel('Days of Trial Participation')

%%
% Scatter plots: HAMD-D on top; Number daily stims (24-hours) on bottom
figure
clf
set(gcf,'Position',[118.6000 332.2000 1272 529.6000])
subplot(2,1,1)
scatter(selectTime_relative,selectSurveyData.hamd_total,10,'k','filled')
hold on
currentYLim = get(gca,'ylim');
currentXLim = get(gca,'xlim');
line([stage2date_start_relative stage2date_start_relative],currentYLim,'Color','red','LineWidth',1,'LineStyle','--')
line([stage3date_start_relative stage3date_start_relative],currentYLim,'Color','red','LineWidth',1,'LineStyle','--')
line([ltfudate_start_relative ltfudate_start_relative],currentYLim,'Color','red','LineWidth',1,'LineStyle','--')
line([0 currentXLim(2)],[7.5 7.5],'Color','black','LineWidth',1,'LineStyle',':')
line([0 currentXLim(2)],[13.5 13.5],'Color','black','LineWidth',1,'LineStyle',':')
line([0 currentXLim(2)],[18.5 18.5],'Color','black','LineWidth',1,'LineStyle',':')
line([0 currentXLim(2)],[22.5 22.5],'Color','black','LineWidth',1,'LineStyle',':')
xlim([0 endLim])
ylabel('HAMD')

subplot(2,1,2)
scatter(allDates_relative,episodeStartsWithRX,10,'k','filled')
hold on
currentYLim = get(gca,'ylim');
line([stage2date_start_relative stage2date_start_relative],currentYLim,'Color','red','LineWidth',1,'LineStyle','--')
line([stage3date_start_relative stage3date_start_relative],currentYLim,'Color','red','LineWidth',1,'LineStyle','--')
line([ltfudate_start_relative ltfudate_start_relative],currentYLim,'Color','red','LineWidth',1,'LineStyle','--')
xlim([0 endLim])
ylabel('Stims/day')
xlabel('Days of Trial Participation')

%%
% Scatter plots: VAS-D on top; Number overnight stimulations on bottom
figure
clf
set(gcf,'Position',[118.6000 332.2000 1272 529.6000])
subplot(2,1,1)
scatter(selectTime_relative,selectSurveyData.vas_depression,10,'k','filled')
hold on
currentYLim = get(gca,'ylim');
line([stage2date_start_relative stage2date_start_relative],currentYLim,'Color','red','LineWidth',1,'LineStyle','--')
line([stage3date_start_relative stage3date_start_relative],currentYLim,'Color','red','LineWidth',1,'LineStyle','--')
line([ltfudate_start_relative ltfudate_start_relative],currentYLim,'Color','red','LineWidth',1,'LineStyle','--')
xlim([0 endLim])
ylim([0 100])
ylabel('VAS Depression')

subplot(2,1,2)
scatter(allDates_relative,episodeStartsWithRX_sleepingHours,10,'k','filled')
hold on
currentYLim = get(gca,'ylim');
line([stage2date_start_relative stage2date_start_relative],currentYLim,'Color','red','LineWidth',1,'LineStyle','--')
line([stage3date_start_relative stage3date_start_relative],currentYLim,'Color','red','LineWidth',1,'LineStyle','--')
line([ltfudate_start_relative ltfudate_start_relative],currentYLim,'Color','red','LineWidth',1,'LineStyle','--')
xlim([0 endLim])
ylabel('Number overnight stimulations')
xlabel('Days of Trial Participation')

%%
% Is there a correlation between number of daytime (6am to 9pm) / overnight (9pm to 6am) / 24-hour stimulations and VAS-D score [or HAMD score]?
stimsToCalculate = episodeStartsWithRX_sleepingHours; % episodeStartsWithRX_sleepingHours / episodeStartsWithRX_wakingHours / episodeStartsWithRX
ylabelText = 'Number of overnight stimulations';
% ylabelText = 'Number of daytime stimulations';
% ylabelText = 'Number of stimulations per 24 hours';
excludeZeroOvernightStims = 1; % 0 = use all available paired days

symptomToPlot = selectSurveyData.hamd_total; % selectSurveyData.vas_depression; selectSurveyData.hamd_total
% symptomLabel = 'VAS-Depression';
symptomLabel = 'HAMD-6';

% Dates for which we have VAS-D measures (selectTime_relative paired with selectSurveyData.vas_depression)
selectTime_relative_floor = floor(selectTime_relative);
uniqueDays_relative = unique(selectTime_relative_floor);

% Dates for which we have detector information (allDates_relative paired with detector info)
allDates_relative_floor = floor(allDates_relative);

% Find shared dates, part 1
clear selectIndices paired_episodeStartsWithRX indicesToUse
selectIndices = find(ismember(allDates_relative_floor,uniqueDays_relative));

if excludeZeroOvernightStims == 1
    % Remove dates for which overnight stims == 0
    indicesToRemove = find(stimsToCalculate == 0);
    indicesToUse = setdiff(selectIndices,indicesToRemove);
else
    indicesToUse = selectIndices;

end
% Find shared dates, part 2
sharedDays_relative = allDates_relative_floor(indicesToUse);

% Number of overnight stims for shared dates
paired_episodeStartsWithRX = stimsToCalculate(indicesToUse);

% Averaged VAS-D for shared dates
clear pairedAverageSymptom
for iDay = 1:length(sharedDays_relative)
    currentDay = sharedDays_relative(iDay);
    currentIndices = find(selectTime_relative_floor == currentDay);
    pairedAverageSymptom(iDay) = mean(symptomToPlot(currentIndices));
end

[R,p] = corrcoef(paired_episodeStartsWithRX,pairedAverageSymptom,'Rows','complete');

figure
clf
scatter(pairedAverageSymptom,paired_episodeStartsWithRX,10,'k','filled')
hold on
h = lsline;
h.Color = 'r';
string = {['R = ' num2str(R(2,1)) '  p = ' num2str(p(2,1))]};
title(string)
xlabel(symptomLabel)
ylabel(ylabelText)

%%
% Is there a correlation between number of daytime (6am to 9pm) / overnight (9pm to 6am) / 24-hour stimulations and MADRS score?
excludeZeroOvernightStims = 1; % 0 = use all available paired days
numDaysStimToAverage = 7; % how many days of overnight stims should be averaged (before the date of the MADRS score)
stimsToCalculate = episodeStartsWithRX_sleepingHours; % episodeStartsWithRX_sleepingHours / episodeStartsWithRX_wakingHours / episodeStartsWithRX
ylabelText = 'Number of overnight stimulations';
% ylabelText = 'Number of daytime stimulations';
% ylabelText = 'Number of stimulations per 24 hours';

% Dates for which we have MADRS measures (MADRS_total is paired with MADRS_relativeTimestamps)
clear selectTimeMADRS_relative_floor uniqueDays_relative
selectTimeMADRS_relative_floor = floor(MADRS_relativeTimestamps);
uniqueDays_relative = unique(selectTimeMADRS_relative_floor);

% Dates for which we have detector information (allDates_relative paired with detector info)
allDates_relative_floor = floor(allDates_relative);

% Find shared dates, part 1
clear selectIndices indicesToUse paired_episodeStartsWithRX
selectIndices = find(ismember(allDates_relative_floor,uniqueDays_relative));

if excludeZeroOvernightStims == 1
    % Remove dates for which overnight stims == 0
    indicesToRemove = find(stimsToCalculate == 0);
    indicesToUse = setdiff(selectIndices,indicesToRemove);
else
    indicesToUse = selectIndices;

end
% Find shared dates, part 2
sharedDays_relative = allDates_relative_floor(indicesToUse);

% Number of overnight stims for shared dates; incorporate numDaysStimToAverage

if numDaysStimToAverage > 1
    for iData = 1:length(indicesToUse)
        clear currentIndices
        currentIndices = indicesToUse(iData) - numDaysStimToAverage + 1:indicesToUse(iData);
        paired_episodeStartsWithRX(iData) = mean(stimsToCalculate(currentIndices));
    end
else
    paired_episodeStartsWithRX = stimsToCalculate(indicesToUse);
end

% MADRS for shared dates
clear pairedAverageMADRS
for iDay = 1:length(sharedDays_relative)
    clear currentIndices
    currentDay = sharedDays_relative(iDay);
    currentIndices = find(selectTimeMADRS_relative_floor == currentDay);
    pairedAverageMADRS(iDay) = MADRS_total(currentIndices);
end

[R,p] = corrcoef(paired_episodeStartsWithRX,pairedAverageMADRS,'Rows','complete');

figure
clf
scatter(pairedAverageMADRS,paired_episodeStartsWithRX,10,'k','filled')
hold on
h = lsline;
h.Color = 'r';
string = {['R = ' num2str(R(2,1)) '  p = ' num2str(p(2,1))]};
title(string)
xlabel('MADRS')
ylabel(ylabelText)

%%
% Correlation between MADRS and VAS-D (over time)

% For each MADRS data point, average the VAS-D from the 7-days preceeding;
% start 7 days after first VAS-D measure

% Replace "selectSurveyData.vas_depression" with "selectSurveyData.hamd_total"
% to calculate and plot comparable for hamd

% selectSurveyData.vas_depression is paired with selectTime_relative_floor
startDate_relative = selectTime_relative_floor(1) + 7;

% MADRS_total is paired with MADRS_relativeTimestamps
MADRS_relativeTimestamps_floor = floor(MADRS_relativeTimestamps);

clear selectIndices
selectIndices = find(MADRS_relativeTimestamps_floor > startDate_relative);
selectMADRSscores = MADRS_total(selectIndices);
selectMADRSrelativeTimestamps = MADRS_relativeTimestamps_floor(selectIndices);

numDaysToAverage = 7;
clear averagedVASD
for iMADRS = 1:length(selectMADRSrelativeTimestamps)
    currentRelativeTime = selectMADRSrelativeTimestamps(iMADRS);

    % Find and average VAS-D scores for 7 days prior to current MADRS
    % relative time (including MADRS date)
    clear currentRangeDatesForVASD currentVASDindices
    currentRangeDatesForVASD = currentRelativeTime-numDaysToAverage-1:currentRelativeTime;
    currentVASDindices = find(ismember(selectTime_relative_floor,currentRangeDatesForVASD));

    averagedVASD(iMADRS) = nanmean(selectSurveyData.vas_depression(currentVASDindices));

end

% Correlation of MADRS vs VAS-D -- adding additional points over time
dateChunks = 7; % days

dateCounter = 1;
accumulatingDays = dateChunks;
clear allR allP plotDates
for iDate = selectMADRSrelativeTimestamps(1):dateChunks:selectMADRSrelativeTimestamps(end)
    currentDateRange = [selectMADRSrelativeTimestamps(1) selectMADRSrelativeTimestamps(1)+accumulatingDays-1];
    % Find indices corrresponding to currentDateRange
    clear currentIndices
    currentIndices = intersect(find(selectMADRSrelativeTimestamps >= currentDateRange(1)),...
        find(selectMADRSrelativeTimestamps <= currentDateRange(2)));

    [R,p] = corrcoef(averagedVASD(currentIndices),selectMADRSscores(currentIndices),'Rows','complete');

    if length(currentIndices) == 1
        allR(dateCounter) = NaN;
        allP(dateCounter) = NaN;
    elseif ~isempty(currentIndices)
        allR(dateCounter) = R(2,1);
        allP(dateCounter) = p(2,1);
    else
        allR(dateCounter) = NaN;
        allP(dateCounter) = NaN;
    end

    plotDates(dateCounter) = currentDateRange(2); % End of date range

    accumulatingDays = accumulatingDays + dateChunks;
    dateCounter = dateCounter + 1;
end

% If p-value is > 0.05, replace allR value with NaN

nonSigIndices = find(allP >= 0.05);
sigR = allR;
sigR(nonSigIndices) = NaN;

sigIndices = find(allP < 0.05);
nonSigR = allR;
nonSigR(sigIndices) = NaN;

figure
clf
set(gcf,'Position',[118.6000 597.8000 1272 264])
plot(plotDates,sigR,'k','LineWidth',1.5)
hold on
% plot(plotDates,nonSigR,'b')
ylim([0 1.1])
currentYLim = get(gca,'ylim');
ylabel('Correlation Coefficient')
xlabel('Days of Trial Participation')
xlim([0 endLim])
line([stage2date_start_relative stage2date_start_relative],currentYLim,'Color','red','LineWidth',1,'LineStyle','--')
line([stage3date_start_relative stage3date_start_relative],currentYLim,'Color','red','LineWidth',1,'LineStyle','--')
line([ltfudate_start_relative ltfudate_start_relative],currentYLim,'Color','red','LineWidth',1,'LineStyle','--')

%%
% Correlation between MADRS and VAS-D during different trial periods

% Correlation for full time period
% Both VASD and selectMADRSscores timestamps == selectMADRSrelativeTimestamps
[R,p] = corrcoef(averagedVASD,selectMADRSscores,'Rows','complete')
figure
clf
set(gcf,'Position',[488 672.2000 321 185.8000])
scatter(averagedVASD,selectMADRSscores,10,'k','filled')
hold on
h = lsline;
h.Color = 'r';
string = {['R = ' num2str(R(2,1)) '  p = ' num2str(p(2,1))]};
title([{'All data: VAS-D and MADRS correlation'};string])
xlabel('VAS-Depression')
ylabel('MADRS')
xlim([0 100])
ylim([0 45])
axis square

% Stage 2
stage2_selectIndices = intersect(find(selectMADRSrelativeTimestamps >= stage2date_start_relative),...
    find(selectMADRSrelativeTimestamps <= stage2date_end_relative));
[R,p] = corrcoef(averagedVASD(stage2_selectIndices),selectMADRSscores(stage2_selectIndices),'Rows','complete')
figure
clf
set(gcf,'Position',[488 672.2000 321 185.8000])
scatter(averagedVASD(stage2_selectIndices),selectMADRSscores(stage2_selectIndices),10,'k','filled')
hold on
h = lsline;
h.Color = 'r';
string = {['R = ' num2str(R(2,1)) '  p = ' num2str(p(2,1))]};
title([{'Stage 2: VAS-D and MADRS correlation'};string])
xlabel('VAS-Depression')
ylabel('MADRS')
xlim([0 100])
ylim([0 45])
axis square

% Stage 3
stage3_selectIndices = intersect(find(selectMADRSrelativeTimestamps >= stage3date_start_relative),...
    find(selectMADRSrelativeTimestamps < ltfudate_start_relative));
[R,p] = corrcoef(averagedVASD(stage3_selectIndices),selectMADRSscores(stage3_selectIndices),'Rows','complete')
figure
clf
set(gcf,'Position',[488 672.2000 321 185.8000])
scatter(averagedVASD(stage3_selectIndices),selectMADRSscores(stage3_selectIndices),10,'k','filled')
hold on
h = lsline;
h.Color = 'r';
string = {['R = ' num2str(R(2,1)) '  p = ' num2str(p(2,1))]};
title([{'Stage 3: VAS-D and MADRS correlation'};string])
xlabel('VAS-Depression')
ylabel('MADRS')
xlim([0 100])
ylim([0 45])
axis square

% LTFU
LTFU_selectIndices = intersect(find(selectMADRSrelativeTimestamps >= ltfudate_start_relative),...
    find(selectMADRSrelativeTimestamps <= days(datetime('now') - enrollment)));
[R,p] = corrcoef(averagedVASD(LTFU_selectIndices),selectMADRSscores(LTFU_selectIndices),'Rows','complete')
figure
clf
set(gcf,'Position',[488 672.2000 321 185.8000])
scatter(averagedVASD(LTFU_selectIndices),selectMADRSscores(LTFU_selectIndices),10,'k','filled')
hold on
h = lsline;
h.Color = 'r';
string = {['R = ' num2str(R(2,1)) '  p = ' num2str(p(2,1))]};
title([{'LTFU: VAS-D and MADRS correlation'};string])
xlabel('VAS-Depression')
ylabel('MADRS')
xlim([0 100])
ylim([0 45])
axis square

%%
% MADRS averages for different trial periods

stage2_MADRSindices = intersect(find(MADRS_timestamps >= stage2date_start),find(MADRS_timestamps <= stage2date_end));
stage3_MADRSindices = intersect(find(MADRS_timestamps >= stage3date_start),find(MADRS_timestamps < ltfudate_start));
ltfu_MADRSindices = find(MADRS_timestamps >= ltfudate_start);

if strcmp(patientID,'PR01')
    % Using 2020-09-16 as the end of the Stage 2 optimization period
    optimizationDate = datetime('2020-09-16');
    stage2_duringOptimization_MADRSindices = intersect(find(MADRS_timestamps >= stage2date_start),find(MADRS_timestamps <= optimizationDate));
    stage2_afterOptimization_MADRSindices = intersect(find(MADRS_timestamps > optimizationDate),find(MADRS_timestamps <= stage2date_end));
end

% MADRS summary statistics
meanMADRS = mean(MADRS_total(stage2_MADRSindices))
stdevMADRS = std(MADRS_total(stage2_MADRSindices))

meanMADRS = mean(MADRS_total(stage2_duringOptimization_MADRSindices))
stdevMADRS = std(MADRS_total(stage2_duringOptimization_MADRSindices))

meanMADRS = mean(MADRS_total(stage2_afterOptimization_MADRSindices))
stdevMADRS = std(MADRS_total(stage2_afterOptimization_MADRSindices))

meanMADRS = mean(MADRS_total(stage3_MADRSindices))
stdevMADRS = std(MADRS_total(stage3_MADRSindices))

meanMADRS = mean(MADRS_total(ltfu_MADRSindices))
stdevMADRS = std(MADRS_total(ltfu_MADRSindices))

%%
% HAMD averages for different trial periods

stage2_RedcapIndices = intersect(find(selectTime >= stage2date_start),find(selectTime <= stage2date_end));
stage3_RedcapIndices = intersect(find(selectTime >= stage3date_start),find(selectTime < ltfudate_start));
ltfu_RedcapIndices = find(selectTime >= ltfudate_start);

if strcmp(patientID,'PR01')
    % Using 2020-09-16 as the end of the Stage 2 optimization period
    optimizationDate = datetime('2020-09-16');
    stage2_duringOptimization_RedcapIndices = intersect(find(selectTime >= stage2date_start),find(selectTime <= optimizationDate));
    stage2_afterOptimization_RedcapIndices = intersect(find(selectTime > optimizationDate),find(selectTime <= stage2date_end));
end

% HAMD summary statistics
meanHAMD = nanmean(selectSurveyData.hamd_total(stage2_RedcapIndices))
stdevHAMD = nanstd(selectSurveyData.hamd_total(stage2_RedcapIndices))

meanHAMD = nanmean(selectSurveyData.hamd_total(stage3_RedcapIndices))
stdevHAMD = nanstd(selectSurveyData.hamd_total(stage3_RedcapIndices))

meanHAMD = nanmean(selectSurveyData.hamd_total(ltfu_RedcapIndices))
stdevHAMD = nanstd(selectSurveyData.hamd_total(ltfu_RedcapIndices))

meanHAMD = nanmean(selectSurveyData.hamd_total(stage2_duringOptimization_RedcapIndices))
stdevHAMD = nanstd(selectSurveyData.hamd_total(stage2_duringOptimization_RedcapIndices))

meanHAMD = nanmean(selectSurveyData.hamd_total(stage2_afterOptimization_RedcapIndices))
stdevHAMD = nanstd(selectSurveyData.hamd_total(stage2_afterOptimization_RedcapIndices))


% VASD summary statistics
meanVASD = nanmean(selectSurveyData.vas_depression(stage2_RedcapIndices))
stdevVASD = nanstd(selectSurveyData.vas_depression(stage2_RedcapIndices))

meanVASD = nanmean(selectSurveyData.vas_depression(stage3_RedcapIndices))
stdevVASD = nanstd(selectSurveyData.vas_depression(stage3_RedcapIndices))

meanVASD = nanmean(selectSurveyData.vas_depression(ltfu_RedcapIndices))
stdevVASD = nanstd(selectSurveyData.vas_depression(ltfu_RedcapIndices))

meanVASD = nanmean(selectSurveyData.vas_depression(stage2_duringOptimization_RedcapIndices))
stdevVASD = nanstd(selectSurveyData.vas_depression(stage2_duringOptimization_RedcapIndices))

meanVASD = nanmean(selectSurveyData.vas_depression(stage2_afterOptimization_RedcapIndices))
stdevVASD = nanstd(selectSurveyData.vas_depression(stage2_afterOptimization_RedcapIndices))

%%
% Total electrical energy delivered (TEED) calculation

% Power in 1 second = (current)^2 x pulsewidth x frequency x resistance
% Units above: Watts, Amps, Seconds, Hz, Ohms
% Ref: https://www.brainstimjrnl.com/article/S1935-861X(20)30218-7/fulltext

pulsewidth = 0.00012;  % in seconds == 120 microseconds
stimFreq = 100; % in Hz
secsPerEpisode = 6;
resistance = 632; % KS: Average of impedances from VCVS3-4 contacts (in Ohms)

% For each day in allDays_relative (days for which we have stim counts),
% create a vector with the stimulation amplitude; most of the time, stim amp = 1
stimAmp = ones(length(allDates_relative),1);

stimAmpChangeDays = [datetime('2022-11-30'),datetime('2023-01-11'),datetime('2023-01-12')];
stimAmpChanges = [1.2 1.5 1];
stimAmpChangeDays_relative = days(stimAmpChangeDays - enrollment);

for iChange = 1:length(stimAmpChangeDays_relative) - 1
    currentChangeDate_relative = stimAmpChangeDays_relative(iChange);
    clear selectIndices
    selectIndices = intersect(find(allDates_relative >= currentChangeDate_relative),find(allDates_relative < stimAmpChangeDays_relative(iChange + 1)));
    stimAmp(selectIndices) = stimAmpChanges(iChange);
end

% Plot the TEED for each day of trial participation
calculatedPower = stimAmp.^2 .* pulsewidth .* stimFreq .* resistance .* secsPerEpisode .* episodeStartsWithRX';

figure
clf
set(gcf,'Position',[118.6000 332.2000 1272 292.8000])
scatter(allDates_relative,calculatedPower,10,'k','filled')
hold on
currentYLim = get(gca,'ylim');
currentXLim = get(gca,'xlim');
line([stage2date_start_relative stage2date_start_relative],currentYLim,'Color','red','LineWidth',1,'LineStyle','--')
line([stage3date_start_relative stage3date_start_relative],currentYLim,'Color','red','LineWidth',1,'LineStyle','--')
line([ltfudate_start_relative ltfudate_start_relative],currentYLim,'Color','red','LineWidth',1,'LineStyle','--')
xlim([0 endLim])
ylabel('TEED')
xlabel('Days of Trial Participation')

%%
% Is there a correlation between TEED and VAS-D score [or HAMD score]?

symptomToPlot = selectSurveyData.hamd_total; % selectSurveyData.vas_depression; selectSurveyData.hamd_total
% symptomLabel = 'VAS-Depression';
symptomLabel = 'HAMD-6';

% Dates for which we have survey measures (selectTime_relative paired with selectSurveyData.vas_depression)
selectTime_relative_floor = floor(selectTime_relative);
uniqueDays_relative = unique(selectTime_relative_floor);

% Dates for which we have stim [TEED] information (allDates_relative paired with detector info)
allDates_relative_floor = floor(allDates_relative);

% Find shared dates, part 1
clear selectIndices paired_episodeStartsWithRX
selectIndices = find(ismember(allDates_relative_floor,uniqueDays_relative));

% Find shared dates, part 2
sharedDays_relative = allDates_relative_floor(selectIndices);

% TEED for shared dates
paired_TEED = calculatedPower(selectIndices);

% Averaged VAS-D/HAMD for shared dates
clear pairedAverageSymptom
for iDay = 1:length(sharedDays_relative)
    currentDay = sharedDays_relative(iDay);
    currentIndices = find(selectTime_relative_floor == currentDay);
    pairedAverageSymptom(iDay) = mean(symptomToPlot(currentIndices));
end

[R,p] = corrcoef(paired_TEED,pairedAverageSymptom,'Rows','complete');

figure
clf
scatter(pairedAverageSymptom,paired_TEED,10,'k','filled')
hold on
h = lsline;
h.Color = 'r';
string = {['R = ' num2str(R(2,1)) '  p = ' num2str(p(2,1))]};
title(string)
xlabel(symptomLabel)
ylabel('TEED')

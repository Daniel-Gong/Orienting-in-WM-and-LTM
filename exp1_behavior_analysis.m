%% Analysis Experiment 1: put together all the analyses required by Experiment 1

% condition matrix

% 1 cueCondition = 1:3 % WM; LTM; baseline
% 2 cueLoc = 1:2
% 3 delayTime = 1:2 % both 1,2 are 800 ms
% 4 taskType = 1:2 %1;perception; 2: memory
% 5-6 wmColorsArray
% 7-10 perceptionCharArray
% 11 perceptionProbeLocations
% 12 testPatternParts
%% prepare
clc;
clear;
close all;
currentwd = pwd();
identifier = '2023Exp1';
dataDir = 'exp1_participant_data/';
figDir = 'figures/';
matDir = 'matfiles/';
% Use filesep instead of hard-coding '/'
dataFiles = dir(fullfile(dataDir, '*mat'));
% Get the names of the files and convert to string
dataFiles = string({dataFiles.name}');
% Trim the strings
dataFiles = strtrim(dataFiles);
nSubjects = length(dataFiles);
rng('Shuffle');
addpath('scripts');
cd(dataDir);
subCountWM = zeros(1,4);subCountLTM = zeros(1,4);
testNumbers = 1:480;
colors = [65, 143, 111;138, 118, 190]/255; 
extra_colors = cbrewer('qual', 'Set2', 10);
%% load the data
for nSub = 1:nSubjects
    load(dataFiles{nSub});
    testRT = testRTs; % * 1000;
    learnRT = learnRTs; % * 1000;
    % sub info
    subNumbers(nSub) = subNumber;
    subGender(nSub) = str2double(subInfo{3});
    subAge(nSub) = str2double(subInfo{4});
    subHandedness(nSub) = str2double(subInfo{5});
    neutralCueTrials(:,nSub) = testNumbers(conditionMatrix(1,:)==3);
    neutralCueLocs(:,nSub) = testTruth(conditionMatrix(1,:)==3);
    
%% cleaning data
    [idx,~,cleanedMemoryRT] = exclude_n_std(testRT(conditionMatrix(4,:)==2),3,1);
    memoryACC = testAccuracy(conditionMatrix(4,:)==2);
    cleanedMemoryACC = memoryACC(idx);
    memoryMatrix = conditionMatrix(:,conditionMatrix(4,:)==2);
    cleanedMemoryMatrix = memoryMatrix(:,idx);
    memoryTestTruth = testTruth(conditionMatrix(4,:)==2);
    cleanedMemoryTestTruth = memoryTestTruth(idx);
    excludedMemoryTrials(nSub) = 240-sum(idx);
    
    [idx,~,cleanedPerceptRT] = exclude_n_std(testRT(conditionMatrix(4,:)==1),3,1);
    perceptACC = testAccuracy(conditionMatrix(4,:)==1);
    cleanedPerceptACC = perceptACC(idx);
    perceptMatrix = conditionMatrix(:,conditionMatrix(4,:)==1);
    cleanedPerceptMatrix = perceptMatrix(:,idx);
    excludedPerceptTrials(nSub) = 240-sum(idx);  
    
    allExcludedTrials(nSub) = excludedMemoryTrials(nSub)+excludedPerceptTrials(nSub);
    allExcludedPercentages(nSub) = allExcludedTrials(nSub)/480;
%% memory data
    % ***** acc at WM & LTM locations - collapsed
    accWM_memo_notCued = cleanedMemoryACC(cleanedMemoryMatrix(4,:)==2 & cleanedMemoryMatrix(1,:)==3 & (cleanedMemoryTestTruth==wmLocs(1) | cleanedMemoryTestTruth==wmLocs(2)));
    accWM_memo_cued = cleanedMemoryACC(cleanedMemoryMatrix(4,:)==2 & cleanedMemoryMatrix(1,:)==1);
    accLTM_memo_notCued = cleanedMemoryACC(cleanedMemoryMatrix(4,:)==2 & cleanedMemoryMatrix(1,:)==3 & (cleanedMemoryTestTruth==ltmLocs(1)  | cleanedMemoryTestTruth==ltmLocs(2)));
    accLTM_memo_cued = cleanedMemoryACC(cleanedMemoryMatrix(4,:)==2 & cleanedMemoryMatrix(1,:)==2);
    
    % AVG acc collapsed
    mean_accWM_memo_notCued(nSub) = mean(accWM_memo_notCued);
    mean_accWM_memo_cued(nSub) = mean(accWM_memo_cued);
    mean_accLTM_memo_notCued(nSub) = mean(accLTM_memo_notCued);
    mean_accLTM_memo_cued(nSub) = mean(accLTM_memo_cued);
    
    % ***** rt at WM & LTM locations - collapsed
    rtWM_memo_notCued = cleanedMemoryRT(cleanedMemoryMatrix(4,:)==2 & cleanedMemoryMatrix(1,:)==3 & (cleanedMemoryTestTruth==wmLocs(1) | cleanedMemoryTestTruth==wmLocs(2)));
    rtWM_memo_cued = cleanedMemoryRT(cleanedMemoryMatrix(4,:)==2 & cleanedMemoryMatrix(1,:)==1);
    rtLTM_memo_notCued = cleanedMemoryRT(cleanedMemoryMatrix(4,:)==2 & cleanedMemoryMatrix(1,:)==3 & (cleanedMemoryTestTruth==ltmLocs(1)  | cleanedMemoryTestTruth==ltmLocs(2)));
    rtLTM_memo_cued = cleanedMemoryRT(cleanedMemoryMatrix(4,:)==2 & cleanedMemoryMatrix(1,:)==2);
    
    % AVG rt collapsed
    mean_rtWM_memo_notCued(nSub) = mean(rtWM_memo_notCued);
    mean_rtWM_memo_cued(nSub) = mean(rtWM_memo_cued);
    mean_rtLTM_memo_notCued(nSub) = mean(rtLTM_memo_notCued);
    mean_rtLTM_memo_cued(nSub) = mean(rtLTM_memo_cued);    
    
%% perception data
    % ***** acc at WM & LTM locations - collapsed
    accWM_percept_notCued = cleanedPerceptACC(cleanedPerceptMatrix(4,:)==1 & cleanedPerceptMatrix(1,:)==3 & (cleanedPerceptMatrix(11,:)==wmLocs(1) | cleanedPerceptMatrix(11,:)==wmLocs(2)));
    accWM_percept_cued = cleanedPerceptACC(cleanedPerceptMatrix(4,:)==1 & cleanedPerceptMatrix(1,:)==1 & ((cleanedPerceptMatrix(2,:)==1 & cleanedPerceptMatrix(11,:)==wmLocs(1)) | (cleanedPerceptMatrix(2,:)==2 & cleanedPerceptMatrix(11,:)==wmLocs(2))));
    accLTM_percept_notCued = cleanedPerceptACC(cleanedPerceptMatrix(4,:)==1 & cleanedPerceptMatrix(1,:)==3 & (cleanedPerceptMatrix(11,:)==ltmLocs(1)  | cleanedPerceptMatrix(11,:)==ltmLocs(2)));
    accLTM_percept_cued = cleanedPerceptACC(cleanedPerceptMatrix(4,:)==1 & cleanedPerceptMatrix(1,:)==2 & ((cleanedPerceptMatrix(2,:)==1 & cleanedPerceptMatrix(11,:)==ltmLocs(1)) | (cleanedPerceptMatrix(2,:)==2 & cleanedPerceptMatrix(11,:)==ltmLocs(2))));

    % AVG acc collapsed
    mean_accWM_percept_notCued(nSub) = mean(accWM_percept_notCued);
    mean_accWM_percept_cued(nSub) = mean(accWM_percept_cued);
    mean_accLTM_percept_notCued(nSub) = mean(accLTM_percept_notCued);
    mean_accLTM_percept_cued(nSub) = mean(accLTM_percept_cued);
    
    % AVG acc benefit
    mean_accBenefit_WM(nSub) = mean_accWM_percept_cued(nSub) - mean_accWM_percept_notCued(nSub);
    mean_accBenefit_LTM(nSub) = mean_accLTM_percept_cued(nSub) - mean_accLTM_percept_notCued(nSub);
 
    % ***** rt at WM & LTM locations - collapsed
    rtWM_percept_notCued = cleanedPerceptRT(cleanedPerceptMatrix(4,:)==1 & cleanedPerceptMatrix(1,:)==3 & (cleanedPerceptMatrix(11,:)==wmLocs(1) | cleanedPerceptMatrix(11,:)==wmLocs(2)));
    rtWM_percept_cued = cleanedPerceptRT(cleanedPerceptMatrix(4,:)==1 & cleanedPerceptMatrix(1,:)==1 & ((cleanedPerceptMatrix(2,:)==1 & cleanedPerceptMatrix(11,:)==wmLocs(1)) | (cleanedPerceptMatrix(2,:)==2 & cleanedPerceptMatrix(11,:)==wmLocs(2))));
    rtLTM_percept_notCued = cleanedPerceptRT(cleanedPerceptMatrix(4,:)==1 & cleanedPerceptMatrix(1,:)==3 & (cleanedPerceptMatrix(11,:)==ltmLocs(1)  | cleanedPerceptMatrix(11,:)==ltmLocs(2)));
    rtLTM_percept_cued = cleanedPerceptRT(cleanedPerceptMatrix(4,:)==1 & cleanedPerceptMatrix(1,:)==2 & ((cleanedPerceptMatrix(2,:)==1 & cleanedPerceptMatrix(11,:)==ltmLocs(1)) | (cleanedPerceptMatrix(2,:)==2 & cleanedPerceptMatrix(11,:)==ltmLocs(2))));
    
    % AVG rt collapsed
    mean_rtWM_percept_notCued(nSub) = mean(rtWM_percept_notCued);
    mean_rtWM_percept_cued(nSub) = mean(rtWM_percept_cued);
    mean_rtLTM_percept_notCued(nSub) = mean(rtLTM_percept_notCued);
    mean_rtLTM_percept_cued(nSub) = mean(rtLTM_percept_cued);
    
    % newly added
    mean_accWM_percept(nSub) = mean([accWM_percept_notCued accWM_percept_cued]);
    mean_accLTM_percept(nSub) = mean([accLTM_percept_notCued accLTM_percept_cued]);
    
%% congruency effects
    temp_rt_1 = cleanedPerceptRT(cleanedPerceptMatrix(4,:)==1 & cleanedPerceptMatrix(1,:)==2 & ((cleanedPerceptMatrix(2,:)==2 & cleanedPerceptMatrix(11,:)==ltmLocs(1)) | (cleanedPerceptMatrix(2,:)==1 & cleanedPerceptMatrix(11,:)==ltmLocs(2))));
    perceptProbeLtmCueTheOtherLtmRT(nSub) = mean(temp_rt_1); 
    temp_acc_1 = cleanedPerceptACC(cleanedPerceptMatrix(4,:)==1 & cleanedPerceptMatrix(1,:)==2 & ((cleanedPerceptMatrix(2,:)==2 & cleanedPerceptMatrix(11,:)==ltmLocs(1)) | (cleanedPerceptMatrix(2,:)==1 & cleanedPerceptMatrix(11,:)==ltmLocs(2))));
    perceptProbeLtmCueTheOtherLtmACC(nSub) = mean(temp_acc_1);
    
    temp_rt_2 = cleanedPerceptRT(cleanedPerceptMatrix(4,:)==1 & cleanedPerceptMatrix(1,:)==1 & (cleanedPerceptMatrix(11,:)==ltmLocs(1) |  cleanedPerceptMatrix(11,:)==ltmLocs(2)));
    perceptProbeLtmCueWmRT(nSub) = mean(temp_rt_2); 
    temp_acc_2 = cleanedPerceptACC(cleanedPerceptMatrix(4,:)==1 & cleanedPerceptMatrix(1,:)==1 & (cleanedPerceptMatrix(11,:)==ltmLocs(1) |  cleanedPerceptMatrix(11,:)==ltmLocs(2)));
    perceptProbeLtmCueWmACC(nSub) = mean(temp_acc_2); 
    
    temp_rt_3 = cleanedPerceptRT(cleanedPerceptMatrix(4,:)==1 & cleanedPerceptMatrix(1,:)==1 & ((cleanedPerceptMatrix(2,:)==2 & cleanedPerceptMatrix(11,:)==wmLocs(1)) | (cleanedPerceptMatrix(2,:)==1 & cleanedPerceptMatrix(11,:)==wmLocs(2))));
    perceptProbeWmCueTheOtherWmRT(nSub) = mean(temp_rt_3); 
    temp_acc_3 = cleanedPerceptACC(cleanedPerceptMatrix(4,:)==1 & cleanedPerceptMatrix(1,:)==1 & ((cleanedPerceptMatrix(2,:)==2 & cleanedPerceptMatrix(11,:)==wmLocs(1)) | (cleanedPerceptMatrix(2,:)==1 & cleanedPerceptMatrix(11,:)==wmLocs(2))));
    perceptProbeWmCueTheOtherWmACC(nSub) = mean(temp_acc_3);
    
    temp_rt_4 = cleanedPerceptRT(cleanedPerceptMatrix(4,:)==1 & cleanedPerceptMatrix(1,:)==2 & (cleanedPerceptMatrix(11,:)==wmLocs(1) |  cleanedPerceptMatrix(11,:)==wmLocs(2)));
    perceptProbeWmCueLtmRT(nSub) = mean(temp_rt_4); 
    temp_acc_4 = cleanedPerceptACC(cleanedPerceptMatrix(4,:)==1 & cleanedPerceptMatrix(1,:)==2 & (cleanedPerceptMatrix(11,:)==wmLocs(1) |  cleanedPerceptMatrix(11,:)==wmLocs(2)));
    perceptProbeWmCueLtmACC(nSub) = mean(temp_acc_4);    
    
    % 29 march added
    perceptProbeWmCueInvalidRT(nSub) = mean([temp_rt_3 temp_rt_4]); 
    perceptProbeWmCueInvalidACC(nSub) = mean([temp_acc_3 temp_acc_4]);
    perceptProbeLtmCueInvalidRT(nSub) = mean([temp_rt_1 temp_rt_2]); 
    perceptProbeLtmCueInvalidACC(nSub) = mean([temp_acc_1 temp_acc_2]);
end

femaleIndex = find(subGender==1);
numFemale = sum(subGender(femaleIndex));
meanAge = mean(subAge);
stdAge = std(subAge);
disp(['numFemale:' num2str(numFemale) ' meanAge:' num2str(meanAge) ' stdAge:' num2str(stdAge)]);
disp(['meanExcludedPercentage:' num2str(mean(allExcludedPercentages)) ' stdExcludedPercentage:' num2str(std(allExcludedPercentages))]);
cd(currentwd);

if nSubjects ~= 1

    % memory recall task
    grandAvg_accWM_memo_notCued = mean(mean_accWM_memo_notCued);
    grandAvg_accWM_memo_cued = mean(mean_accWM_memo_cued);
    grandAvg_accLTM_memo_notCued = mean(mean_accLTM_memo_notCued);
    grandAvg_accLTM_memo_cued = mean(mean_accLTM_memo_cued);
    
    acc_memo_Normalized = Normalization([mean_accWM_memo_notCued' mean_accWM_memo_cued' ...
        mean_accLTM_memo_notCued' mean_accLTM_memo_cued']);
    
    se_accWM_memo_notCued = acc_memo_Normalized(1);
    se_accWM_memo_cued = acc_memo_Normalized(2);
    se_accLTM_memo_notCued = acc_memo_Normalized(3);
    se_accLTM_memo_cued = acc_memo_Normalized(4);
    
    % testing data - rt collapsed
    grandAvg_rtWM_memo_notCued = mean(mean_rtWM_memo_notCued);
    grandAvg_rtWM_memo_cued = mean(mean_rtWM_memo_cued);
    grandAvg_rtLTM_memo_notCued = mean(mean_rtLTM_memo_notCued);
    grandAvg_rtLTM_memo_cued = mean(mean_rtLTM_memo_cued);
    
    rt_memo_Normalized = Normalization([mean_rtWM_memo_notCued' mean_rtWM_memo_cued' ...
        mean_rtLTM_memo_notCued' mean_rtLTM_memo_cued']);
    
    se_rtWM_memo_notCued = rt_memo_Normalized(1);
    se_rtWM_memo_cued = rt_memo_Normalized(2);
    se_rtLTM_memo_notCued = rt_memo_Normalized(3);
    se_rtLTM_memo_cued = rt_memo_Normalized(4);

    % perceptual task
    grandAvg_accWM_percept_notCued = mean(mean_accWM_percept_notCued);
    grandAvg_accWM_percept_cued = mean(mean_accWM_percept_cued);
    grandAvg_accLTM_percept_notCued = mean(mean_accLTM_percept_notCued);
    grandAvg_accLTM_percept_cued = mean(mean_accLTM_percept_cued);
    
    acc_percept_Normalized = Normalization([mean_accWM_percept_notCued' mean_accWM_percept_cued' ...
        mean_accLTM_percept_notCued' mean_accLTM_percept_cued']);
    
    se_accWM_percept_notCued = acc_percept_Normalized(1);
    se_accWM_percept_cued = acc_percept_Normalized(2);
    se_accLTM_percept_notCued = acc_percept_Normalized(3);
    se_accLTM_percept_cued = acc_percept_Normalized(4);
    
    % perceptual task - accBenefit
    grandAvg_accBenefit_WM = mean(mean_accBenefit_WM);
    grandAvg_accBenefit_LTM = mean(mean_accBenefit_LTM);
    se_accBenefit_WM = std(mean_accBenefit_WM)/sqrt(length(mean_accBenefit_WM));
    se_accBenefit_LTM = std(mean_accBenefit_LTM)/sqrt(length(mean_accBenefit_LTM));
    
    % perceptual task - rt collapsed
    grandAvg_rtWM_percept_notCued = mean(mean_rtWM_percept_notCued);
    grandAvg_rtWM_percept_cued = mean(mean_rtWM_percept_cued);
    grandAvg_rtLTM_percept_notCued = mean(mean_rtLTM_percept_notCued);
    grandAvg_rtLTM_percept_cued = mean(mean_rtLTM_percept_cued);
    
    rt_percept_Normalized = Normalization([mean_rtWM_percept_notCued' mean_rtWM_percept_cued' ...
        mean_rtLTM_percept_notCued' mean_rtLTM_percept_cued']);
    
    se_rtWM_percept_notCued = rt_percept_Normalized(1);
    se_rtWM_percept_cued = rt_percept_Normalized(2);
    se_rtLTM_percept_notCued = rt_percept_Normalized(3);
    se_rtLTM_percept_cued = rt_percept_Normalized(4);
    
    % newly added
    grandAvg_accWM_percept = mean(mean_accWM_percept);
    grandAvg_accLTM_percept = mean(mean_accLTM_percept);
    se_accWM_percept = std(mean_accWM_percept)/sqrt(length(mean_accWM_percept));
    se_accLTM_percept = std(mean_accLTM_percept)/sqrt(length(mean_accLTM_percept));
%% congruency effects
    congruencyProbeLtmAccNorm = Normalization([mean_accLTM_percept_notCued' mean_accLTM_percept_cued' perceptProbeLtmCueTheOtherLtmACC' perceptProbeLtmCueWmACC']);
    congruencyProbeLtmRtNorm = Normalization([mean_rtLTM_percept_notCued' mean_rtLTM_percept_cued' perceptProbeLtmCueTheOtherLtmRT' perceptProbeLtmCueWmRT']);
    congruencyProbeWmAccNorm = Normalization([mean_accWM_percept_notCued' mean_accWM_percept_cued' perceptProbeWmCueTheOtherWmACC' perceptProbeWmCueLtmACC']);
    congruencyProbeWmRtNorm = Normalization([mean_rtWM_percept_notCued' mean_rtWM_percept_cued' perceptProbeWmCueTheOtherWmRT' perceptProbeWmCueLtmRT']);

    % 29 March New
    ValidtyAccNorm = Normalization([mean_accWM_percept_notCued' perceptProbeWmCueInvalidACC' mean_accLTM_percept_notCued' perceptProbeLtmCueInvalidACC']);
    ValidtyRtNorm = Normalization([mean_rtWM_percept_notCued' perceptProbeWmCueInvalidRT' mean_rtLTM_percept_notCued' perceptProbeLtmCueInvalidRT']);
end

%% Memory ACC
f = figure;
figureStartup;
subplot(2,2,2);
f.WindowState = 'maximized';
hold on
b = bar(100*[grandAvg_accWM_memo_notCued,grandAvg_accWM_memo_cued; grandAvg_accLTM_memo_notCued, grandAvg_accLTM_memo_cued], 'EdgeColor', 'none', 'BarWidth', 0.8);%, grandAvg_accWM_memo_cued];%,'s-','Color',colors(3,:),'MarkerFaceColor',colors(3,:),'LineWidth',2,'MarkerSize',7.3);
b(2).FaceColor = 'flat';
b(2).CData = [colors(1,:);colors(2,:)];
hold on
bb1 = bar(100*[grandAvg_accWM_memo_notCued,NaN; NaN, NaN],'BarWidth', 0.8);
bb1(1).EdgeColor = colors(1,:);
bb1(1).FaceColor = [1 1 1];
bb1(1).LineWidth = 1;
hold on
bb2 = bar(100*[NaN,NaN; grandAvg_accLTM_memo_notCued,NaN],'BarWidth', 0.8);
bb2(1).EdgeColor = colors(2,:);
bb2(1).FaceColor = [1 1 1];
bb2(1).LineWidth = 1;
hold on
for i = 1:nSubjects
    % WM
    plot([b(1).XEndPoints(1) b(2).XEndPoints(1)],100*[mean_accWM_memo_notCued(i) mean_accWM_memo_cued(i)],'-','Color',extra_colors(10,:));
%     scatter([b(1).XEndPoints(1) b(2).XEndPoints(1)],[memoryNoCueProbeWmACC(i) memoryCueWmACC(i)],[],extra_colors(10,:),'filled');
    hold on
    % LTM
    plot([b(1).XEndPoints(2) b(2).XEndPoints(2)],100*[mean_accLTM_memo_notCued(i) mean_accLTM_memo_cued(i)],'-','Color',extra_colors(10,:));
%     scatter([b(1).XEndPoints(2) b(2).XEndPoints(2)],[memoryNoCueProbeLtmACC(i) memoryCueLtmACC(i)],[],extra_colors(10,:),'filled');
    hold on
end
h = ploterr([b(1).XEndPoints b(2).XEndPoints] , [b(1).YEndPoints b(2).YEndPoints], [],100*[se_accWM_memo_notCued, se_accLTM_memo_notCued,se_accWM_memo_cued,se_accLTM_memo_cued], 'k.', 'abshhxy', 0);
set(h(1), 'marker', 'none','LineWidth',2);
hold on   
set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'WM','LTM'}); 
xlim([0.5 2.5]);
ylim([70 100]);
ylabel('Accuracy (%)');

CC = [mean_accWM_memo_notCued'  mean_accWM_memo_cued' mean_accLTM_memo_notCued' mean_accLTM_memo_cued'];
Ncompare = nchoosek(1:4,2);
for n = 1:size(Ncompare,1)
    tstats{n} = mes(CC(:,Ncompare(n,2)),CC(:,Ncompare(n,1)),'hedgesg','isDep',1);
    p_temp(n) = tstats{n}.t.p;
end
groups = mat2cell(nchoosek([b(1).XEndPoints(1) b(2).XEndPoints(1) b(1).XEndPoints(2) b(2).XEndPoints(2)],2),ones(1,size(Ncompare,1)));
% H=sigstar(groups,p_temp*2,1);
list_comparisons = [1,6]; % 1,2,5,6  are the interested comparisons

for i = list_comparisons
    mysigstar(gca,groups{i},myRange(ylim)*0.05+findMinY(groups{i}),p_temp(i));
end

% legend(b,{'LTM item','WM item'},'Location','eastoutside');
% legend boxoff
% saveas(f,[figDir identifier  '-memory-acc.emf']);
%% memory ACC two way anova
%first factor: WM/LTM; second factor: not cued / cued

groups=[repmat([1,1],nSubjects,1);
    repmat([1,2],nSubjects,1);
    repmat([2,1],nSubjects,1);
    repmat([2,2],nSubjects,1)];
anova_twoway_data = [mean_accWM_memo_notCued';
    mean_accWM_memo_cued';
    mean_accLTM_memo_notCued';
    mean_accLTM_memo_cued'];
[anovastats_memo_acc,anovatable_memo_acc]=mes2way(anova_twoway_data,groups,'partialeta2','isDep',[1 1]);
ttest_wm_memo_acc = mes(mean_accWM_memo_notCued,mean_accWM_memo_cued,'hedgesg','isDep',1);
ttest_ltm_memo_acc = mes(mean_accLTM_memo_notCued,mean_accLTM_memo_cued,'hedgesg','isDep',1);
%% Memory RT
% f = figure;
% figureStartup;
subplot(2,2,1);
f.WindowState = 'maximized';
hold on
ylim([0 2]);
% axis tight;
xlim([0.5 2.5]);
b = bar([grandAvg_rtWM_memo_notCued,grandAvg_rtWM_memo_cued; grandAvg_rtLTM_memo_notCued, grandAvg_rtLTM_memo_cued], 'EdgeColor', 'none', 'BarWidth', 0.8);%0.6
b(2).FaceColor = 'flat';
b(2).CData = [colors(1,:);colors(2,:)];
hold on
bb1 = bar([grandAvg_rtWM_memo_notCued,NaN; NaN, NaN],'BarWidth', 0.8);
bb1(1).EdgeColor = colors(1,:);
bb1(1).FaceColor = [1 1 1];
bb1(1).LineWidth = 1;
hold on
bb2 = bar([NaN,NaN; grandAvg_rtLTM_memo_notCued,NaN],'BarWidth', 0.8);
bb2(1).EdgeColor = colors(2,:);
bb2(1).FaceColor = [1 1 1];
bb2(1).LineWidth = 1;
hold on
for i = 1:nSubjects
    % WM
    plot([b(1).XEndPoints(1) b(2).XEndPoints(1)],[mean_rtWM_memo_notCued(i) mean_rtWM_memo_cued(i)],'-','Color',extra_colors(10,:));
%     scatter([b(1).XEndPoints(1) b(2).XEndPoints(1)],[memoryNoCueProbeWmACC(i) memoryCueWmACC(i)],[],extra_colors(10,:),'filled');
    hold on
    % LTM
    plot([b(1).XEndPoints(2) b(2).XEndPoints(2)],[mean_rtLTM_memo_notCued(i) mean_rtLTM_memo_cued(i)],'-','Color',extra_colors(10,:));
%     scatter([b(1).XEndPoints(2) b(2).XEndPoints(2)],[memoryNoCueProbeLtmACC(i) memoryCueLtmACC(i)],[],extra_colors(10,:),'filled');
    hold on
end

h = ploterr([b(1).XEndPoints b(2).XEndPoints] , [b(1).YEndPoints b(2).YEndPoints], [],[se_rtWM_memo_notCued, se_rtLTM_memo_notCued, se_rtWM_memo_cued, se_rtLTM_memo_cued], 'k.', 'abshhxy', 0);
set(h(1), 'marker', 'none','LineWidth',2);
hold on   
set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'WM','LTM'}); 
ylabel('RT (s)');

CC = [mean_rtWM_memo_notCued' mean_rtWM_memo_cued' mean_rtLTM_memo_notCued' mean_rtLTM_memo_cued'];
Ncompare = nchoosek(1:4,2);
for n = 1:size(Ncompare,1)
    tstats{n} = mes(CC(:,Ncompare(n,2)),CC(:,Ncompare(n,1)),'hedgesg','isDep',1);
    p_temp(n) = tstats{n}.t.p;
end
groups = mat2cell(nchoosek([b(1).XEndPoints(1) b(2).XEndPoints(1) b(1).XEndPoints(2) b(2).XEndPoints(2)],2),ones(1,size(Ncompare,1)));
% H=sigstar(groups,p_temp*2,1);
list_comparisons = [1,6];

for i = list_comparisons
    mysigstar(gca,groups{i},myRange(ylim)*0.05+findMinY(groups{i}),p_temp(i));
end

% legend(b,{'LTM item','WM item'},'Location','eastoutside');
% legend boxoff
% saveas(f,[figDir identifier  '-memory-rt.emf']);
%% memory RT two way anova
%first factor: WM/LTM; second factor: not cued / cued

groups=[repmat([1,1],nSubjects,1);
    repmat([1,2],nSubjects,1);
    repmat([2,1],nSubjects,1);
    repmat([2,2],nSubjects,1)];
anova_twoway_data = [mean_rtWM_memo_notCued';
    mean_rtWM_memo_cued';
    mean_rtLTM_memo_notCued';
    mean_rtLTM_memo_cued'];
[anovastats_memo_rt,anovatable_memo_rt]=mes2way(anova_twoway_data,groups,'partialeta2','isDep',[1 1]);

ttest_wm_memo_rt = mes(mean_rtWM_memo_notCued,mean_rtWM_memo_cued,'hedgesg','isDep',1);
ttest_ltm_memo_rt = mes(mean_rtLTM_memo_notCued,mean_rtLTM_memo_cued,'hedgesg','isDep',1);
ttest_memo_rt_benefit = mes(mean_rtLTM_memo_notCued-mean_rtLTM_memo_cued,mean_rtWM_memo_notCued-mean_rtWM_memo_cued,'hedgesg','isDep',1);
%% Perception ACC
% f = figure;
% figureStartup;
subplot(2,2,4);
f.WindowState = 'maximized';
hold on
b1 = bar(100*[NaN,mean(grandAvg_accWM_percept_cued); NaN, NaN], 'BarWidth', 0.8);
hatchfill2(b1(2),'single','HatchAngle',45,'hatchcolor',colors(1,:),'HatchLineWidth',1);
b1(2).FaceColor = 'none';
b1(2).EdgeColor = colors(1,:);
b1(2).LineWidth = 1;
hold on
b2 = bar(100*[NaN,NaN; NaN, mean(grandAvg_accLTM_percept_cued)], 'BarWidth', 0.8);
hatchfill2(b2(2),'single','HatchAngle',45,'hatchcolor',colors(2,:),'HatchLineWidth',1,'HatchOffset',0.5);
b2(2).FaceColor = 'none';
b2(2).EdgeColor = colors(2,:);
b2(2).LineWidth = 1;
hold on
bb1 = bar(100*[grandAvg_accWM_percept_notCued,NaN; NaN, NaN],'BarWidth', 0.8);
bb1(1).EdgeColor = colors(1,:);
bb1(1).FaceColor = [1 1 1];
bb1(1).LineWidth = 1;
hold on
bb2 = bar(100*[NaN,NaN; grandAvg_accLTM_percept_notCued,NaN],'BarWidth', 0.8);
bb2(1).EdgeColor = colors(2,:);
bb2(1).FaceColor = [1 1 1];
bb2(1).LineWidth = 1;
hold on
for i = 1:nSubjects
    % WM
    plot([b(1).XEndPoints(1) b(2).XEndPoints(1)],100*[mean_accWM_percept_notCued(i) mean_accWM_percept_cued(i)],'-','Color',extra_colors(10,:));
%     scatter([b(1).XEndPoints(1) b(2).XEndPoints(1)],[memoryNoCueProbeWmACC(i) memoryCueWmACC(i)],[],extra_colors(10,:),'filled');
    hold on
    % LTM
    plot([b(1).XEndPoints(2) b(2).XEndPoints(2)],100*[mean_accLTM_percept_notCued(i) mean_accLTM_percept_cued(i)],'-','Color',extra_colors(10,:));
%     scatter([b(1).XEndPoints(2) b(2).XEndPoints(2)],[memoryNoCueProbeLtmACC(i) memoryCueLtmACC(i)],[],extra_colors(10,:),'filled');
    hold on
end
h = ploterr([bb1(1).XEndPoints(1) b1(2).XEndPoints(1) bb2(1).XEndPoints(2) b2(2).XEndPoints(2)] , [bb1(1).YEndPoints(1) b1(2).YEndPoints(1) bb2(1).YEndPoints(2) b2(2).YEndPoints(2)], [],100*acc_percept_Normalized([1 2 3 4]), 'k.', 'abshhxy', 0);
set(h(1), 'marker', 'none','LineWidth',2);
hold on   
set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'WM','LTM'}); 
xlim([0.5 2.5]);
% ylim([50 80]);
ylabel('Accuracy (%)');

CC = [mean_accWM_percept_notCued' mean_accWM_percept_cued' mean_accLTM_percept_notCued'  mean_accLTM_percept_cued'];
Ncompare = nchoosek(1:4,2);
for n = 1:size(Ncompare,1)
    tstats{n} = mes(CC(:,Ncompare(n,2)),CC(:,Ncompare(n,1)),'hedgesg','isDep',1);
    p_temp(n) = tstats{n}.t.p;
end
groups = mat2cell(nchoosek([b(1).XEndPoints(1) b(2).XEndPoints(1) b(1).XEndPoints(2) b(2).XEndPoints(2)],2),ones(1,size(Ncompare,1)));
% H=sigstar(groups,p_temp*2,1);

list_comparisons = [1,6]; % 1,2,5,6  are the interested comparisons
for i = list_comparisons
    mysigstar(gca,groups{i},myRange(ylim)*0.05+findMinY(groups{i}),p_temp(i));
end

% legend(b,{'LTM location','WM location'},'Location','eastoutside');
% legend boxoff
% saveas(f,[figDir identifier  '-perception-acc.emf']);
%% perception ACC two way anova
%first factor: WM/LTM; second factor: not cued / cued

ttest_baseline = mes(mean_accWM_percept_notCued',mean_accLTM_percept_notCued','hedgesg','isDep',1);
ttest_cued = mes(mean_accWM_percept_cued',mean_accLTM_percept_cued','hedgesg','isDep',1);
ttest_uplift = mes(mean_accWM_percept_cued',mean_accLTM_percept_notCued','hedgesg','isDep',1);

groups=[repmat([1,1],nSubjects,1);
    repmat([1,2],nSubjects,1);
    repmat([2,1],nSubjects,1);
    repmat([2,2],nSubjects,1)];
anova_twoway_data = [mean_accWM_percept_notCued';
    mean_accWM_percept_cued';
    mean_accLTM_percept_notCued';
    mean_accLTM_percept_cued'];
[anovastats_percept_acc,anovatable_percept_acc]=mes2way(anova_twoway_data,groups,'partialeta2','isDep',[1 1]);

ttest_wm_percept_acc = mes(mean_accWM_percept_notCued',mean_accWM_percept_cued','hedgesg','isDep',1);
ttest_ltm_percept_acc = mes(mean_accLTM_percept_notCued',mean_accLTM_percept_cued','hedgesg','isDep',1);
ttest_nocue_wm_ltm_acc = mes(mean_accWM_percept_notCued',mean_accLTM_percept_notCued','hedgesg','isDep',1);
ttest_cued_wm_ltm_acc = mes(mean_accWM_percept_cued',mean_accLTM_percept_cued','hedgesg','isDep',1);


%%  Perception RT
% f = figure;
% figureStartup;
subplot(2,2,3);
f.WindowState = 'maximized';
hold on
% ylim([700 1200]);
% axis tight;
xlim([0.5 2.5]);
b1 = bar([NaN,mean(grandAvg_rtWM_percept_cued); NaN, NaN], 'BarWidth', 0.8);
hatchfill2(b1(2),'single','HatchAngle',45,'hatchcolor',colors(1,:),'HatchLineWidth',1);
b1(2).FaceColor = 'none';
b1(2).EdgeColor = colors(1,:);
b1(2).LineWidth = 1;
hold on
b2 = bar([NaN,NaN; NaN, mean(grandAvg_rtLTM_percept_cued)], 'BarWidth', 0.8);
hatchfill2(b2(2),'single','HatchAngle',45,'hatchcolor',colors(2,:),'HatchLineWidth',1,'HatchOffset',0.5);
b2(2).FaceColor = 'none';
b2(2).EdgeColor = colors(2,:);
b2(2).LineWidth = 1;
hold on
bb1 = bar([grandAvg_rtWM_percept_notCued,NaN; NaN, NaN],'BarWidth', 0.8);
bb1(1).EdgeColor = colors(1,:);
bb1(1).FaceColor = [1 1 1];
bb1(1).LineWidth = 1;
hold on
bb2 = bar([NaN,NaN; grandAvg_rtLTM_percept_notCued,NaN],'BarWidth', 0.8);
bb2(1).EdgeColor = colors(2,:);
bb2(1).FaceColor = [1 1 1];
bb2(1).LineWidth = 1;
hold on
  
for i = 1:nSubjects
    % WM
    plot([bb1(1).XEndPoints(1) b1(2).XEndPoints(1)],[mean_rtWM_percept_notCued(i) mean_rtWM_percept_cued(i)],'-','Color',extra_colors(10,:));
%     scatter([bb1(1).XEndPoints(1) b1(2).XEndPoints(1)],100*[perceptNoCueProbeWmACC(i) perceptCueWmProbeWmACC(i)],[],extra_colors(10,:),'filled');
    hold on
    % LTM
    plot([bb2(1).XEndPoints(2) b2(2).XEndPoints(2)],[mean_rtLTM_percept_notCued(i) mean_rtLTM_percept_cued(i)],'-','Color',extra_colors(10,:));
%     scatter([bb2(1).XEndPoints(2) b2(2).XEndPoints(2)],100*[perceptNoCueProbeLtmACC(i) perceptCueLtmProbeLtmACC(i)],[],extra_colors(10,:),'filled');
    hold on
end

h = ploterr([bb1(1).XEndPoints(1) b1(2).XEndPoints(1) bb2(1).XEndPoints(2) b2(2).XEndPoints(2)] , [bb1(1).YEndPoints(1) b1(2).YEndPoints(1) bb2(1).YEndPoints(2) b2(2).YEndPoints(2)], [],rt_percept_Normalized([1 2 3 4]), 'k.', 'abshhxy', 0);
set(h(1), 'marker', 'none','LineWidth',2);
hold on   
set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'WM','LTM'}); 
ylabel('RT (ms)');

CC = [mean_rtWM_percept_notCued' mean_rtWM_percept_cued' mean_rtLTM_percept_notCued'  mean_rtLTM_percept_cued'];
Ncompare = nchoosek(1:4,2);
for n = 1:size(Ncompare,1)
    tstats{n} = mes(CC(:,Ncompare(n,2)),CC(:,Ncompare(n,1)),'hedgesg','isDep',1);
    p_temp(n) = tstats{n}.t.p;
end
groups = mat2cell(nchoosek([b(1).XEndPoints(1) b(2).XEndPoints(1) b(1).XEndPoints(2) b(2).XEndPoints(2)],2),ones(1,size(Ncompare,1)));
% H=sigstar(groups,p_temp*2,1);
list_comparisons = [1,6];

for i = list_comparisons
    mysigstar(gca,groups{i},myRange(ylim)*0.05+findMinY(groups{i}),p_temp(i));
end

% legend(b,{'LTM location','WM location'},'Location','eastoutside');
% legend boxoff
% saveas(gca,[figDir identifier  '-perception-rt.emf']);
saveas(gca,[figDir identifier  '-2-2.svg']);
%% perception RT two way anova
%first factor: WM/LTM; second factor: not cued / cued

groups=[repmat([1,1],nSubjects,1);
    repmat([1,2],nSubjects,1);
    repmat([2,1],nSubjects,1);
    repmat([2,2],nSubjects,1)];
anova_twoway_data = [mean_rtWM_percept_notCued';
    mean_rtWM_percept_cued';
    mean_rtLTM_percept_notCued';
    mean_rtLTM_percept_cued'];
[anovastats_percept_rt,anovatable_percept_rt]=mes2way(anova_twoway_data,groups,'partialeta2','isDep',[1 1]);

ttest_wm_percept_rt = mes(mean_rtWM_percept_notCued',mean_rtWM_percept_cued','hedgesg','isDep',1);
ttest_ltm_percept_rt = mes(mean_rtLTM_percept_notCued',mean_rtLTM_percept_cued','hedgesg','isDep',1);
%% ACC benefit
f = figure;
figureStartup;
subplot(1,2,1); hold on;
f.WindowState = 'maximized';
% rather than a square plot, make it thinner
handles = violinPlot(100*mean_accBenefit_LTM', 'histOri', 'right', 'widthDiv', [2 2], 'showMM', 4, ...
    'color',  mat2cell(colors(2, : ), 1));
hatchfill2(handles{1},'single','HatchAngle',45,'HatchColor',colors(2, : ),'HatchLineWidth',1);
hold on
handles = violinPlot(100*mean_accBenefit_WM', 'histOri', 'left', 'widthDiv', [2 1], 'showMM', 4, ...
    'color',  mat2cell(colors(1, : ), 1));
hatchfill2(handles{1},'single','HatchAngle',45,'HatchColor',colors(1, : ),'HatchLineWidth',1,'HatchOffset',0.5);
set(gca, 'XTick', [0.6 1.4], 'XTickLabel', {'WM cue', 'LTM cue'});
xlim([0.2 1.8]);
ylim([-40,70]);
ylabel('Accuracy benefit (%)');
hold on
yline(0,'--k','LineWidth',1);

if nSubjects ~= 1
    ttest_accBenefit = mes(mean_accBenefit_WM',mean_accBenefit_LTM','hedgesg','isDep',1);
    ttest_onesample_WM = mes(mean_accBenefit_WM',0,'g1');
    ttest_onesample_LTM = mes(mean_accBenefit_LTM',0,'g1');
end

%add significance stars for each bar
xticks = get(gca, 'xtick');
ypos = -30; % plot below
% mysigstar(gca, xticks(1), ypos, ttest_onesample_LTM.t.p);
% mysigstar(gca, xticks(2), ypos, ttest_onesample_WM.t.p);
% mysigstar(gca, xticks, 65, ttest_accBenefit.t.p);
saveas(gca,[figDir identifier  '-perception-accBenefit.svg']);
%}


function [idx,data_mean,data_whole] = exclude_n_std(data,nstd,ntail)
    std_data = std(data);
    mean_data = mean(data);
    idx = NaN(size(data));
    if ntail == 1
        for i = 1:length(data)
            item = data(i);
            if item > mean_data + nstd * std_data
                idx(i) = 0;
            else
                idx(i) = 1;
            end
        end
    elseif ntail == 2
        for i = 1:length(data)
            item = data(i);
            if item > mean_data + nstd * std_data || item < mean_data - nstd * std_data
                idx(i) = 0;
            else
                idx(i) = 1;
            end
        end
    end
    idx = logical(idx);
    data_whole = data(idx);
    data_mean = mean(data(idx));
end    
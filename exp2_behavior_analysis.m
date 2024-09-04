%% Jan 19th, by Dongyu Gong

%% prepare
clc;
clear;
close all;
rng('Shuffle');
currentwd = pwd();
identifier = '2023Exp3';
dataDir = 'exp3_participant_data/';
figDir = 'figures/';
matDir = 'matfiles/';
% Use filesep instead of hard-coding '/'
dataFiles = dir(fullfile(dataDir, '*mat'));
% Get the names of the files and convert to string
dataFiles = string({dataFiles.name}');
% Trim the strings
dataFiles = strtrim(dataFiles);
nSubjects = length(dataFiles);
filenames = cell(1,nSubjects);

addpath('scripts');
testNumbers = 1:600;
colors = [65, 143, 111;138, 118, 190]/255; 
extra_colors = cbrewer('qual', 'Set2', 10);

%%

memoryCueWmRT_pooled = [];
memoryCueLtmRT_pooled = [];
memoryNoCueProbeWmRT_pooled = [];
memoryNoCueProbeLtmRT_pooled = [];
perceptCueWmProbeWmRT_pooled = [];
perceptCueLtmProbeLtmRT_pooled = [];
perceptNoCueProbeWmRT_pooled = [];
perceptNoCueProbeLtmRT_pooled = [];

memoryCueWmACC_pooled = [];
memoryCueLtmACC_pooled = [];
memoryNoCueProbeWmACC_pooled = [];
memoryNoCueProbeLtmACC_pooled = [];

subCountWM = zeros(1,4);subCountLTM = zeros(1,4);

for nSub = 1:nSubjects
    filenames{nSub} = ['exp3_result_' subString{nSub} '.mat'];
    load(dataFiles{nSub});
    disp(dataFiles{nSub});

    % sub info
    subGender(nSub) = str2double(subInfo{3});
    subAge(nSub) = str2double(subInfo{4});
    subHandedness(nSub) = str2double(subInfo{5}); 
    subVision(nSub) = str2double(subInfo{6});
    
    % load data
    learnMatrix{nSub} = learn.Matrix;
    % row 1: thisTrialProbeType(1=color, 2=shape);
    % row 2: thisTrialProbeShape (ltmShapes, which is an index);
    learnRT{nSub} = learn.LongRTs;
    learnAcc{nSub} = learn.Accuracy;
    
    testMatrix{nSub} = conditionMatrix;
    % row 1: cueLoc (1-4: four locations; 5:no cue)
    % row 2: taskType (1=perception, 2=memory)
    % row 3: probeArrow (1-4: four locations)
    testRT{nSub} = test.LongRTs;
    testAcc{nSub} = test.Accuracy;
    allWMLocs{nSub} = wmLocs;
    allLTMLocs{nSub} = ltmLocs;
    
    test.LongRTs = test.LongRTs; % 1000 * test.LongRTs;
    
    %% organizing testing data
    [idx,~,cleanedMemoryRT] = exclude_n_std(test.LongRTs(conditionMatrix(2,:)==2),3,1);
    memoryACC = test.Accuracy(conditionMatrix(2,:)==2);
    cleanedMemoryACC = memoryACC(idx);
    memoryMatrix = conditionMatrix(:,conditionMatrix(2,:)==2);
    cleanedMemoryMatrix = memoryMatrix(:,idx);
    memoryNoCueProbeColors = test.NoCueProbeColors(conditionMatrix(2,:)==2);
    cleanedNoCueProbeColors = memoryNoCueProbeColors(idx);
    excludedMemoryTrials(nSub) = 300-sum(idx);
    
    [idx,~,cleanedPerceptRT] = exclude_n_std(test.LongRTs(conditionMatrix(2,:)==1),3,1);
    perceptACC = test.Accuracy(conditionMatrix(2,:)==1);
    cleanedPerceptACC = perceptACC(idx);
    perceptMatrix = conditionMatrix(:,conditionMatrix(2,:)==1);
    cleanedPerceptMatrix = perceptMatrix(:,idx);
    excludedPerceptTrials(nSub) = 300-sum(idx);   
    
    allExcludedTrials(nSub) = excludedMemoryTrials(nSub)+excludedPerceptTrials(nSub);
    allExcludedPercentages(nSub) = allExcludedTrials(nSub)/600;    

    % memory task
    temp_rt = cleanedMemoryRT(cleanedMemoryMatrix(2,:)==2 & ismember(cleanedMemoryMatrix(1,:),wmLocs));
    memoryCueWmRT(nSub) = mean(temp_rt); 
    temp_acc = cleanedMemoryACC(cleanedMemoryMatrix(2,:)==2 & ismember(cleanedMemoryMatrix(1,:),wmLocs));
    memoryCueWmACC(nSub) = mean(temp_acc);
    memoryCueWmACC_pooled = [memoryCueWmACC_pooled temp_acc];
    memoryCueWmRT_pooled = [memoryCueWmRT_pooled; temp_rt'];
    
    temp_rt = cleanedMemoryRT(cleanedMemoryMatrix(2,:)==2 & ismember(cleanedMemoryMatrix(1,:),ltmLocs));
    memoryCueLtmRT(nSub) = mean(temp_rt); 
    temp_acc = cleanedMemoryACC(cleanedMemoryMatrix(2,:)==2 & ismember(cleanedMemoryMatrix(1,:),ltmLocs));
    memoryCueLtmACC(nSub) = mean(temp_acc);
    memoryCueLtmACC_pooled = [memoryCueLtmACC_pooled temp_acc];
    memoryCueLtmRT_pooled = [memoryCueLtmRT_pooled; temp_rt'];
    
    temp_rt = cleanedMemoryRT(cleanedMemoryMatrix(2,:)==2 & cleanedMemoryMatrix(1,:)==5);
    memoryNoCueRT(nSub) = mean(temp_rt); 
    temp_acc = cleanedMemoryACC(cleanedMemoryMatrix(2,:)==2 & cleanedMemoryMatrix(1,:)==5);
    memoryNoCueACC(nSub) = mean(temp_acc);
    
    temp_rt = cleanedMemoryRT(cleanedMemoryMatrix(2,:)==2 & cleanedMemoryMatrix(1,:)==5 & ismember(cleanedNoCueProbeColors,wmColors));
    memoryNoCueProbeWmRT(nSub) = mean(temp_rt); 
    temp_acc = cleanedMemoryACC(cleanedMemoryMatrix(2,:)==2 & cleanedMemoryMatrix(1,:)==5 & ismember(cleanedNoCueProbeColors,wmColors));
    memoryNoCueProbeWmACC(nSub) = mean(temp_acc);
    memoryNoCueProbeWmACC_pooled = [memoryNoCueProbeWmACC_pooled temp_acc];
    memoryNoCueProbeWmRT_pooled = [memoryNoCueProbeWmRT_pooled; temp_rt'];
    
    temp_rt = cleanedMemoryRT(cleanedMemoryMatrix(2,:)==2 & cleanedMemoryMatrix(1,:)==5 & ismember(cleanedNoCueProbeColors,ltmColors));
    memoryNoCueProbeLtmRT(nSub) = mean(temp_rt); 
    temp_acc = cleanedMemoryACC(cleanedMemoryMatrix(2,:)==2 & cleanedMemoryMatrix(1,:)==5 & ismember(cleanedNoCueProbeColors,ltmColors));
    memoryNoCueProbeLtmACC(nSub) = mean(temp_acc);
    memoryNoCueProbeLtmACC_pooled = [memoryNoCueProbeLtmACC_pooled temp_acc];
    memoryNoCueProbeLtmRT_pooled = [memoryNoCueProbeLtmRT_pooled; temp_rt'];
    
    % perception task
    temp_rt = cleanedPerceptRT(cleanedPerceptMatrix(2,:)==1 & cleanedPerceptMatrix(1,:)==5 & ismember(cleanedPerceptMatrix(3,:),wmLocs));
    perceptNoCueProbeWmRT(nSub) = mean(temp_rt); 
    temp_acc = cleanedPerceptACC(cleanedPerceptMatrix(2,:)==1 & cleanedPerceptMatrix(1,:)==5 & ismember(cleanedPerceptMatrix(3,:),wmLocs));
    perceptNoCueProbeWmACC(nSub) = mean(temp_acc);
    perceptNoCueProbeWmRT_pooled = [perceptNoCueProbeWmRT_pooled; temp_rt'];
    
    temp_rt = cleanedPerceptRT(cleanedPerceptMatrix(2,:)==1 & cleanedPerceptMatrix(1,:)==5 & ismember(cleanedPerceptMatrix(3,:),ltmLocs));
    perceptNoCueProbeLtmRT(nSub) = mean(temp_rt); 
    temp_acc = cleanedPerceptACC(cleanedPerceptMatrix(2,:)==1 & cleanedPerceptMatrix(1,:)==5 & ismember(cleanedPerceptMatrix(3,:),ltmLocs));
    perceptNoCueProbeLtmACC(nSub) = mean(temp_acc);
    perceptNoCueProbeLtmRT_pooled = [perceptNoCueProbeLtmRT_pooled; temp_rt'];
    
    temp_rt = cleanedPerceptRT(cleanedPerceptMatrix(2,:)==1 & ((cleanedPerceptMatrix(1,:)==wmLocs(1) & cleanedPerceptMatrix(3,:)==wmLocs(1)) | (cleanedPerceptMatrix(1,:)==wmLocs(2) & cleanedPerceptMatrix(3,:)==wmLocs(2))));
    perceptCueWmProbeWmRT(nSub) = mean(temp_rt); 
    temp_acc = cleanedPerceptACC(cleanedPerceptMatrix(2,:)==1 & ((cleanedPerceptMatrix(1,:)==wmLocs(1) & cleanedPerceptMatrix(3,:)==wmLocs(1)) | (cleanedPerceptMatrix(1,:)==wmLocs(2) & cleanedPerceptMatrix(3,:)==wmLocs(2))));
    perceptCueWmProbeWmACC(nSub) = mean(temp_acc);
    perceptCueWmProbeWmRT_pooled = [perceptCueWmProbeWmRT_pooled; temp_rt'];
    
    temp_rt = cleanedPerceptRT(cleanedPerceptMatrix(2,:)==1 & ((cleanedPerceptMatrix(1,:)==ltmLocs(1) & cleanedPerceptMatrix(3,:)==ltmLocs(1)) | (cleanedPerceptMatrix(1,:)==ltmLocs(2) & cleanedPerceptMatrix(3,:)==ltmLocs(2))));
    perceptCueLtmProbeLtmRT(nSub) = mean(temp_rt); 
    temp_acc = cleanedPerceptACC(cleanedPerceptMatrix(2,:)==1 & ((cleanedPerceptMatrix(1,:)==ltmLocs(1) & cleanedPerceptMatrix(3,:)==ltmLocs(1)) | (cleanedPerceptMatrix(1,:)==ltmLocs(2) & cleanedPerceptMatrix(3,:)==ltmLocs(2))));
    perceptCueLtmProbeLtmACC(nSub) = mean(temp_acc);
    perceptCueLtmProbeLtmRT_pooled = [perceptCueLtmProbeLtmRT_pooled; temp_rt'];
    
    % cue validity effect
    temp_rt = cleanedPerceptRT(cleanedPerceptMatrix(2,:)==1 & ((cleanedPerceptMatrix(1,:)==ltmLocs(2) & cleanedPerceptMatrix(3,:)==ltmLocs(1)) | (cleanedPerceptMatrix(1,:)==ltmLocs(1) & cleanedPerceptMatrix(3,:)==ltmLocs(2))));
    perceptProbeLtmCueTheOtherLtmRT(nSub) = mean(temp_rt); 
    temp_acc = cleanedPerceptACC(cleanedPerceptMatrix(2,:)==1 & ((cleanedPerceptMatrix(1,:)==ltmLocs(2) & cleanedPerceptMatrix(3,:)==ltmLocs(1)) | (cleanedPerceptMatrix(1,:)==ltmLocs(1) & cleanedPerceptMatrix(3,:)==ltmLocs(2))));
    perceptProbeLtmCueTheOtherLtmACC(nSub) = mean(temp_acc);
    
    temp_rt = cleanedPerceptRT(cleanedPerceptMatrix(2,:)==1 & ((ismember(cleanedPerceptMatrix(1,:),wmLocs) & cleanedPerceptMatrix(3,:)==ltmLocs(1)) | (ismember(cleanedPerceptMatrix(1,:),wmLocs) & cleanedPerceptMatrix(3,:)==ltmLocs(2))));
    perceptProbeLtmCueWmRT(nSub) = mean(temp_rt); 
    temp_acc = cleanedPerceptACC(cleanedPerceptMatrix(2,:)==1 & ((ismember(cleanedPerceptMatrix(1,:),wmLocs) & cleanedPerceptMatrix(3,:)==ltmLocs(1)) | (ismember(cleanedPerceptMatrix(1,:),wmLocs) & cleanedPerceptMatrix(3,:)==ltmLocs(2))));
    perceptProbeLtmCueWmACC(nSub) = mean(temp_acc); 
    
    temp_rt = cleanedPerceptRT(cleanedPerceptMatrix(2,:)==1 & ((cleanedPerceptMatrix(1,:)==wmLocs(2) & cleanedPerceptMatrix(3,:)==wmLocs(1)) | (cleanedPerceptMatrix(1,:)==wmLocs(1) & cleanedPerceptMatrix(3,:)==wmLocs(2))));
    perceptProbeWmCueTheOtherWmRT(nSub) = mean(temp_rt); 
    temp_acc = cleanedPerceptACC(cleanedPerceptMatrix(2,:)==1 & ((cleanedPerceptMatrix(1,:)==wmLocs(2) & cleanedPerceptMatrix(3,:)==wmLocs(1)) | (cleanedPerceptMatrix(1,:)==wmLocs(1) & cleanedPerceptMatrix(3,:)==wmLocs(2))));
    perceptProbeWmCueTheOtherWmACC(nSub) = mean(temp_acc);
    
    temp_rt = cleanedPerceptRT(cleanedPerceptMatrix(2,:)==1 & ((ismember(cleanedPerceptMatrix(1,:),ltmLocs) & cleanedPerceptMatrix(3,:)==wmLocs(1)) | (ismember(cleanedPerceptMatrix(1,:),ltmLocs) & cleanedPerceptMatrix(3,:)==wmLocs(2))));
    perceptProbeWmCueLtmRT(nSub) = mean(temp_rt); 
    temp_acc = cleanedPerceptACC(cleanedPerceptMatrix(2,:)==1 & ((ismember(cleanedPerceptMatrix(1,:),ltmLocs) & cleanedPerceptMatrix(3,:)==wmLocs(1)) | (ismember(cleanedPerceptMatrix(1,:),ltmLocs) & cleanedPerceptMatrix(3,:)==wmLocs(2))));
    perceptProbeWmCueLtmACC(nSub) = mean(temp_acc);     
    
    perceptCueWmBenefit(nSub) =  perceptCueWmProbeWmACC(nSub) - perceptNoCueProbeWmACC(nSub);
    perceptCueLtmBenefit(nSub) =  perceptCueLtmProbeLtmACC(nSub) - perceptNoCueProbeLtmACC(nSub);
    
    % 29 March New
    temp_rt = cleanedPerceptRT(cleanedPerceptMatrix(2,:)==1 & ((cleanedPerceptMatrix(3,:)==ltmLocs(1) & cleanedPerceptMatrix(1,:)~=ltmLocs(1) & cleanedPerceptMatrix(1,:)~=5) | (cleanedPerceptMatrix(3,:)==ltmLocs(2) & cleanedPerceptMatrix(1,:)~=ltmLocs(2) & cleanedPerceptMatrix(1,:)~=5)));
    perceptProbeLtmCueInvalidRT(nSub) = mean(temp_rt); 
    temp_acc = cleanedPerceptACC(cleanedPerceptMatrix(2,:)==1 & ((cleanedPerceptMatrix(3,:)==ltmLocs(1) & cleanedPerceptMatrix(1,:)~=ltmLocs(1) & cleanedPerceptMatrix(1,:)~=5) | (cleanedPerceptMatrix(3,:)==ltmLocs(2) & cleanedPerceptMatrix(1,:)~=ltmLocs(2) & cleanedPerceptMatrix(1,:)~=5)));
    perceptProbeLtmCueInvalidACC(nSub) = mean(temp_acc);   

    temp_rt = cleanedPerceptRT(cleanedPerceptMatrix(2,:)==1 & ((cleanedPerceptMatrix(3,:)==wmLocs(1) & cleanedPerceptMatrix(1,:)~=wmLocs(1) & cleanedPerceptMatrix(1,:)~=5) | (cleanedPerceptMatrix(3,:)==wmLocs(2) & cleanedPerceptMatrix(1,:)~=wmLocs(2) & cleanedPerceptMatrix(1,:)~=5)));
    perceptProbeWmCueInvalidRT(nSub) = mean(temp_rt); 
    temp_acc = cleanedPerceptACC(cleanedPerceptMatrix(2,:)==1 & ((cleanedPerceptMatrix(3,:)==wmLocs(1) & cleanedPerceptMatrix(1,:)~=wmLocs(1) & cleanedPerceptMatrix(1,:)~=5) | (cleanedPerceptMatrix(3,:)==wmLocs(2) & cleanedPerceptMatrix(1,:)~=wmLocs(2) & cleanedPerceptMatrix(1,:)~=5)));
    perceptProbeWmCueInvalidACC(nSub) = mean(temp_acc);

end

%% subject info
femaleIndex = find(subGender==1);
numFemale = sum(subGender(femaleIndex));
meanAge = mean(subAge);
stdAge = std(subAge);
disp(['numFemale:' num2str(numFemale) ' meanAge:' num2str(meanAge) ' stdAge:' num2str(stdAge)]);
disp(['meanExcludedPercentage:' num2str(mean(allExcludedPercentages)) ' stdExcludedPercentage:' num2str(std(allExcludedPercentages))]);

if nSubjects ~= 1

    memoryAccNorm3 = Normalization([memoryCueWmACC' memoryCueLtmACC' memoryNoCueACC']);
    memoryAccNorm4 = Normalization([memoryNoCueProbeWmACC' memoryCueWmACC' memoryNoCueProbeLtmACC' memoryCueLtmACC']);
    
    memoryRtNorm3 = Normalization([memoryCueWmRT' memoryCueLtmRT' memoryNoCueRT']);
    memoryRtNorm4 = Normalization([memoryNoCueProbeWmRT' memoryCueWmRT' memoryNoCueProbeLtmRT' memoryCueLtmRT']);
   
    perceptAccNorm = Normalization([perceptNoCueProbeWmACC' perceptCueWmProbeWmACC' perceptNoCueProbeLtmACC' perceptCueLtmProbeLtmACC']);
    perceptRtNorm = Normalization([perceptNoCueProbeWmRT' perceptCueWmProbeWmRT' perceptNoCueProbeLtmRT' perceptCueLtmProbeLtmRT']);
    
    congruencyProbeLtmAccNorm = Normalization([perceptNoCueProbeLtmACC' perceptCueLtmProbeLtmACC' perceptProbeLtmCueTheOtherLtmACC' perceptProbeLtmCueWmACC']);
    congruencyProbeLtmRtNorm = Normalization([perceptNoCueProbeLtmRT' perceptCueLtmProbeLtmRT' perceptProbeLtmCueTheOtherLtmRT' perceptProbeLtmCueWmRT']);
    congruencyProbeWmAccNorm = Normalization([perceptNoCueProbeWmACC' perceptCueWmProbeWmACC' perceptProbeWmCueTheOtherWmACC' perceptProbeWmCueLtmACC']);
    congruencyProbeWmRtNorm = Normalization([perceptNoCueProbeWmRT' perceptCueWmProbeWmRT' perceptProbeWmCueTheOtherWmRT' perceptProbeWmCueLtmRT']);
    % newly added:
    mean_perceptCueWmBenefit = mean(perceptCueWmBenefit);
    mean_perceptCueLtmBenefit = mean(perceptCueLtmBenefit);
    se_perceptCueWmBenefit = std(perceptCueWmBenefit)/length(perceptCueWmBenefit);
    se_perceptCueLtmBenefit = std(perceptCueLtmBenefit)/length(perceptCueLtmBenefit);
end

    % 29 March New
    ValidtyAccNorm = Normalization([perceptNoCueProbeWmACC' perceptProbeWmCueInvalidACC' perceptNoCueProbeLtmACC' perceptProbeLtmCueInvalidACC']);
    ValidtyRtNorm = Normalization([perceptNoCueProbeWmRT' perceptProbeWmCueInvalidRT' perceptNoCueProbeLtmRT' perceptProbeLtmCueInvalidRT']);
%     ValidtyProbeWmAccNorm = Normalization([perceptNoCueProbeWmACC' perceptCueWmProbeWmACC' perceptProbeWmCueInvalidACC']);
%     ValidtyProbeWmRtNorm = Normalization([perceptNoCueProbeWmRT' perceptCueWmProbeWmRT' perceptProbeWmCueInvalidRT']);    
%% Congruency RT

% LTM
f = figure;
f.WindowState = 'maximized';
figureStartup;
sgtitle('Experiment 2','FontSize',25,'FontWeight','bold');
subplot(1,2,1);
hold on
b = bar([mean(perceptNoCueProbeLtmRT),mean(perceptCueLtmProbeLtmRT), mean(perceptProbeLtmCueTheOtherLtmRT), mean(perceptProbeLtmCueWmRT)], 'EdgeColor', 'none', 'BarWidth', 0.8,'FaceColor',colors(2,:));%0.6
hold on
h = ploterr(b.XEndPoints , b.YEndPoints, [],congruencyProbeLtmRtNorm, 'k.', 'abshhxy', 0);
set(h(1), 'marker', 'none','LineWidth',2);
hold on
set(gca,'XTick',1:4);
set(gca,'XTickLabel',{'Neutral','Congruent',sprintf('Incongruent LTM'),sprintf('Incongruent WM')}); 
ylabel('RT (s)');
ylim([0.8 1.1]);
xlim([0.5 4.5]);
title('Probe LTM locations');
xtickangle(45);
% WM

subplot(1,2,2);
hold on
b = bar([mean(perceptNoCueProbeWmRT),mean(perceptCueWmProbeWmRT), mean(perceptProbeWmCueTheOtherWmRT), mean(perceptProbeWmCueLtmRT)], 'EdgeColor', 'none', 'BarWidth', 0.8,'FaceColor',colors(1,:));%0.6
hold on
h = ploterr(b.XEndPoints , b.YEndPoints, [],congruencyProbeWmRtNorm, 'k.', 'abshhxy', 0);
set(h(1), 'marker', 'none','LineWidth',2);
hold on
set(gca,'XTick',1:4);
set(gca,'XTickLabel',{'Neutral','Congruent',sprintf('Incongruent WM'),sprintf('Incongruent LTM')}); 
ylabel('RT (s)');
ylim([0.8 1.1]);
xlim([0.5 4.5]);
xtickangle(45);
title('Probe WM locations');

%% Congruency ACC
f = figure;
f.WindowState = 'maximized';
figureStartup;
sgtitle('Experiment 2','FontSize',25,'FontWeight','bold');
subplot(1,2,1);
hold on
b = bar(100*[mean(perceptNoCueProbeLtmACC),mean(perceptCueLtmProbeLtmACC), mean(perceptProbeLtmCueTheOtherLtmACC), mean(perceptProbeLtmCueWmACC)], 'EdgeColor', 'none', 'BarWidth', 0.8, 'FaceColor',colors(2,:));%0.6
hold on
h = ploterr(b.XEndPoints , b.YEndPoints, [],100*congruencyProbeLtmAccNorm, 'k.', 'abshhxy', 0);
set(h(1), 'marker', 'none','LineWidth',2);
hold on
set(gca,'XTick',1:4);
tickLabels = {'Neutral' 'Congruent' 'Incongruent LTM' 'Incongruent WM'};
set(gca,'XTickLabel',tickLabels);
ylabel('Accuracy (%)');
% xlabel('Retrocue category');
ylim([45 75]);
xlim([0.5 4.5]);
title('Probe LTM locations');
xtickangle(45);

subplot(1,2,2);
hold on
b = bar(100*[mean(perceptNoCueProbeWmACC),mean(perceptCueWmProbeWmACC), mean(perceptProbeWmCueTheOtherWmACC), mean(perceptProbeWmCueLtmACC)], 'EdgeColor', 'none', 'BarWidth', 0.8,'FaceColor',colors(1,:));%0.6
hold on
h = ploterr(b.XEndPoints , b.YEndPoints, [],100*congruencyProbeWmAccNorm, 'k.', 'abshhxy', 0);
set(h(1), 'marker', 'none','LineWidth',2);
hold on
set(gca,'XTick',1:4);
tickLabels = {'Neutral' 'Congruent' 'Incongruent WM' 'Incongruent LTM'};
set(gca,'XTickLabel',tickLabels);
% xlabel('Retrocue category');
ylabel('Accuracy (%)');
ylim([45 75]);
xlim([0.5 4.5]);
title('Probe WM locations');
xtickangle(45);

%%
if nSubjects ~= 1
    [anovastats_congruency_ltm_acc,anovatable_congruency_ltm_acc]=mes1way([perceptNoCueProbeLtmACC' perceptCueLtmProbeLtmACC' perceptProbeLtmCueTheOtherLtmACC' perceptProbeLtmCueWmACC'],{'partialeta2','g_psi','psibysd'},'isDep',1,'tDenom','sd','cWeight',[-1 1 0 0; 2 0 -1 -1; 1 0 -1 0; 1 0 0 -1; 0 1 -1 0; 0 1 0 -1]);
    [anovastats_congruency_wm_acc,anovatable_congruency_wm_acc]=mes1way([perceptNoCueProbeWmACC' perceptCueWmProbeWmACC' perceptProbeWmCueTheOtherWmACC' perceptProbeWmCueLtmACC'],{'partialeta2','g_psi','psibysd'},'isDep',1,'tDenom','sd','cWeight',[-1 1 0 0; 2 0 -1 -1; 1 0 -1 0; 1 0 0 -1; 0 1 -1 0; 0 1 0 -1]);
    [anovastats_congruency_ltm_rt,anovatable_congruency_ltm_rt]=mes1way([perceptNoCueProbeLtmRT' perceptCueLtmProbeLtmRT' perceptProbeLtmCueTheOtherLtmRT' perceptProbeLtmCueWmRT'],{'partialeta2','g_psi','psibysd'},'isDep',1,'tDenom','sd','cWeight',[-1 1 0 0; 2 0 -1 -1; 1 0 -1 0; 1 0 0 -1; 0 1 -1 0; 0 1 0 -1]);
    [anovastats_congruency_wm_rt,anovatable_congruency_wm_rt]=mes1way([perceptNoCueProbeWmRT' perceptCueWmProbeWmRT' perceptProbeWmCueTheOtherWmRT' perceptProbeWmCueLtmRT'],{'partialeta2','g_psi','psibysd'},'isDep',1,'tDenom','sd','cWeight',[-1 1 0 0; 2 0 -1 -1; 1 0 -1 0; 1 0 0 -1; 0 1 -1 0; 0 1 0 -1]);
end

%% NEW Congruency plots
%% Perception RT
f = figure;
figureStartup;
subplot(2,2,1);
f.WindowState = 'maximized';
hold on

b11 = bar([NaN,mean(perceptProbeWmCueInvalidRT); NaN,NaN], 'BarWidth', 0.8);
hatchfill2(b11(2),'cross','HatchAngle',45,'hatchcolor',colors(1,:),'HatchLineWidth',1);
b11(2).FaceColor = 'none';
b11(2).EdgeColor = colors(1,:);
b11(2).LineWidth = 1;
hold on

b22 = bar([NaN,NaN; NaN,mean(perceptProbeLtmCueInvalidRT)], 'BarWidth', 0.8);
hatchfill2(b22(2),'cross','HatchAngle',45,'hatchcolor',colors(2,:),'HatchLineWidth',1,'HatchOffset',0.5);
b22(2).FaceColor = 'none';
b22(2).EdgeColor = colors(2,:);
b22(2).LineWidth = 1;
hold on
bb1 = bar([mean(perceptNoCueProbeWmRT),NaN; NaN, NaN],'BarWidth', 0.8);
bb1(1).EdgeColor = colors(1,:);
bb1(1).FaceColor = [1 1 1];
bb1(1).LineWidth = 1;
hold on
bb2 = bar([NaN,NaN; mean(perceptNoCueProbeLtmRT),NaN],'BarWidth', 0.8);
bb2(1).EdgeColor = colors(2,:);
bb2(1).FaceColor = [1 1 1];
bb2(1).LineWidth = 1;
hold on

% for i = 1:nSubjects
%     % WM
%     plot([bb1(1).XEndPoints(1) b1(2).XEndPoints(1)],[perceptNoCueProbeWmRT(i) perceptCueWmProbeWmRT(i)],'-','Color',extra_colors(10,:));
% %     scatter([bb1(1).XEndPoints(1) b1(2).XEndPoints(1)],[perceptNoCueProbeWmRT(i) perceptCueWmProbeWmRT(i)],[],extra_colors(10,:),'filled');
%     hold on
%     % LTM
%     plot([bb2(1).XEndPoints(2) b2(2).XEndPoints(2)],[perceptNoCueProbeLtmRT(i) perceptCueLtmProbeLtmRT(i)],'-','Color',extra_colors(10,:));
% %     scatter([bb2(1).XEndPoints(2) b2(2).XEndPoints(2)],[perceptNoCueProbeLtmRT(i) perceptCueLtmProbeLtmRT(i)],[],extra_colors(10,:),'filled');
%     hold on
% end


h = ploterr([bb1(1).XEndPoints(1) b11(2).XEndPoints(1) bb2(1).XEndPoints(2) b22(2).XEndPoints(2)] , [bb1(1).YEndPoints(1) b11(2).YEndPoints(1) bb2(1).YEndPoints(2) b22(2).YEndPoints(2)], [],ValidtyRtNorm, 'k.', 'abshhxy', 0);
set(h(1), 'marker', 'none','LineWidth',2);
hold on
set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'WM','LTM'});
xlim([0.5 2.5]);
ylim([0.8 1.2]);
ylabel('RT (s)');

CC = [perceptNoCueProbeWmRT' perceptCueWmProbeWmRT' perceptProbeWmCueInvalidRT' perceptNoCueProbeLtmRT'  perceptCueLtmProbeLtmRT' perceptProbeLtmCueInvalidRT'];
Ncompare = nchoosek(1:6,2);
for n = 1:size(Ncompare,1)
    tstats{n} = mes(CC(:,Ncompare(n,2)),CC(:,Ncompare(n,1)),'hedgesg','isDep',1);
    p_temp(n) = tstats{n}.t.p;
end
% groups = mat2cell(nchoosek([b(1).XEndPoints(1) b(2).XEndPoints(1) b(1).XEndPoints(2) b(2).XEndPoints(2)],2),ones(1,size(Ncompare,1)));
% H=sigstar(groups,p_temp*2,1);

list_comparisons = [1,6]; % 1,2,5,6  are the interested comparisons
for i = list_comparisons
%     mysigstar(gca,groups{i},myRange(ylim)*0.05+findMinY(groups{i}),p_temp(i));
end

% legend(b,{'LTM location','WM location'},'Location','eastoutside');
% legend boxoff
% saveas(f,[figDir identifier  '-perception-rt.emf']);
% saveas(gca,[figDir identifier  '-2-2.svg']);
%%
if nSubjects ~= 1   
    %first factor: WM/LTM; second factor: not cued / cued
    groups=[repmat([1,1],nSubjects,1);
        repmat([1,2],nSubjects,1);
        repmat([2,1],nSubjects,1);
        repmat([2,2],nSubjects,1)];
    anova_twoway_data = [perceptNoCueProbeWmRT';
        perceptProbeWmCueInvalidRT';
        perceptNoCueProbeLtmRT';
        perceptProbeLtmCueInvalidRT'];
    [anovastats_twoway_validty_rt,anovatable_twoway_validty_rt]=mes2way(anova_twoway_data,groups,'partialeta2','isDep',[1 1]);
    ttest_wm_validty_rt_1 = mes(perceptNoCueProbeWmRT',perceptProbeWmCueInvalidRT','hedgesg','isDep',1);
    ttest_wm_validty_rt_2 = mes(perceptProbeWmCueInvalidRT',perceptCueWmProbeWmRT','hedgesg','isDep',1);
    ttest_ltm_validty_rt_1 = mes(perceptNoCueProbeLtmRT',perceptProbeLtmCueInvalidRT','hedgesg','isDep',1);
    ttest_ltm_validty_rt_2 = mes(perceptProbeLtmCueInvalidRT',perceptCueLtmProbeLtmRT','hedgesg','isDep',1);
end


% f = figure;
% figureStartup;
subplot(2,2,2);
f.WindowState = 'maximized';
hold on

b11 = bar(100*[NaN,mean(perceptProbeWmCueInvalidACC); NaN, NaN], 'BarWidth', 0.8);
hatchfill2(b11(2),'cross','HatchAngle',45,'hatchcolor',colors(1,:),'HatchLineWidth',1);
b11(2).FaceColor = 'none';
b11(2).EdgeColor = colors(1,:);
b11(2).LineWidth = 1;
hold on

b22 = bar(100*[NaN,NaN; NaN,mean(perceptProbeLtmCueInvalidACC)], 'BarWidth', 0.8);
hatchfill2(b22(2),'cross','HatchAngle',45,'hatchcolor',colors(2,:),'HatchLineWidth',1,'HatchOffset',0.5);
b22(2).FaceColor = 'none';
b22(2).EdgeColor = colors(2,:);
b22(2).LineWidth = 1;
hold on
bb1 = bar(100*[mean(perceptNoCueProbeWmACC),NaN; NaN, NaN],'BarWidth', 0.8);
bb1(1).EdgeColor = colors(1,:);
bb1(1).FaceColor = [1 1 1];
bb1(1).LineWidth = 1;
hold on
bb2 = bar(100*[NaN,NaN; mean(perceptNoCueProbeLtmACC),NaN],'BarWidth', 0.8);
bb2(1).EdgeColor = colors(2,:);
bb2(1).FaceColor = [1 1 1];
bb2(1).LineWidth = 1;
hold on

% for i = 1:nSubjects
%     % WM
%     plot([bb1(1).XEndPoints(1) b1(2).XEndPoints(1)],[perceptNoCueProbeWmRT(i) perceptCueWmProbeWmRT(i)],'-','Color',extra_colors(10,:));
% %     scatter([bb1(1).XEndPoints(1) b1(2).XEndPoints(1)],[perceptNoCueProbeWmRT(i) perceptCueWmProbeWmRT(i)],[],extra_colors(10,:),'filled');
%     hold on
%     % LTM
%     plot([bb2(1).XEndPoints(2) b2(2).XEndPoints(2)],[perceptNoCueProbeLtmRT(i) perceptCueLtmProbeLtmRT(i)],'-','Color',extra_colors(10,:));
% %     scatter([bb2(1).XEndPoints(2) b2(2).XEndPoints(2)],[perceptNoCueProbeLtmRT(i) perceptCueLtmProbeLtmRT(i)],[],extra_colors(10,:),'filled');
%     hold on
% end


h = ploterr([bb1(1).XEndPoints(1) b11(2).XEndPoints(1) bb2(1).XEndPoints(2) b22(2).XEndPoints(2)] , [bb1(1).YEndPoints(1) b11(2).YEndPoints(1) bb2(1).YEndPoints(2) b22(2).YEndPoints(2)], [],100*ValidtyAccNorm, 'k.', 'abshhxy', 0);
set(h(1), 'marker', 'none','LineWidth',2);
hold on
set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'WM','LTM'});
xlim([0.5 2.5]);
ylim([45 75]);
ylabel('Accuracy (%)');

CC = [perceptNoCueProbeWmACC' perceptCueWmProbeWmACC' perceptProbeWmCueInvalidACC' perceptNoCueProbeLtmACC'  perceptCueLtmProbeLtmACC' perceptProbeLtmCueInvalidACC'];
Ncompare = nchoosek(1:6,2);
for n = 1:size(Ncompare,1)
    tstats{n} = mes(CC(:,Ncompare(n,2)),CC(:,Ncompare(n,1)),'hedgesg','isDep',1);
    p_temp(n) = tstats{n}.t.p;
end
% groups = mat2cell(nchoosek([b(1).XEndPoints(1) b(2).XEndPoints(1) b(1).XEndPoints(2) b(2).XEndPoints(2)],2),ones(1,size(Ncompare,1)));
% H=sigstar(groups,p_temp*2,1);

list_comparisons = [1,6]; % 1,2,5,6  are the interested comparisons
for i = list_comparisons
%     mysigstar(gca,groups{i},myRange(ylim)*0.05+findMinY(groups{i}),p_temp(i));
end

% legend(b,{'LTM location','WM location'},'Location','eastoutside');
% legend boxoff
% saveas(f,[figDir identifier  '-perception-rt.emf']);
% saveas(gca,[figDir identifier  '-2-2.svg']);
%%
if nSubjects ~= 1   
    %first factor: WM/LTM; second factor: not cued / cued
    groups=[repmat([1,1],nSubjects,1);
        repmat([1,2],nSubjects,1);
        repmat([2,1],nSubjects,1);
        repmat([2,2],nSubjects,1)];
    anova_twoway_data = [perceptNoCueProbeWmACC';
        perceptProbeWmCueInvalidACC';
        perceptNoCueProbeLtmACC';
        perceptProbeLtmCueInvalidACC'];
    [anovastats_twoway_validty_acc,anovatable_twoway_validty_acc]=mes2way(anova_twoway_data,groups,'partialeta2','isDep',[1 1]);
    ttest_wm_validty_acc_1 = mes(perceptNoCueProbeWmACC',perceptProbeWmCueInvalidACC','hedgesg','isDep',1);
    ttest_wm_validty_acc_2 = mes(perceptProbeWmCueInvalidACC',perceptCueWmProbeWmACC','hedgesg','isDep',1);
%     ttest_wm_validty_acc_3 = mes(perceptNoCueProbeWmACC',perceptCueWmProbeWmACC','hedgesg','isDep',1);
    ttest_ltm_validty_acc_1 = mes(perceptNoCueProbeLtmACC',perceptProbeLtmCueInvalidACC','hedgesg','isDep',1);
    ttest_ltm_validty_acc_2 = mes(perceptProbeLtmCueInvalidACC',perceptCueLtmProbeLtmACC','hedgesg','isDep',1);
%     ttest_ltm_validty_acc_3 = mes(perceptNoCueProbeLtmACC',perceptCueLtmProbeLtmACC','hedgesg','isDep',1);
    
end


%% Memory RT

f = figure;
figureStartup;
subplot(2,2,1);
f.WindowState = 'maximized';

hold on
ylim([2 6]);
% axis tight;
xlim([0.5 2.5]);
b = bar([mean(memoryNoCueProbeWmRT),mean(memoryCueWmRT); mean(memoryNoCueProbeLtmRT), mean(memoryCueLtmRT)], 'EdgeColor', 'none', 'BarWidth', 0.8);%0.6
b(2).FaceColor = 'flat';
b(2).CData = [colors(1,:);colors(2,:)];
hold on
bb1 = bar([mean(memoryNoCueProbeWmRT),NaN; NaN, NaN],'BarWidth', 0.8);
bb1(1).EdgeColor = colors(1,:);
bb1(1).FaceColor = [1 1 1];
bb1(1).LineWidth = 1;
hold on
bb2 = bar([NaN,NaN; mean(memoryNoCueProbeLtmRT),NaN],'BarWidth', 0.8);
bb2(1).EdgeColor = colors(2,:);
bb2(1).FaceColor = [1 1 1];
bb2(1).LineWidth = 1;
hold on

for i = 1:nSubjects
    % WM
    plot([b(1).XEndPoints(1) b(2).XEndPoints(1)],[memoryNoCueProbeWmRT(i) memoryCueWmRT(i)],'-','Color',extra_colors(10,:));
%     scatter([b(1).XEndPoints(1) b(2).XEndPoints(1)],[memoryNoCueProbeWmRT(i) memoryCueWmRT(i)],[],extra_colors(10,:),'filled');
    hold on
    % LTM
    plot([b(1).XEndPoints(2) b(2).XEndPoints(2)],[memoryNoCueProbeLtmRT(i) memoryCueLtmRT(i)],'-','Color',extra_colors(10,:));
%     scatter([b(1).XEndPoints(2) b(2).XEndPoints(2)],[memoryNoCueProbeLtmRT(i) memoryCueLtmRT(i)],[],extra_colors(10,:),'filled');
    hold on
end

h = ploterr([b(1).XEndPoints b(2).XEndPoints] , [b(1).YEndPoints b(2).YEndPoints], [],memoryRtNorm4([1 3 2 4]), 'k.', 'abshhxy', 0);
set(h(1), 'marker', 'none','LineWidth',2);
hold on
set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'WM','LTM'}); 
ylabel('RT (s)');

CC = [memoryNoCueProbeWmRT' memoryCueWmRT' memoryNoCueProbeLtmRT' memoryCueLtmRT'];
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
%%
if nSubjects ~= 1   
    %first factor: WM/LTM; second factor: not cued / cued
    groups=[repmat([1,1],nSubjects,1);
        repmat([1,2],nSubjects,1);
        repmat([2,1],nSubjects,1);
        repmat([2,2],nSubjects,1)];
    anova_twoway_data = [memoryNoCueProbeWmRT';
        memoryCueWmRT';
        memoryNoCueProbeLtmRT';
        memoryCueLtmRT'];
    [anovastats_twoway_memo_rt,anovatable_twoway_memo_rt]=mes2way(anova_twoway_data,groups,'partialeta2','isDep',[1 1]);

    ttest_wm_memo_rt = mes(memoryNoCueProbeWmRT',memoryCueWmRT','hedgesg','isDep',1);
    ttest_ltm_memo_rt = mes(memoryNoCueProbeLtmRT',memoryCueLtmRT','hedgesg','isDep',1);
    ttest_nocue_memo_rt = mes(memoryNoCueProbeWmRT',memoryNoCueProbeLtmRT','hedgesg','isDep',1);
    ttest_cued_memo_rt = mes(memoryCueWmRT',memoryCueLtmRT','hedgesg','isDep',1);
    ttest_memo_rt_benefit = mes(memoryNoCueProbeWmRT'-memoryCueWmRT',memoryNoCueProbeLtmRT'-memoryCueLtmRT','hedgesg','isDep',1);
end
%% Memory ACC
% f = figure;
% figureStartup;
subplot(2,2,2);
f.WindowState = 'maximized';
hold on
ylim([10 100]);
% axis tight;
xlim([0.5 2.5]);
b = bar([mean(memoryNoCueProbeWmACC),mean(memoryCueWmACC); mean(memoryNoCueProbeLtmACC), mean(memoryCueLtmACC)], 'EdgeColor', 'none', 'BarWidth', 0.8);%0.6
b(2).FaceColor = 'flat';
b(2).CData = [colors(1,:);colors(2,:)];
hold on
bb1 = bar([mean(memoryNoCueProbeWmACC),NaN; NaN, NaN],'BarWidth', 0.8);
bb1(1).EdgeColor = colors(1,:);
bb1(1).FaceColor = [1 1 1];
bb1(1).LineWidth = 1;
hold on
bb2 = bar([NaN,NaN; mean(memoryNoCueProbeLtmACC),NaN],'BarWidth', 0.8);
bb2(1).EdgeColor = colors(2,:);
bb2(1).FaceColor = [1 1 1];
bb2(1).LineWidth = 1;
hold on

for i = 1:nSubjects
    % WM
    plot([b(1).XEndPoints(1) b(2).XEndPoints(1)],[memoryNoCueProbeWmACC(i) memoryCueWmACC(i)],'-','Color',extra_colors(10,:));
%     scatter([b(1).XEndPoints(1) b(2).XEndPoints(1)],[memoryNoCueProbeWmACC(i) memoryCueWmACC(i)],[],extra_colors(10,:),'filled');
    hold on
    % LTM
    plot([b(1).XEndPoints(2) b(2).XEndPoints(2)],[memoryNoCueProbeLtmACC(i) memoryCueLtmACC(i)],'-','Color',extra_colors(10,:));
%     scatter([b(1).XEndPoints(2) b(2).XEndPoints(2)],[memoryNoCueProbeLtmACC(i) memoryCueLtmACC(i)],[],extra_colors(10,:),'filled');
    hold on
end
h = ploterr([b(1).XEndPoints b(2).XEndPoints] , [b(1).YEndPoints b(2).YEndPoints], [],memoryAccNorm4([1 3 2 4]), 'k.', 'abshhxy', 0);
set(h(1), 'marker', 'none','LineWidth',2);
hold on   
set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'WM','LTM'}); 
ylabel('Reproduction error (°)');

CC = [memoryNoCueProbeWmACC' memoryCueWmACC' memoryNoCueProbeLtmACC' memoryCueLtmACC'];
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
% saveas(f,[figDir identifier  '-memory-acc.emf']);

%%
if nSubjects ~= 1   
    %first factor: WM/LTM; second factor: not cued / cued
    groups=[repmat([1,1],nSubjects,1);
        repmat([1,2],nSubjects,1);
        repmat([2,1],nSubjects,1);
        repmat([2,2],nSubjects,1)];
    anova_twoway_data = [memoryNoCueProbeWmACC';
        memoryCueWmACC';
        memoryNoCueProbeLtmACC';
        memoryCueLtmACC'];
    [anovastats_twoway_memo_acc,anovatable_twoway_memo_acc]=mes2way(anova_twoway_data,groups,'partialeta2','isDep',[1 1]);
    ttest_wm_memo_acc = mes(memoryNoCueProbeWmACC',memoryCueWmACC','hedgesg','isDep',1);
    ttest_ltm_memo_acc = mes(memoryNoCueProbeLtmACC',memoryCueLtmACC','hedgesg','isDep',1);
    ttest_nocue_memo_acc = mes(memoryNoCueProbeWmACC',memoryNoCueProbeLtmACC','hedgesg','isDep',1);
    ttest_cued_memo_acc = mes(memoryCueWmACC',memoryCueLtmACC','hedgesg','isDep',1);
end

%% Perception RT
% f = figure;
% figureStartup;
subplot(2,2,3);
f.WindowState = 'maximized';
hold on
b1 = bar([NaN,mean(perceptCueWmProbeWmRT); NaN, NaN], 'BarWidth', 0.8);
hatchfill2(b1(2),'single','HatchAngle',45,'hatchcolor',colors(1,:),'HatchLineWidth',1);
b1(2).FaceColor = 'none';
b1(2).EdgeColor = colors(1,:);
b1(2).LineWidth = 1;
hold on
b2 = bar([NaN,NaN; NaN, mean(perceptCueLtmProbeLtmRT)], 'BarWidth', 0.8);
hatchfill2(b2(2),'single','HatchAngle',45,'hatchcolor',colors(2,:),'HatchLineWidth',1,'HatchOffset',0.5);
b2(2).FaceColor = 'none';
b2(2).EdgeColor = colors(2,:);
b2(2).LineWidth = 1;
hold on
bb1 = bar([mean(perceptNoCueProbeWmRT),NaN; NaN, NaN],'BarWidth', 0.8);
bb1(1).EdgeColor = colors(1,:);
bb1(1).FaceColor = [1 1 1];
bb1(1).LineWidth = 1;
hold on
bb2 = bar([NaN,NaN; mean(perceptNoCueProbeLtmRT),NaN],'BarWidth', 0.8);
bb2(1).EdgeColor = colors(2,:);
bb2(1).FaceColor = [1 1 1];
bb2(1).LineWidth = 1;
hold on

for i = 1:nSubjects
    % WM
    plot([bb1(1).XEndPoints(1) b1(2).XEndPoints(1)],[perceptNoCueProbeWmRT(i) perceptCueWmProbeWmRT(i)],'-','Color',extra_colors(10,:));
%     scatter([bb1(1).XEndPoints(1) b1(2).XEndPoints(1)],[perceptNoCueProbeWmRT(i) perceptCueWmProbeWmRT(i)],[],extra_colors(10,:),'filled');
    hold on
    % LTM
    plot([bb2(1).XEndPoints(2) b2(2).XEndPoints(2)],[perceptNoCueProbeLtmRT(i) perceptCueLtmProbeLtmRT(i)],'-','Color',extra_colors(10,:));
%     scatter([bb2(1).XEndPoints(2) b2(2).XEndPoints(2)],[perceptNoCueProbeLtmRT(i) perceptCueLtmProbeLtmRT(i)],[],extra_colors(10,:),'filled');
    hold on
end
h = ploterr([bb1(1).XEndPoints(1) b1(2).XEndPoints(1) bb2(1).XEndPoints(2) b2(2).XEndPoints(2)] , [bb1(1).YEndPoints(1) b1(2).YEndPoints(1) bb2(1).YEndPoints(2) b2(2).YEndPoints(2)], [],perceptRtNorm([1 2 3 4]), 'k.', 'abshhxy', 0);
set(h(1), 'marker', 'none','LineWidth',2);
hold on
set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'WM','LTM'});
xlim([0.5 2.5]);
ylim([0 4]);
ylabel('RT (s)');

CC = [perceptNoCueProbeWmRT' perceptCueWmProbeWmRT' perceptNoCueProbeLtmRT'  perceptCueLtmProbeLtmRT'];
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
% saveas(f,[figDir identifier  '-perception-rt.emf']);
saveas(gca,[figDir identifier  '-2-2.svg']);
%%
if nSubjects ~= 1   
    %first factor: WM/LTM; second factor: not cued / cued
    groups=[repmat([1,1],nSubjects,1);
        repmat([1,2],nSubjects,1);
        repmat([2,1],nSubjects,1);
        repmat([2,2],nSubjects,1)];
    anova_twoway_data = [perceptNoCueProbeWmRT';
        perceptCueWmProbeWmRT';
        perceptNoCueProbeLtmRT';
        perceptCueLtmProbeLtmRT'];
    [anovastats_twoway_percept_rt,anovatable_twoway_percept_rt]=mes2way(anova_twoway_data,groups,'partialeta2','isDep',[1 1]);
    ttest_wm_percept_rt = mes(perceptNoCueProbeWmRT',perceptCueWmProbeWmRT','hedgesg','isDep',1);
    ttest_ltm_percept_rt = mes(perceptNoCueProbeLtmRT',perceptCueLtmProbeLtmRT','hedgesg','isDep',1);
end

%% Perception ACC
% f = figure;
% figureStartup;
subplot(2,2,4);
f.WindowState = 'maximized';
hold on
b1 = bar(100*[NaN,mean(perceptCueWmProbeWmACC); NaN, NaN], 'BarWidth', 0.8);
hatchfill2(b1(2),'single','HatchAngle',45,'hatchcolor',colors(1,:),'HatchLineWidth',1);
b1(2).FaceColor = 'none';
b1(2).EdgeColor = colors(1,:);
b1(2).LineWidth = 1;
hold on
b2 = bar(100*[NaN,NaN; NaN, mean(perceptCueLtmProbeLtmACC)], 'BarWidth', 0.8);
hatchfill2(b2(2),'single','HatchAngle',45,'hatchcolor',colors(2,:),'HatchLineWidth',1,'HatchOffset',0.5);
b2(2).FaceColor = 'none';
b2(2).EdgeColor = colors(2,:);
b2(2).LineWidth = 1;
hold on
bb1 = bar(100*[mean(perceptNoCueProbeWmACC),NaN; NaN, NaN],'BarWidth', 0.8);
bb1(1).EdgeColor = colors(1,:);
bb1(1).FaceColor = [1 1 1];
bb1(1).LineWidth = 1;
hold on
bb2 = bar(100*[NaN,NaN; mean(perceptNoCueProbeLtmACC),NaN],'BarWidth', 0.8);
bb2(1).EdgeColor = colors(2,:);
bb2(1).FaceColor = [1 1 1];
bb2(1).LineWidth = 1;
hold on
  
for i = 1:nSubjects
    % WM
    plot([bb1(1).XEndPoints(1) b1(2).XEndPoints(1)],100*[perceptNoCueProbeWmACC(i) perceptCueWmProbeWmACC(i)],'-','Color',extra_colors(10,:));
%     scatter([bb1(1).XEndPoints(1) b1(2).XEndPoints(1)],100*[perceptNoCueProbeWmACC(i) perceptCueWmProbeWmACC(i)],[],extra_colors(10,:),'filled');
    hold on
    % LTM
    plot([bb2(1).XEndPoints(2) b2(2).XEndPoints(2)],100*[perceptNoCueProbeLtmACC(i) perceptCueLtmProbeLtmACC(i)],'-','Color',extra_colors(10,:));
%     scatter([bb2(1).XEndPoints(2) b2(2).XEndPoints(2)],100*[perceptNoCueProbeLtmACC(i) perceptCueLtmProbeLtmACC(i)],[],extra_colors(10,:),'filled');
    hold on
end
h = ploterr([bb1(1).XEndPoints(1) b1(2).XEndPoints(1) bb2(1).XEndPoints(2) b2(2).XEndPoints(2)] , [bb1(1).YEndPoints(1) b1(2).YEndPoints(1) bb2(1).YEndPoints(2) b2(2).YEndPoints(2)], [],100*perceptAccNorm([1 2 3 4]), 'k.', 'abshhxy', 0);
set(h(1), 'marker', 'none','LineWidth',2);
hold on 
set(gca,'XTick',1:2);
set(gca,'XTickLabel',{'WM','LTM'}); 
xlim([0.5 2.5]);
ylim([30 100]);
ylabel('Accuracy (%)');

CC = [perceptNoCueProbeWmACC' perceptCueWmProbeWmACC' perceptNoCueProbeLtmACC'  perceptCueLtmProbeLtmACC'];
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

%%
if nSubjects ~= 1   
    %first factor: WM/LTM; second factor: not cued / cued
    groups=[repmat([1,1],nSubjects,1);
        repmat([1,2],nSubjects,1);
        repmat([2,1],nSubjects,1);
        repmat([2,2],nSubjects,1)];
    anova_twoway_data = [perceptNoCueProbeWmACC';
        perceptCueWmProbeWmACC';
        perceptNoCueProbeLtmACC';
        perceptCueLtmProbeLtmACC'];
    [anovastats_twoway_percept_acc,anovatable_twoway_percept_acc]=mes2way(anova_twoway_data,groups,'partialeta2','isDep',[1 1]);
    ttest_wm_percept_acc = mes(perceptNoCueProbeWmACC,perceptCueWmProbeWmACC,'hedgesg','isDep',1);
    ttest_ltm_percept_acc = mes(perceptNoCueProbeLtmACC,perceptCueLtmProbeLtmACC,'hedgesg','isDep',1);
    ttest_nocue_wm_ltm_acc = mes(perceptNoCueProbeWmACC,perceptNoCueProbeLtmACC,'hedgesg','isDep',1);
    ttest_cued_wm_ltm_acc = mes(perceptCueWmProbeWmACC,perceptCueLtmProbeLtmACC,'hedgesg','isDep',1);
end



%% ACC benefit
f = figure;
figureStartup;
subplot(1,2,1); hold on;
f.WindowState = 'maximized';
% rather than a square plot, make it thinner
handles=violinPlot(100*perceptCueLtmBenefit', 'histOri', 'right', 'widthDiv', [2 2], 'showMM', 4, ...
    'color',  mat2cell(colors(2, : ), 1)); %%%%%%MODIFIED Line 678 in violinplot.m
% handles{1}.EdgeColor = colors(2, : );
hatchfill2(handles{1},'single','HatchAngle',45,'HatchColor',colors(2, : ),'HatchLineWidth',1);
hold on
handles=violinPlot(100*perceptCueWmBenefit', 'histOri', 'left', 'widthDiv', [2 1], 'showMM', 4, ...
    'color',  mat2cell(colors(1, : ), 1));
% handles{1}.EdgeColor = colors(1, : );
hatchfill2(handles{1},'single','HatchAngle',45,'HatchColor',colors(1, : ),'HatchLineWidth',1,'HatchOffset',0.5);
set(gca, 'XTick', [0.6 1.4], 'XTickLabel', {'WM cue', 'LTM cue'});
xlim([0.2 1.8]);
ylim([-40,70]);
ylabel('Accuracy benefit (%)');

hold on
yline(0,'--k','LineWidth',1);

if nSubjects ~= 1
    ttest_accBenefit = mes(perceptCueWmBenefit',perceptCueLtmBenefit','hedgesg','isDep',1);
    ttest_onesample_WM = mes(perceptCueWmBenefit',0,'g1');
    ttest_onesample_LTM = mes(perceptCueLtmBenefit',0,'g1');
end

%add significance stars for each bar
xtick = get(gca, 'xtick');
ypos = -30; % plot below
% mysigstar(gca, xtick(1), ypos, ttest_onesample_LTM.t.p);
% mysigstar(gca, xtick(2), ypos, ttest_onesample_WM.t.p);
% mysigstar(gca, xtick, 65, ttest_accBenefit.t.p);
saveas(gca,[figDir identifier  '-perception-accBenefit.svg']);

%}

%%
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
    
    
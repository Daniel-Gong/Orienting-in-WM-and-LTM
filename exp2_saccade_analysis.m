%% It's always good to start with a clean sheet
clc, clear, close all, warning('off','all')

%% set project name
projectname = 'detect gaze shift (microsaccade) from gaze position';

%% prepare Toolbox
parent_folder = '/Users/dongyu/Documents/Github/spontaneousMicrosaccade_EEGAlpha/';
package_folder = '/Users/dongyu/Documents/Github/spontaneousMicrosaccade_EEGAlpha/customToolbox/';
data_folder = '/Users/dongyu/MATLAB-Drive/exp2_eye_data';

cfg= [];
cfg.package_folder = package_folder;
cfg.futureTB_path = [package_folder 'futureToolBox'];
cd(cfg.futureTB_path);
futureTB_setup(cfg);

%% Set directions
write_dir_eye = creatDir([parent_folder 'eye_results_probe']);
write_dir_fig = creatDir([write_dir_eye filesep 'figures']);

%% get gaze shift
readDir = [write_dir_eye filesep 'epoched_data'];
nonNaNDir = [write_dir_eye filesep 'trial_ok_noNan'];
subjList =  get_subFiles(readDir);
nonNaNList = get_subFiles(nonNaNDir);
output_dir = creatDir([write_dir_eye filesep 'gaze_shift']);
for subjInd = 1:length(subjList)
    load(subjList{subjInd})
    load(nonNaNList{subjInd})
    % get data in x axis
    data_shift =[];
    if strcmp(subjList{subjInd}(end-7:end-4),'pp49') % deal with the missing data (this participant doesn't have the first 80 trials of learning phase)
        startTrial = 1;
        disp('start from trial 1');
    else
        startTrial = 81;
        disp('start from trial 81');
    end
    eye_dataX = [];
    for trial = startTrial:length(eye_data.trial)
        trial_dataX = eye_data.trial{trial}(1,:);
        eye_dataX = [eye_dataX;trial_dataX];
    end
    nonNaNtrials = event.sel(startTrial:length(eye_data.trial));
    eye_dataX = eye_dataX(nonNaNtrials,:);
    % get gaze shift 
    cfg = [];
    cfg.threshold = 3;
    [eye_shift,time_shift] = PBlab_gazepos2shift_1D(cfg, eye_dataX, eye_data.time{1});
    eye_shift = 2*atand(eye_shift/(2*95*28.3604));

    % remove saccade whose size > 1 degree of visual angle
    eye_shift(abs(eye_shift)>1) = 0;

    data_shift.shift = eye_shift;
    data_shift.time = time_shift;
    trialInfo = eye_data.trialinfo(startTrial:length(eye_data.trial));
    data_shift.trialinfo = trialInfo(nonNaNtrials);
    
    % save data
    save([output_dir filesep subjList{subjInd}(end-7:end-4)] ,'data_shift')
end

%% plot gaze magnitude x shift size figure: LTM
sublist_shift = get_subFiles([write_dir_eye filesep 'gaze_shift']);

for subjInd = 1:length(sublist_shift)
    load(sublist_shift{subjInd})
    GA_struct = [];
    cfg = [];
    cfg.size_range = [0.05 5.5];
    cfg.binWin = 0.25;
    cfg.binstep = 0.05;
    subNumber = str2double(sublist_shift{subjInd}(end-5:end-4));
    disp(['subNumber:' num2str(subNumber)]);
    if mod(subNumber,2) == 0
        cfg.trigs_left = 721;
        cfg.trigs_right = 722;
    elseif mod(subNumber,2) == 1
        cfg.trigs_left = 722;
        cfg.trigs_right = 721;
    end
    [rate_size,bin_range]=gazeShiftRateOverSize(cfg, data_shift.shift, data_shift.trialinfo);
    
    GA_struct.toward(subjInd,:,:) = rate_size.toward;
    GA_struct.away(subjInd,:,:) = rate_size.away;
    GA_struct.diff(subjInd,:,:) = rate_size.diff;
end

GA_struct.bin_range = bin_range;
GA_struct.time = data_shift.time;

% save the group data
save([write_dir_eye filesep 'GA_shift_rateAndsize_LTM'] ,'GA_struct')

%% plot gaze magnitude x shift size figure: WM
sublist_shift = get_subFiles([write_dir_eye filesep 'gaze_shift']);

for subjInd = 1:length(sublist_shift)
    load(sublist_shift{subjInd})
    GA_struct = [];
    cfg = [];
    cfg.size_range = [0.05 5.5];
    cfg.binWin = 0.25;
    cfg.binstep = 0.05;
    subNumber = str2double(sublist_shift{subjInd}(end-5:end-4));
    disp(['subNumber:' num2str(subNumber)]);
    if mod(subNumber,2) == 0
        cfg.trigs_left = 712;
        cfg.trigs_right = 711;
    elseif mod(subNumber,2) == 1
        cfg.trigs_left = 711;
        cfg.trigs_right = 712;
    end
    [rate_size,bin_range]=gazeShiftRateOverSize(cfg, data_shift.shift, data_shift.trialinfo);
    
    GA_struct.toward(subjInd,:,:) = rate_size.toward;
    GA_struct.away(subjInd,:,:) = rate_size.away;
    GA_struct.diff(subjInd,:,:) = rate_size.diff;
end

GA_struct.bin_range = bin_range;
GA_struct.time = data_shift.time;

% save the group data
save([write_dir_eye filesep 'GA_shift_rateAndsize_WM'] ,'GA_struct')

%% plotting the figure: LTM
load([write_dir_eye filesep 'GA_shift_rateAndsize_LTM']);
cmap = brewermap([],'*RdBu');
xli = [-0.3 2];
figure('position', [100 100 1200 300])
subplot(1,3,1)
% difference figure

hz2plot= squeeze(nanmean(GA_struct.diff,1));
contourf(GA_struct.time,GA_struct.bin_range,hz2plot,50,'linecolor','none')
maxValue = max(max(max(hz2plot)), abs(min(min(hz2plot)))); 
% maxValue = 0.1;
caxis([-maxValue maxValue])
colorbar
colormap(cmap);
xlabel('Time since cue onset (s)')
ylabel('Shift size (degree)')
title('Toward - Away')
hold on
plot([0 0], [GA_struct.bin_range(1), GA_struct.bin_range(end)], '--k')
plot([1 1], [GA_struct.bin_range(1), GA_struct.bin_range(end)], '--k')
xlim(xli)

% toward
subplot(1,3,2)
hz2plot= squeeze(nanmean(GA_struct.toward,1));
contourf(GA_struct.time,GA_struct.bin_range,hz2plot,50,'linecolor','none')
maxValue = max([max(max(squeeze(nanmean(GA_struct.toward,1)))) max(max(squeeze(nanmean(GA_struct.away,1))))]);
% maxValue = 0.1;
caxis([-maxValue maxValue])
colorbar
colormap(cmap);
xlabel('Time since cue onset (s)')
ylabel('Shift size (degree)')
title('Toward')
hold on
plot([0 0], [GA_struct.bin_range(1), GA_struct.bin_range(end)], '--k')
plot([1 1], [GA_struct.bin_range(1), GA_struct.bin_range(end)], '--k')
xlim(xli)

% away
subplot(1,3,3)
hz2plot= squeeze(nanmean(GA_struct.away,1));
contourf(GA_struct.time,GA_struct.bin_range,hz2plot,50,'linecolor','none')
% maxValue = 0.1;
caxis([-maxValue maxValue])
colorbar
colormap(cmap);
xlabel('Time since cue onset (s)')
ylabel('Shift size (degree)')
title('away')
hold on
plot([0 0], [GA_struct.bin_range(1), GA_struct.bin_range(end)], '--k')
plot([1 1], [GA_struct.bin_range(1), GA_struct.bin_range(end)], '--k')
xlim(xli)



%% plotting the figure: WM
load([write_dir_eye filesep 'GA_shift_rateAndsize_WM']);
cmap = brewermap([],'*RdBu');
xli = [-0.3 2];
figure('position', [100 100 1200 300])
subplot(1,3,1)
% difference figure

hz2plot= squeeze(nanmean(GA_struct.diff,1));
contourf(GA_struct.time,GA_struct.bin_range,hz2plot,50,'linecolor','none')
maxValue = max(max(max(hz2plot)), abs(min(min(hz2plot)))); 
% maxValue = 0.1;
caxis([-maxValue maxValue])
colorbar
colormap(cmap);
xlabel('Time since cue onset (s)')
ylabel('Shift size (degree)')
title('Toward - Away')
hold on
plot([0 0], [GA_struct.bin_range(1), GA_struct.bin_range(end)], '--k')
plot([1 1], [GA_struct.bin_range(1), GA_struct.bin_range(end)], '--k')
xlim(xli)

% toward
subplot(1,3,2)
hz2plot= squeeze(nanmean(GA_struct.toward,1));
contourf(GA_struct.time,GA_struct.bin_range,hz2plot,50,'linecolor','none')
maxValue = max([max(max(squeeze(nanmean(GA_struct.toward,1)))) max(max(squeeze(nanmean(GA_struct.away,1))))]);
% maxValue = 0.1;
caxis([-maxValue maxValue])
colorbar
colormap(cmap);
xlabel('Time since cue onset (s)')
ylabel('Shift size (degree)')
title('Toward')
hold on
plot([0 0], [GA_struct.bin_range(1), GA_struct.bin_range(end)], '--k')
plot([1 1], [GA_struct.bin_range(1), GA_struct.bin_range(end)], '--k')
xlim(xli)

% away
subplot(1,3,3)
hz2plot= squeeze(nanmean(GA_struct.away,1));
contourf(GA_struct.time,GA_struct.bin_range,hz2plot,50,'linecolor','none')
% maxValue = 0.1;
caxis([-maxValue maxValue])
colorbar
colormap(cmap);
xlabel('Time since cue onset (s)')
ylabel('Shift size (degree)')
title('away')
hold on
plot([0 0], [GA_struct.bin_range(1), GA_struct.bin_range(end)], '--k')
plot([1 1], [GA_struct.bin_range(1), GA_struct.bin_range(end)], '--k')
xlim(xli)


%% plot average toward rate and away rate: LTM
toward_mean = [];
away_mean = [];
sublist_shift = get_subFiles([write_dir_eye filesep 'gaze_shift']);
slideWin4rate = 50;
for subjInd = 1:length(sublist_shift)
    load(sublist_shift{subjInd})
    subNumber = str2double(sublist_shift{subjInd}(end-5:end-4));
    disp(['subNumber:' num2str(subNumber)]);
    if mod(subNumber,2) == 0
        trigs_left = 721;
        trigs_right = 722;
    elseif mod(subNumber,2) == 1
        trigs_left = 722;
        trigs_right = 721;
    end

    % select left and right trial 
    sel_left = ismember(data_shift.trialinfo,trigs_left);
    sel_right = ismember(data_shift.trialinfo,trigs_right);
    
    % select left and right shift 
    shift_left = data_shift.shift < 0;
    shift_right = data_shift.shift > 0;
    
    % get toward and away saccade
    toward = (mean(shift_left(sel_left,:)) + mean(shift_right(sel_right,:))) ./ 2;
    away = (mean(shift_left(sel_right,:)) + mean(shift_right(sel_left,:))) ./ 2;
    
    % get the shift rate through movemean method
    toward_rate = smoothdata(toward,2,'movmean',slideWin4rate)*1000;
    away_rate = smoothdata(away,2,'movmean',slideWin4rate)*1000;

    toward_mean = [toward_mean;toward_rate];
    away_mean = [away_mean;away_rate];
end

% Calculate the mean and 95% CI
mean_toward = mean(toward_mean, 1); % Mean across subjects, size (1, 1149)
mean_away = mean(away_mean, 1); % Mean across subjects, size (1, 1149)

% Standard error of the mean
sem_toward = std(toward_mean, 0, 1) / sqrt(size(toward_mean, 1));
sem_away = std(away_mean, 0, 1) / sqrt(size(away_mean, 1));

% 95% Confidence Interval
ci_toward = 1.96 * sem_toward;
ci_away = 1.96 * sem_away;

% Time points
times = 1:size(toward_mean, 2); %(1:2149)
real_time_points = times - 250; %(-250:1900)
toi = 50:1250; %(-200:1000)
times_to_plot = real_time_points(toi);
mean_toward_plot = mean_toward(toi);
ci_toward_plot = ci_toward(toi);
mean_away_plot = mean_away(toi);
ci_away_plot = ci_away(toi);

toward_mean_plot = toward_mean(:,toi);
away_mean_plot = away_mean(:,toi);

% Plot
f = figure;
figureStartup;
hold on;

% Plotting toward
fill([times_to_plot, fliplr(times_to_plot)], ...
    [mean_toward_plot + ci_toward_plot, fliplr(mean_toward_plot - ci_toward_plot)], ...
    'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Shaded area for CI

l1 = plot(times_to_plot, mean_toward_plot, 'b', 'LineWidth', 3); % Mean line

% Plotting away
fill([times_to_plot, fliplr(times_to_plot)], ...
    [mean_away_plot + ci_away_plot, fliplr(mean_away_plot - ci_away_plot)], ...
    'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Shaded area for CI
l2 = plot(times_to_plot, mean_away_plot, 'r', 'LineWidth', 3); % Mean line

cfg = [];
cfg.xax = -200:1000;
cfg.npermutations = 5000;
cfg.clusterStatEvalaluationAlpha= 0.05;
cfg.nsub=44;
cfg.statMethod = 'montecarlo';
cfg.time_i = [-100 600];
state_t = cluster_perm_1D(cfg,toward_mean_plot,away_mean_plot);

mask_xxx = double(state_t.mask); mask_xxx(mask_xxx==0) = nan;
plot(-100:600, mask_xxx * 0.7, 'k', 'LineWidth', 3);

% change x sticks:
xticks(0:400:800);
yticks(0:0.5:2);

% Labels and legend
xlabel('Time relative to cue onset (ms)');
ylabel('Rate (Hz)');
legend([l1,l2],'Toward','Away');
legend box off
hold off;

%% plot average toward rate and away rate: WM
toward_mean = [];
away_mean = [];
sublist_shift = get_subFiles([write_dir_eye filesep 'gaze_shift']);
slideWin4rate = 50;
for subjInd = 1:length(sublist_shift)
    load(sublist_shift{subjInd})
    subNumber = str2double(sublist_shift{subjInd}(end-5:end-4));
    disp(['subNumber:' num2str(subNumber)]);
    if mod(subNumber,2) == 0
        trigs_left = 712;
        trigs_right = 711;
    elseif mod(subNumber,2) == 1
        trigs_left = 711;
        trigs_right = 712;
    end

    % select left and right trial 
    sel_left = ismember(data_shift.trialinfo,trigs_left);
    sel_right = ismember(data_shift.trialinfo,trigs_right);
    
    % select left and right shift 
    shift_left = data_shift.shift < 0;
    shift_right = data_shift.shift >0;
    
    % get toward and away saccade
    toward = (mean(shift_left(sel_left,:)) + mean(shift_right(sel_right,:))) ./ 2;
    away = (mean(shift_left(sel_right,:)) + mean(shift_right(sel_left,:))) ./ 2;
    
    % get the shift rate through movemean method
    toward_rate = smoothdata(toward,2,'movmean',slideWin4rate)*1000;
    away_rate = smoothdata(away,2,'movmean',slideWin4rate)*1000;

    toward_mean = [toward_mean;toward_rate];
    away_mean = [away_mean;away_rate];
end

% Calculate the mean and 95% CI
mean_toward = mean(toward_mean, 1); % Mean across subjects, size (1, 1149)
mean_away = mean(away_mean, 1); % Mean across subjects, size (1, 1149)

% Standard error of the mean
sem_toward = std(toward_mean, 0, 1) / sqrt(size(toward_mean, 1));
sem_away = std(away_mean, 0, 1) / sqrt(size(away_mean, 1));

% 95% Confidence Interval
ci_toward = 1.96 * sem_toward;
ci_away = 1.96 * sem_away;

% Time points
times = 1:size(toward_mean, 2); %(1:2149)
real_time_points = times - 250; %(-250:1900)
toi = 50:1250; %(-200:1000)
times_to_plot = real_time_points(toi);
mean_toward_plot = mean_toward(toi);
ci_toward_plot = ci_toward(toi);
mean_away_plot = mean_away(toi);
ci_away_plot = ci_away(toi);
toward_mean_plot = toward_mean(:,toi);
away_mean_plot = away_mean(:,toi);

% Plot
f = figure;
figureStartup;
hold on;

% Plotting toward
fill([times_to_plot, fliplr(times_to_plot)], ...
    [mean_toward_plot + ci_toward_plot, fliplr(mean_toward_plot - ci_toward_plot)], ...
    'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Shaded area for CI

l1 = plot(times_to_plot, mean_toward_plot, 'b', 'LineWidth', 3); % Mean line

% Plotting away
fill([times_to_plot, fliplr(times_to_plot)], ...
    [mean_away_plot + ci_away_plot, fliplr(mean_away_plot - ci_away_plot)], ...
    'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Shaded area for CI
l2 = plot(times_to_plot, mean_away_plot, 'r', 'LineWidth', 3); % Mean line

cfg = [];
cfg.xax = -200:1000;
cfg.npermutations = 5000;
cfg.clusterStatEvalaluationAlpha= 0.05;
cfg.nsub = 44;
cfg.statMethod = 'montecarlo';
cfg.time_i = [-100 600];
state_t = cluster_perm_1D(cfg,toward_mean_plot,away_mean_plot);

mask_xxx = double(state_t.mask); mask_xxx(mask_xxx==0) = nan;
plot(-100:600, mask_xxx * 0.7, 'k', 'LineWidth', 5);

% change x sticks:
xticks(0:400:800);
yticks(0:0.5:2);

% Labels and legend
xlabel('Time relative to cue onset (ms)');
ylabel('Rate (Hz)');
legend([l1,l2],'Toward','Away');
legend box off
hold off;
% plot figures for three lesion groups
clearvars; clc; close all
%   Rows:
%       first row: eGFP
%       second row: vOFC
%       third row: BLA
%       fourth row: ACC 
%   Columns: all phases
%       first: D (100/0)
%       second: R1(100/0)
%       third:  R2(100/0)
%       fourth: R3(90/10)
%       fifth:  R4(90/10)

% load data files
output_dir = "output/plot_data/";    % specify directory

load(output_dir+"block_output.mat",'block_output'); % whole-block stats
load(output_dir+"trial_output100.mat",'trial_output'); % trials stats (cumulative/running avg)

disp("Output files loaded.");

%% settings

% Choose task:
% taskType = "SO";        % Stimulus-Outcome task
taskType = "AO";      % Action-Outcome task

font_size = 12;

% CNO group colors
red = [.9 .0 .0];  % red
blue = [.1 .3 .8];  % blue
pink = [255 130 121]./255; % pink

group_colors = {red, blue};

% Choose metric to plot:
%     MET.name = "pbetter";       MET.label = "P(Correct)";
%     MET.name = "pwin";          MET.label = "P(Win)";
%     MET.name = "pstay";          MET.label = "P(Stay)";
%     MET.name = "RI_BW";          MET.label = "RI_{BW}";
%     MET.name = "RI_B";          MET.label = "RI_{Better}";
%     MET.name = "RI_W";          MET.label = "RI_{Worse}";
%     MET.name = "winstay";       MET.label = "WinStay";
    MET.name = "loseswitch";    MET.label = "LoseSwitch";
%     MET.name = "H_str";         MET.label = "H(Str)";
%     MET.name = "ERDS";          MET.label = "ERDS";
%         MET.name = "MIRS";          MET.label = "MIRS";

% Exists only for SO task:
%     MET.name = "pstay_Loc";          MET.label = "P(Stay_{Action})";
%     MET.name = "ERDS_Loc";      MET.label = "ERDS_{Action}";
%         MET.name = "MIRS_Loc";      MET.label = "MIRS_{Action}";
%     MET.name = "RI_LR";          MET.label = "RI_{LR}";
%     MET.name = "RI_L";          MET.label = "RI_{Left}";
%     MET.name = "RI_R";          MET.label = "RI_{Right}";

% 0. Cut-offs based on Discrimination performance
plot_summary_fig = 1;
[block_output, trial_output] = apply_cutoff_D_performance(taskType, block_output, trial_output, plot_summary_fig);

%% Fig1. First five days in each phase
close all

[outVals, subLabels] = figure_first_five_days(block_output, taskType, MET, group_colors, font_size);

%% Fig2. Transition: Days around reversal
close all

[outVals, subLabels] = figure_days_around_rev(block_output, taskType, MET, group_colors, font_size);

%% Fig3. First N=100 trials in each phase
close all

smooth_window = 10;     % additional smoothing window size

% Choose stats version
% plot_type1 = "Cum";    % cumulative average
% plot_type2 = "Run";    % running window average

[outValsCum, subLabels] = figure_first_N_trials(trial_output, "Cum", taskType, MET, group_colors, font_size, 1);   % no smoothing for cumulative stats
[outValsRun, ~] = figure_first_N_trials(trial_output, "Run", taskType, MET, group_colors, font_size, smooth_window);   % smooth using 10-point moving avg.

%% Fig4. Transition: Trials around reversal
close all

smooth_window = 10;     % additional smoothing window size    

% Choose stats version
% plot_type1 = "Cum";    % cumulative average
% plot_type2 = "Run";    % running window average

[outValCum, subLabels] = figure_trials_around_rev(trial_output, "Cum", taskType, MET, group_colors, font_size, smooth_window);     % no smoothing for cumulative stats
% [outValRun, ~] = figure_trials_around_rev(trial_output, "Run", taskType, MET, group_colors, font_size, smooth_window);      % smooth using 10-point moving avg.



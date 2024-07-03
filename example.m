% Add utils path
clc; clear; close all;  % Clear command window, workspace, and close all figures
addpath("./utils/");      % Add the utils directory to the search path
addpath("./group_analysis/");      % Add the directory to the search path
%% Data Initialization
% Define the number of subjects, timepoints, components, and sampling time
num_subjects = 5;       % Number of subjects
num_timepoints = 150;   % Number of timepoints
num_components = 10;    % Number of components
Tr = 2;                 % repetition time

% Generate example post-processed fMRI timecourses
subTcs = randn(num_subjects, num_timepoints, num_components);

% Define groupings for subjects (example, 1 is schizophrenia, 2 is controls)
groupings = randi([1 2], num_subjects, 1);

% Obtain warp elasticity
% win_size = 45;
% order = 16;
% Fpass = 0.12;
% Fstop = 0.15;
% [~, ~, we] = calculate_warp_elasticity(subTcs, win_size, order, Fpass, Fstop, Tr);

% Generate example warp elasticity 4D matrix: subjects x time x component x
% component
we = randi([-100 100], num_subjects, num_timepoints, num_components, num_components) / 100;

%% Preparing for k-means clustering
% Concatenate subjects and time of the trFNC to get a 3D matrix
we_reshape = reshape(we, num_subjects*num_timepoints, num_components, num_components);

% Since the component pairs (num_components x num_components) are
% symmetric, vectorize it by selecting only the lower triangle (or upper)
we_k = icatb_mat2vec(we_reshape);

%% Perform K-means
% Define the number of clusters
cluster_num = 4;
replications = 20;

% %Perform k-means
% [cluster_idx, cluster_c, cluster_sumd, cluster_D] = calculate_kmeans(we_k, cluster_num, 'cityblock', replications);
% cluster_idx = reshape(cluster_idx, subjects_cnt, []); 

% Generate random cluster indices for warp elasticity to mimic kmeans index
% output
cluster_idx = randi([1 cluster_num], num_subjects * num_timepoints, 1);
% Reshape into subjects x timepoints
cluster_idx = reshape(cluster_idx, num_subjects, num_timepoints);
%% Extract groups based on cluster indices
group1 = []; % Initialize group1
group2 = []; % Initialize group2

% Separate subjects into two groups based on cluster index
group1 = find(groupings == 1);
group2 = find(groupings == 2);

% Count the number of subjects in each group
num_group1 = length(group1);
num_group2 = length(group2);

% Extract trFNC data for each group
trFNC_group1 = we(group1, :, :, :);
trFNC_group2 = we(group2, :, :, :);

% Reshape cluster indices for each group
cluster_idx_group1 = cluster_idx(group1, :);
cluster_idx_group2 = cluster_idx(group2, :);

%% Extract subject states for each group
% Calculate group centroids for each group
sub_centroids_group1 = group_sub_centroids(trFNC_group1, cluster_idx_group1, cluster_num);
sub_centroids_group2 = group_sub_centroids(trFNC_group2, cluster_idx_group2, cluster_num);

% Reshape and vectorize the centroids for statistical analysis
sub_centroids_group1 = reshape(icatb_mat2vec(reshape(sub_centroids_group1, num_group1*cluster_num, num_components, num_components)), num_group1, cluster_num, []);
sub_centroids_group2 = reshape(icatb_mat2vec(reshape(sub_centroids_group2, num_group2*cluster_num, num_components, num_components)), num_group2, cluster_num, []);

%% Perform statistical analysis on group states
[p_values, fdr_p_values, s_values, sign_values] = compute_cluster_stats(sub_centroids_group1, sub_centroids_group2, cluster_num);

%% Analyze groups' cluster dynamics
% Compute cluster analysis metrics for each subject in each group
[mdt_group1, fr_group1, tm_group1] = trFNC_cluster_analysis(cluster_idx_group1, cluster_num);
[mdt_group2, fr_group2, tm_group2] = trFNC_cluster_analysis(cluster_idx_group2, cluster_num);

%% Perform t-tests on mean dwell time (mdt), fractional occupancy (fr), and transition matrix (tm)
% Note: Assuming the implementation of t-tests will follow here

% Example of performing t-test (uncomment and adapt as needed)
% [h_mdt, p_mdt] = ttest2(mdt_group1, mdt_group2);
% [h_fr, p_fr] = ttest2(fr_group1, fr_group2);
% [h_tm, p_tm] = ttest2(mdt_group1, mdt_group2);



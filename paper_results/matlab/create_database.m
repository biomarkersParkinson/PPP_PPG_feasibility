%% This script creates the total database in matlab for the feasibility study 
%%--- Database HR measurements---%%
% Containing:
% 1. POMU number (string)
% 2. Age (numeric)
% 3. Sex (boolean)
% 4. Years since diagnosis (numeric)
% 5a. Specify ECG removal (string) --> nog niet in database


%% Load all relevant clinical data in memory
clearvars; close all; clc;
general_data = readtable('path_to_general_data.csv');       % load this file since it is required
concatenateFunction = @(x) ['POMU' x(1:min(16, numel(x)))];
general_data.id = cellfun(concatenateFunction, general_data.id, 'UniformOutput', false);

scopa_data = readtable('path_to_scopa_data.csv');
updrs1b_data = readtable('path_to_updrs1b_data.csv');
updrs3_data = readtable('path_to_updrs3_data.csv');
updrs4_data = readtable('path_to_updrs4_data.csv');
pase_data = readtable('path_to_pase_data.csv');
ecg_data = readtable('path_to_ecg_data.csv');
%% Clinical database
rel_general_columns = general_data(:, {'id', 'Age', 'Gender', 'MonthSinceDiag', 'Up3OfRAmpArmYesDev', 'Up3OfRAmpArmNonDev', 'Up3OfRAmpJaw', 'Up3OfRAmpLegYesDev','Up3OfRAmpLegNonDev', 'Up3OfHoeYah', 'Up3OnHoeYah'});
rel_scopa_columns = scopa_data(:, {'id', 'ScopaAut_tot'});
rel_updrs1b_columns = updrs1b_data(:, {'id', 'UPDRS_1b_total'});
database_clinical = outerjoin(rel_general_columns, rel_scopa_columns, "Keys", 'id', 'MergeKeys',true);
database_clinical = outerjoin(database_clinical, rel_updrs1b_columns, "Keys", 'id', 'MergeKeys',true);
database_clinical = outerjoin(database_clinical, updrs3_data, "Keys", 'id', 'MergeKeys',true);
database_clinical = outerjoin(database_clinical, updrs4_data, "Keys", 'id', 'MergeKeys',true);
database_clinical = outerjoin(database_clinical, pase_data, "Keys", 'id', 'MergeKeys', true);

%% SQA database calculation
week_array = [0, 1];
sqa_data_root = "path_to_sqa_output\";

threshold = 0.5;  % Threshold for posterior probability
unix_ticks_ms = 1000;
timezone = 'Europe/Amsterdam';  % Everything is stored in uinx until this!
day_hrs = [6, 23];
night_hrs = [0, 5];

fs_ppg = 1;  % sampling rate for sample_prob_final --> now return is 1 per s --> sufficient

database_sqa_week0 = table;
database_sqa_week1 = table;


for week = week_array
    sqa_data_path = fullfile(sqa_data_root, ['WatchData.PPG.Week' num2str(week)]);
    sqa_data_list = dir(fullfile(sqa_data_path)); sqa_data_list(ismember( {sqa_data_list.name}, {'.', '..'})) = [];   
    
    id_list = {};
    mean_sync_time_list = [];
    mean_sync_time_day_list= [];
    mean_sync_time_night_list = [];
    mean_hq_time_list = [];
    mean_hq_time_day_list = [];
    mean_hq_time_night_list = [];
        
    for K = 1:length(sqa_data_list)
    
        hour_counts = [];
        hq_hour_counts = [];
        hq_count_day = [];
        hq_count_night = [];
    
        sqa_output_list = dir(fullfile(sqa_data_path, sqa_data_list(K).name,'*meta.json'));   % create the segment list for every different SQA segment
            
        meta_path_sqa = fullfile(sqa_output_list.folder, sqa_output_list.name); % Get the name of the SQA output file
        [metadata_list_sqa, data_list_sqa] = load_tsdf_metadata_from_path(meta_path_sqa); % Load the metadata of the SQA output file
    
        ppg_prob_idx = tsdf_values_idx(metadata_list_sqa, 'ppg'); % Get the index of the ppg field in the SQA output file
        ppg_post_prob_total = data_list_sqa{ppg_prob_idx};    % Get the PPG probability data of the SQA output file
    
        imu_idx = tsdf_values_idx(metadata_list_sqa, 'sqa_imu'); % Get the index of the imu field in the SQA output file
        imu_label = data_list_sqa{imu_idx};    % Get the IMU label data of the SQA output file
    
        time_idx = tsdf_values_idx(metadata_list_sqa, 'time'); % Get the index of the time field in the SQA output file
        time_class_unix = data_list_sqa{time_idx};    % Get the time data of the SQA output file
        
        sync_idx = tsdf_values_idx(metadata_list_sqa, 'sync'); % Get the index of the sync field in the SQA output file
        data_sync = data_list_sqa{sync_idx};    % Get the sync data of the SQA output file
        data_sync_zeros = all(data_sync == 0, 2); % Find rows containing only zeros
        data_sync(data_sync_zeros, :) = [];  % Remove rows containing only zeros
        n_segments_sync = size(data_sync,1);  % Get the number of segments in the SQA output file
    
        for i = 1:length(data_sync(:,4))
            if i == 1
                start_end_indices(i, 1) = 1;
            else
                start_end_indices(i, 1) = start_end_indices(i-1, 2) + 1;
            end
            start_end_indices(i, 2) = sum(data_sync(1:i,4));
        end
        
        for n = 1:n_segments_sync
            class_start = start_end_indices(n,1);
            class_end = start_end_indices(n,2); 
            class_ppg_segment = ppg_post_prob_total(class_start:class_end);
            class_acc_segment = imu_label(class_start:class_end);
            class_time_segment = time_class_unix(class_start:class_end);

            ppg_post_prob = sample_prob_final(class_ppg_segment, fs_ppg, class_acc_segment);
            %ppg_post_prob = sample_prob_final(class_ppg_segment, fs_ppg);
            ppg_post_prob = ppg_post_prob(1:end-5);  % update with overlap so we don't have to overwrite sample_prob_final

            sqa_datetimes = datetime(class_time_segment, 'ConvertFrom', 'posixtime', 'TimeZone', timezone);      % convert to datetime

            sqa_hour = hour(sqa_datetimes);               % extract hour from datetimes!
            hour_counts(n,:) = sum_hours(sqa_hour);
            hq_hour_counts(n,:) = sum_hours(sqa_hour(ppg_post_prob>=threshold));

        end
    
        total_hours = sum(hour_counts,1);
        total_hq_hours = sum(hq_hour_counts,1);
        
        mean_sync_time = sum(total_hours)/(7*60*60);            % mean wearing time per day in hours
        mean_sync_time_day = sum(total_hours(1,day_hrs(1)+1:day_hrs(2)+1))/(7*60*60);
        mean_sync_time_night = sum(total_hours(1,night_hrs(1)+1:night_hrs(2)+ 1))/(7*60*60);
        mean_hq_time = sum(total_hq_hours)/(7*60*60);
        mean_hq_time_day = sum(total_hq_hours(1,day_hrs(1)+1:day_hrs(2)+ 1))/(7*60*60);
        mean_hq_time_night = sum(total_hq_hours(1,night_hrs(1)+1:night_hrs(2)+ 1))/(7*60*60);

        id_list{K} = sqa_data_list(K).name;
        mean_sync_time_list(K) = mean_sync_time;
        mean_sync_time_day_list(K) = mean_sync_time_day;
        mean_sync_time_night_list(K) = mean_sync_time_night;
        mean_hq_time_list(K) = mean_hq_time;
        mean_hq_time_day_list(K) = mean_hq_time_day;
        mean_hq_time_night_list(K) = mean_hq_time_night;

        if mod(K, 10) == 0
            disp(['Processed week ' num2str(week) ', Finished sqa iterations: ' num2str(K)]);
        end
        
    end
    
   if week == 0
        database_sqa_week0 = table(id_list', mean_sync_time_list', mean_sync_time_day_list', mean_sync_time_night_list', mean_hq_time_list', mean_hq_time_day_list', mean_hq_time_night_list',...
            'VariableNames', {'id', ['mean_sync_time_week' num2str(week)], ['mean_sync_time_day_week' num2str(week)], ['mean_sync_time_night_week' num2str(week)], ['mean_hq_time_week' num2str(week)], ['mean_hq_time_day_week' num2str(week)], ['mean_hq_time_night_week' num2str(week)]});
   else
        database_sqa_week1 = table(id_list', mean_sync_time_list', mean_sync_time_day_list', mean_sync_time_night_list', mean_hq_time_list', mean_hq_time_day_list', mean_hq_time_night_list',...
            'VariableNames', {'id', ['mean_sync_time_week' num2str(week)], ['mean_sync_time_day_week' num2str(week)], ['mean_sync_time_night_week' num2str(week)], ['mean_hq_time_week' num2str(week)], ['mean_hq_time_day_week' num2str(week)], ['mean_hq_time_night_week' num2str(week)]});
    end
    
end


%% HR parameter calculation
week_array = [0, 1];
hr_data_root = "path_to_hr_output\";
timezone = 'Europe/Amsterdam';
day_hrs = [8, 22];
night_hrs = [0, 5];

database_hr_week0 = table;
database_hr_week1 = table;

for week = week_array
    hr_data_path = fullfile(hr_data_root, ['WatchData.PPG.Week' num2str(week)]);
    hr_data_list = dir(fullfile(hr_data_path)); hr_data_list(ismember( {hr_data_list.name}, {'.', '..'})) = [];

    id_list = {};
    HR_rest_day_list = [];
    HR_rest_night_list = [];
    HR_max_day_list = [];
    HR_max_night_list = [];

    for K = 1:length(hr_data_list)
        hr_output_list = dir(fullfile(hr_data_path, hr_data_list(K).name, '*meta.json'));   % create the segment list for every different HR segment      
        meta_path_hr = fullfile(hr_output_list.folder, hr_output_list.name); % Get the name of the SQA output file
        [metadata_list_hr, data_list_hr] = load_tsdf_metadata_from_path(meta_path_hr); % Load the metadata of the SQA output file
        
        hr_values_idx = tsdf_values_idx(metadata_list_hr, 'values'); % Get the index of the ppg field in the SQA output file
        hr_values = data_list_hr{hr_values_idx};    % Get the PPG probability data of the SQA output file
        
        hr_time_idx = tsdf_values_idx(metadata_list_hr, 'time'); % Get the index of the imu field in the SQA output file
        hr_time_unix = data_list_hr{hr_time_idx};    % Get the IMU label data of the SQA output file
        
        hr_param = calc_hr_param(hr_time_unix, hr_values, day_hrs, night_hrs, timezone);
       
        id_list{K} = hr_data_list(K).name;
        HR_rest_day_list(K) = hr_param(1);
        HR_rest_night_list(K) = hr_param(2);
        HR_max_day_list(K) = hr_param(3);
        HR_max_night_list(K) = hr_param(4);
        HR_rest_day_smooth_list(K) = hr_param(5);
        HR_rest_night_smooth_list(K) = hr_param(6);

        % Display information every 10 iterations
        if mod(K, 10) == 0
            disp(['Processed week ' num2str(week) ', Finished hr iterations: ' num2str(K)]);
        end

    end
    
    if week == 0
        database_hr_week0 = table(id_list', HR_rest_day_list', HR_rest_night_list', HR_max_day_list', HR_max_night_list', HR_rest_day_smooth_list', HR_rest_night_smooth_list', ...
            'VariableNames', {'id', ['HR_rest_day_week' num2str(week)], ['HR_rest_night_week' num2str(week)], ['HR_max_day_week' num2str(week)], ['HR_max_night_week' num2str(week)], ['HR_rest_day_smooth_week' num2str(week)], ['HR_rest_night_smooth_week' num2str(week)]});
    else
        database_hr_week1 = table(id_list', HR_rest_day_list', HR_rest_night_list', HR_max_day_list', HR_max_night_list', HR_rest_day_smooth_list', HR_rest_night_smooth_list', ...
            'VariableNames', {'id', ['HR_rest_day_week' num2str(week)], ['HR_rest_night_week' num2str(week)], ['HR_max_day_week' num2str(week)], ['HR_max_night_week' num2str(week)], ['HR_rest_day_smooth_week' num2str(week)], ['HR_rest_night_smooth_week' num2str(week)]});
    end

end

%writetable(database_hr_week0, "C:\Users\z863160\Documents\AI4P\PPG\Artikel feasibility\data\final database\20240515_week0_HR_unadjusted.csv")
%% Combine clinical, sqa and hr databases
final_database = database_clinical;
final_database = outerjoin(final_database, database_sqa_week0, "Keys", 'id', 'MergeKeys',true);
final_database = outerjoin(final_database, database_hr_week0, "Keys", 'id', 'MergeKeys',true);
final_database = outerjoin(final_database, database_sqa_week1, "Keys", 'id', 'MergeKeys',true);
final_database = outerjoin(final_database, database_hr_week1, "Keys", 'id', 'MergeKeys',true);


%% Add inclusion criteria
final_database_incl = database_clinical;
min_hour_mean_24h = 12;
min_hour_mean_night = 3;
min_hq_hour_mean_24h = 4;
min_hq_hour_mean_day = 2;

inclusion_table = [];

for week = week_array
    incl_wear_time = {};
    incl_hq_time = {};
    id_list = {};

    for K = 1:height(final_database)
        str_time = ['mean_sync_time_week' num2str(week)];
        str_time_night = ['mean_sync_time_night_week' num2str(week)];
        mean_sync_time = table2array(final_database(K, str_time));
        mean_sync_time_night = table2array(final_database(K, str_time_night));
       
        if mean_sync_time < min_hour_mean_24h || mean_sync_time_night < min_hour_mean_night || isnan(mean_sync_time) || isnan(mean_sync_time_night)
            incl_wear_time{K} = 0;
        else
            incl_wear_time{K} = 1;
        end

        str_hq_time = ['mean_hq_time_week' num2str(week)];
        str_hq_time_night = ['mean_hq_time_night_week' num2str(week)];
        mean_hq_time = table2array(final_database(K, str_time));
        mean_hq_time_night = table2array(final_database(K, str_time_night));

        if mean_hq_time < min_hq_hour_mean_24h || mean_hq_time_night < min_hq_hour_mean_day || isnan(mean_hq_time) || isnan(mean_hq_time_night)
            incl_hq_time{K} = 0;
        else
            incl_hq_time{K} = 1;
        end

        id_list{K} = final_database_incl.id{K};

    
    end
    if week == 0
        inclusion_table_week0 = table(id_list', incl_wear_time', incl_hq_time', ...
            'VariableNames', {'id', ['Incl_wear_time_week' num2str(week)], ['Incl_hq_time' num2str(week)]});
    else
        inclusion_table_week1 = table(id_list', incl_wear_time', incl_hq_time', ...
            'VariableNames', {'id', ['Incl_wear_time_week' num2str(week)], ['Incl_hq_time' num2str(week)]});
    end

end


%%----- ECG inclusion ----- %%

inclusion_table_cardiac = ecg_data(:, {'id', 'beta_user', 'exclusion_rhythm'});
%% Combine clinical, sqa and hr databases
final_database_incl = outerjoin(final_database_incl, inclusion_table_cardiac, "Keys", 'id', 'MergeKeys', true);
final_database_incl = outerjoin(final_database_incl, inclusion_table_week0, "Keys", 'id', 'MergeKeys', true);
final_database_incl = outerjoin(final_database_incl, inclusion_table_week1, "Keys", 'id', 'MergeKeys', true);
final_database_incl = outerjoin(final_database_incl, database_sqa_week0, "Keys", 'id', 'MergeKeys', true);
final_database_incl = outerjoin(final_database_incl, database_hr_week0, "Keys", 'id', 'MergeKeys', true);
final_database_incl = outerjoin(final_database_incl, database_sqa_week1, "Keys", 'id', 'MergeKeys', true);
final_database_incl = outerjoin(final_database_incl, database_hr_week1, "Keys", 'id', 'MergeKeys', true);

%% Label annotation subjects
directory = "path_to_annotation_data";
all_items = dir(directory);

% Filter out only folders (excluding '.' and '..')
subfolders = {all_items([all_items.isdir] & ~ismember({all_items.name}, {'.', '..'})).name};

% Initialize a cell array to store the refactored names
refactored_names_annotation = cell(size(subfolders));

% Define your custom refactoring function
refactor_function = @(x) ['POMU' x(1:min(16, numel(x)))];

% Loop through each subfolder name
for i = 1:length(subfolders)
    old_name = subfolders{i};
    
    % Apply the custom refactoring function to the folder name
    new_name = refactor_function(old_name);
    
    % Store the refactored name in the cell array
    refactored_names_annotation{i} = new_name;
    
    % Display the old and new folder names (optional)
    fprintf('Refactoring folder: %s -> %s\n', old_name, new_name);
    
end
final_database_incl.training_set(:) = 0;
[~, row_index] = ismember(refactored_names_annotation, final_database_incl.id);
row_index(row_index==0) = [];
final_database_incl.training_set(row_index) = 1;

%% Store csv file

final_database_root = "path_to_final_database";
filename_database = 'feasibility_database.csv';
location_final = fullfile(final_database_root, filename_database);
writetable(final_database_incl, location_final);

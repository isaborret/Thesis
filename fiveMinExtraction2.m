%% extract_all_subjects.m
% Extracts 5 minutes of EEG and ECG immediately BEFORE the awakening
% for all subjects in DREAM dataset 12, except subject 1.
%
% Uses awakening clock time from Records.csv + EDF header start time
% to locate the correct pre-awakening window.
%
% Output per subject: one .mat file, e.g. sub002_noDream_REM_last5min.mat
%
% Each .mat contains:
%   EEG_data     [19 x 76800]   EEG signals
%   ECG_data     [ 1 x 76800]   ECG signal
%   EEG_labels   {1x19}         channel names
%   fs           256             sampling rate (Hz)
%   t            [1 x 76800]    time vector in seconds
%   subject_id   string
%   dream        logical         true = dream recall
%   sleep_stage  string          'REM' or 'N2'
%   n_samples    scalar          (= 76800)
%
% -------------------------------------------------------------------------
% SET THESE PATHS BEFORE RUNNING
% -------------------------------------------------------------------------

data_dir    = 'C:\Users\Isa\OneDrive\Documenten\Unif 2e master\Masterproef_II\data\raw data\';
output_dir = 'C:\Users\Isa\OneDrive\Documenten\Unif 2e master\Masterproef_II\data\derivatives\5_minuteExtractionRawEEGandECG\';
records_csv = 'C:\Users\Isa\OneDrive\Documenten\Unif 2e master\Masterproef_II\data\raw data\Dataset 12 info\Description & docs dataset\Records.csv'; 

% -------------------------------------------------------------------------

%% Setup
fs        = 256;
n_extract = fs * 5 * 60;   % 76800 samples

eeg_ch_idx = 1:19;
ecg_ch_idx = 24;

EEG_labels = {'Fp1','Fp2','F3','F4','Fz','F7','F8', ...
               'C3','C4','Cz','P3','P4','Pz','T3','T4','T5','T6','O1','O2'};

if ~exist(output_dir, 'dir'), mkdir(output_dir); end

%% Load Records.csv and find awakening time column
records   = readtable(records_csv, 'VariableNamingRule', 'preserve');
col_names = records.Properties.VariableNames;
awake_col = find(contains(col_names, 'awakening', 'IgnoreCase', true) | ...
                 contains(col_names, 'Time', 'IgnoreCase', true) & ...
                 ~contains(col_names, 'Age', 'IgnoreCase', true), 1);
fprintf('Awakening time column: "%s"\n\n', col_names{awake_col});

%% Find all EDF files and filter out hidden files
edf_files = dir(fullfile(data_dir, '*.edf'));
files = files(~[files.isdir]);
files = files(~startsWith({files.name}, '.'));
fprintf('Found %d EDF files.\n\n', n_files);

%% Preallocate log
log_entries = cell(n_files, 4);  % filename | dream | stage | status

%% Ignore subjects
ignoresubs = [1]

%% Main loop
for i = 1:n_files

    fname          = edf_files(i).name;
    fpath          = fullfile(edf_files(i).folder, fname);
    [~, file_stem] = fileparts(fname);

    fprintf('[%02d/%02d] %s\n', i, n_files, file_stem);

    %% --- Get subject number and skip subject 1 ---
    subnum_str = regexp(file_stem, '\d+', 'match'); % searches the string file_stem for any sequence of digits (\d+ means "one or more digits"). For 'sub002_noDream_REM' it returns {'002'}. Then str2double(subnum_str{1}) converts that string '002' to the number 2. So it's just extracting the subject number from the filename.
    subnum     = str2double(subnum_str{1});

    if subnum == [ignoresubs]
        fprintf('  Skipping subject 1 (insufficient pre-awakening data).\n');
        log_entries(i,:) = {file_stem, '', '', 'SKIPPED'};
        continue
    end

    %% --- Parse dream status and sleep stage from filename ---
    parts = strsplit(file_stem, '_');
    subject_id = parts{1};

    dream_token = '';
    for p = 1:numel(parts)
        if strcmpi(parts{p}, 'noDream')
            dream_token = 'noDream';
        elseif strcmpi(parts{p}, 'dream') && isempty(dream_token)
            dream_token = 'dream';
        end
    end
    dream = strcmpi(dream_token, 'dream');

    sleep_stage = '';
    for p = 1:numel(parts)
        if strcmpi(parts{p}, 'REM'),    sleep_stage = 'REM';
        elseif strcmpi(parts{p}, 'N2'), sleep_stage = 'N2'; end
    end

    fprintf('  Subject: %s | Dream: %d | Stage: %s\n', subject_id, dream, sleep_stage);

    %% --- Compute awakening sample index ---
    % Step 1: get awakening clock time from Records.csv and convert to seconds
    awake_val   = records{subnum, awake_col}; % fetches the cell at row 2 (for subject 2, or row 7 for subject 7), awakening-time column for subject 2, you mean? — this gives something like {['5:23:00']}
    if iscell(awake_val), awake_val = awake_val{1}; end
    awake_parts = strsplit(strtrim(char(awake_val)), ':'); % splits on the colon, giving {'5', '23', '00'}
    awake_sec   = str2double(awake_parts{1})*3600 + ...
                  str2double(awake_parts{2})*60   + ...
                  str2double(awake_parts{3}); % converts to total seconds: 5×3600 + 23×60 + 0 = 19380 seconds

    % Step 2: get recording start time from EDF header and convert to seconds
    hdr       = edfinfo(fpath);
    start_str = strrep(char(hdr.StartTime), '.', ':');
    sp        = strsplit(start_str, ':');
    start_sec = str2double(sp{1})*3600 + ...
                str2double(sp{2})*60   + ...
                str2double(sp{3});

    % Step 3: elapsed time = difference; correct for midnight crossover
    elapsed_sec      = awake_sec - start_sec;
    if elapsed_sec < 0, elapsed_sec = elapsed_sec + 86400; end
    awakening_sample = round(elapsed_sec * fs);

    fprintf('  Awakening at %.1f min into recording\n', elapsed_sec/60);

    %% --- Load EDF signal data ---
    try
        tbl       = edfread(fpath); %Loads the EDF file into a timetable. 
                                    % Each row is one data record (1 second),each column is one channel, 
                                    % each cell contains a vector of 256 samples.
        tbl_arr   = table2array(tbl); % Converts the timetable to a plain cell array of size [n_records × n_channels]
        n_ch_file = size(tbl_arr, 2); % Counts how many channels are in the file (should be 24).

        first_ch  = cell2mat(tbl_arr(:, 1)); %Takes the entire first column (all records for channel 1) and concatenates all the 1-second vectors into one long vector.
        n_total   = numel(first_ch); %counts its total number of samples, giving us the full recording length.
        data_mat  = zeros(n_ch_file, n_total); %Pre-allocates the full [24 × n_total] matrix with zeros,
        data_mat(1,:) = first_ch(:)'; %fills in the already-computed channel 1
        for ch = 2:n_ch_file %Repeats the same concatenation for every remaining channel. After this loop, data_mat is the complete [24 × n_total] numeric matrix of the entire recording
            data_mat(ch,:) = cell2mat(tbl_arr(:, ch))';
        end
        fprintf('  Loaded: %d ch x %d samples (%.1f min)\n', ...
                n_ch_file, n_total, n_total/fs/60);
    catch ME
        fprintf('  ERROR loading: %s\n  Skipping.\n', ME.message);
        log_entries(i,:) = {file_stem, dream_token, sleep_stage, 'LOAD ERROR'};
        continue
    end

    %% --- Extract 5 minutes ending at awakening sample ---
    end_s    = min(awakening_sample, size(data_mat, 2));
    start_s  = end_s - n_extract + 1;

    EEG_data = data_mat(eeg_ch_idx, start_s:end_s);   % [19 x 76800]
    ECG_data = data_mat(ecg_ch_idx, start_s:end_s);   % [ 1 x 76800]
    t        = (0:n_extract-1) / fs; % t is a time vector where t(1) = 0 seconds (start of the 5-minute window) and t(end) ≈ 300 seconds (end of the window, just before awakening). 
    n_samples = n_extract;

    fprintf('  Window: %.1f to %.1f min from recording start\n', ...
            start_s/fs/60, end_s/fs/60);

    %% --- Save ---
    out_name = fullfile(output_dir, [file_stem '_last5min.mat']);
    save(out_name, 'EEG_data', 'ECG_data', 'EEG_labels', ...
         'fs', 't', 'n_samples', 'subject_id', 'dream', 'sleep_stage', '-v7.3');

    fprintf('  Saved: %s\n', [file_stem '_last5min.mat']);
    log_entries(i,:) = {file_stem, dream_token, sleep_stage, 'OK'};

end

%% Summary
fprintf('\n========== BATCH COMPLETE ==========\n');
ok_count = sum(strcmp(log_entries(:,4), 'OK'));
fprintf('Extracted: %d / %d\n', ok_count, n_files);

log_table = cell2table(log_entries(1:n_files,:), ...
    'VariableNames', {'Filename','Dream','SleepStage','Status'});
writetable(log_table, fullfile(output_dir, 'extraction_log.csv'));
fprintf('Log saved to extraction_log.csv\n');


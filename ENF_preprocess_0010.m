%%%% 
% Jonathan Weinrib
% Date of Origination: 10/17/18
% Date of continuation, and date of origination with this name: 12/25/18
% File number for ease of remembering, as my titles aren't going to stay
% constant: File #0010


%% Summary and explanation of file.
% This file takes in a .mat file. This .mat file was created by reading in
% a wav file (done by the scripts in the Reading_Dataset folder) and
% converting to a .mat.

% Each .mat represents a time domain signal.


% In this file, we first load the time domain signal.
% Next, we bandpass filter the signals.
% Following this we normalize the signals (still in time-domain).

% The next step is to take the STFT of the signal.
% From this we get the estimated ENF (by taking the frequency bin with the
% most power at a given time window to be what we consider the ENF).
% We can also extract other frequency-domain features 




%% Initial Preparation

% Clear Workspace
clc, clear all, close all

% Variable for plotting
plotOn = false;
tic

%% 1) Load time-domain ENF recordings:
num_recordings_per_grid = [9,10,11,11,11,8,11,11,11];
mean_freq_array         = [60,50,60,50,60,50,50,50,60];
grid_letter_list = 'ABCDEFGHI';
%% For testing only: Remove after.
% For testing (so as not to use way too much memory and not to take so
% long), I will use 4 recordings from each of the first 5 grids.
num_recordings_per_grid = [2,2,2,2,2];
mean_freq_array = mean_freq_array(1:length(num_recordings_per_grid));
grid_letter_list = 'ABCDE';
%grid_letter_list = 'BCDFG'; % for this round of testing, since grids A,E, and I are 60 hz, stay with 50hz wones

%%
num_grids = length(num_recordings_per_grid);

%ENF_Files_Array is an 2 dimensional cell array.
%The first dimension is the grid letter
%The second dimension is the recording number.
%So, ENF_Files_Array{grid}{rec_num} will contain the pure ENF recording.
ENF_Files_Array = cell(length(grid_letter_list),max(num_recordings_per_grid));
%ENF_Filter_Freq_Array will for continuity sake stay the same dimensions and
%sizings and stuff as ENF_Files_Array (unless I decide to change that
%later), and will contain the median frequency for that grid.
ENF_Filter_Freq_Array = cell(length(grid_letter_list),max(num_recordings_per_grid));

% Define variables for bandpass filtering:
freq_deviation = 5; %5 Hz allowance to start
fs = 1000; % Sampling frequency
fn = fs/2; % Nyquist sampling rate

%design the butterworth filter
bttr_filt_ord = 3;

for grid = 1:num_grids % Loop through the grids
    
    grid_letter = grid_letter_list(grid); % Set the grid letter variable
    
    for rec_num = 1:num_recordings_per_grid(grid) % loop through each recording in the grid
        recording_name = [grid_letter '_P' num2str(rec_num) '.mat']
        curr_sig = load(recording_name);
        
        % Bandpass filter the current signal
        if grid_letter == 'A' | grid_letter == 'E'|grid_letter == 'I' % JON!! check that it's really grid E that has 60 hz mean)
            mean_freq = 60;
        else 
            mean_freq = 50;
        end
        
        % Regardless of what we chose as the mean_freq, it's now time to
        % save that in ENF_Filter_Freq_Array
        ENF_Filter_Freq_Array{grid,rec_num} = mean_freq;
        %ENF_Files_Array
        
        % Define more variables for BP filtering 
        % cutoff frequencies in Hz
        f1 = mean_freq - freq_deviation; f2 = mean_freq + freq_deviation;
        w1 = f1/fn; w2 = f2/fn; %cutoff frequencies (nyquist)
        
        % Design the butterworth filter 
        [b,a] = butter(bttr_filt_ord,[w1 w2]);
        
        % Now filter the signal
        filtered_sig = filter(b,a,curr_sig.x);
        
        % Now normalize the signal
        scaled_sig = filtered_sig./max(abs(filtered_sig));
        
        % Save the normalized filtered signal into a variable
        ENF_Files_Array{grid,rec_num} = scaled_sig;
        
    end % end for rec_num = 1:num_recordings_per_grid(grid)
    
    
end %end for grid = 1:num_grids




% Now to make looping through the signals Easier
num_sigs = sum(num_recordings_per_grid); % Know how many signals we have. 
snb = 'sig_'; % Sig Name Base
sne = '_final'; %Sig Name End

% initalize array of signal names
sig_name_list = {};
% Loop through all desired signals to create a 
% string for that signal name
for sigNum = 1:num_sigs
  snn = num2str(sigNum); % Signal Name Number
  curr_sig_name = [snb snn sne]; % Name of current signal
  curr_sig_string = string(curr_sig_name);
  sig_name_list{sigNum} = curr_sig_string; % Add that signals name to the list
end
toc

%% Saving Variables:
clear scaled_sig;
save('Saved_#0010_Data'); % Previously titled Preprocessed_clean_recordings_and_extracted_ground_truth
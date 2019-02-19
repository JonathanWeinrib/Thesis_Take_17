%%%% Function to preprocess Signal and do ENF estimation %%%%%
%%%% Jonathan Weinrib 
%%%% Date of Origination: 10/17/2018
% File Number: #0012

function [processed_ENF, processed_noise] = preprocess_noisy_rec(Signal_Array,noiseRec,freq_dev, mean_freq_list)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%Inputs
% Signal_Array = cell array of time-domain ENF signals with noise added
% noiseRec = cell array of time-domain purely noise signals
% freq_dev = frequency (in Hz) that we are allowed to be above or below the
% 50hz mean
% scaleOn = should we be scaling the recording. almost always yes. \
% mean_freq_list is a list containing which frequency to filter at for an
% individual recording. this is made in ENF_preprocess_0011

% maybe to add
% signal_name = string of the signal's name?
% save_dir = location to put saved data. e.g. 'Extracted_ENF_frame_4s/est_enf_'

%RETURNS
% processed_ENF = cell array of preprocessed ENF(+noise) recordings
% processed_noise = cell array of preprocessed pure noise recordings





% Arrays to hold the filtered ENF and pure_noise signals
num_ENF_recs = length(Signal_Array);
filtered_ENF = cell(1,num_ENF_recs);

num_noise_recs = length(noiseRec);
filtered_noise = cell(1,num_noise_recs);


% Arrays to hold scaled signals
scaled_ENF = cell(1,num_ENF_recs);
scaled_noise = cell(1,num_noise_recs);

% filter the ENF recordings

% Bandpassing the signals:
for numRec = 1:num_ENF_recs
    %%
    % Bandpass filter the current signal
    %% UMMMMM, but i need to know the grid letter then?
%     if grid_letter == 'A' | grid_letter == 'E'|grid_letter == 'I' % JON!! check that it's really grid E that has 60 hz mean)
%         mean_freq = 60;
%     else 
%         mean_freq = 50;
%     end
    %mean_freq = 50;
    mean_freq = mean_freq_list(numRec);


    %mean_freq = 50; % later change
    % cutoff frequencies in Hz
    f1 = mean_freq - freq_dev; f2 = mean_freq + freq_dev;
    fs = 1000; % Sampling frequency
    fn = fs/2; % Nyquist sampling rate
    w1 = f1/fn; w2 = f2/fn; %cutoff frequencies (nyquist)

    % for jonathan: do i need to downsample????

    %design the butterworth filter
    bttr_filt_ord = 3;
    [b,a] = butter(bttr_filt_ord,[w1 w2]);

    filtered_ENF{numRec} = filter(b,a,Signal_Array{numRec});
    
end % for numRec = 1:num_ENF_recs


%filter the noise recordings the same way
% note. I haven't accounted for the 50 to 60 hz difference for the noise
% recordings....... 
% update: 2/19/19: done i believe for above, but not below

% ask noah: i don't think this is too important, b/c white noise should be
% the same in all regions. oh actually it is, b/c we're putting it into a
% training set!
for numRec = 1:num_noise_recs
    

    mean_freq = mean_freq_list(numRec);

    %mean_freq = 50; % later change
    % cutoff frequencies in Hz
    f1 = mean_freq - freq_dev; f2 = mean_freq + freq_dev;
    fs = 1000; % Sampling frequency
    fn = fs/2; % Nyquist sampling rate
    w1 = f1/fn; w2 = f2/fn; %cutoff frequencies (nyquist)

    % for jonathan: do i need to downsample????

    %design the butterworth filter
    bttr_filt_ord = 3;
    [b,a] = butter(bttr_filt_ord,[w1 w2]);

    filtered_noise{numRec} = filter(b,a,noiseRec{numRec});
end 



%% Now, we normalize each signal.
% So, again, 
for numRec = 1:num_ENF_recs
    scaled_ENF{numRec} = filtered_ENF{numRec}./max(abs(filtered_ENF{numRec}));
end

%normalize the noise recordings the same way
for numRec = 1:num_noise_recs
    scaled_noise{numRec} = filtered_noise{numRec}./max(abs(filtered_noise{numRec}));
end

% Set the return variables
processed_ENF = scaled_ENF;
processed_noise = scaled_noise;


end % end ENF_STFT
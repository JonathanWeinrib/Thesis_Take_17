% Getting the ENF from the noised up signals
% Author: Jonathan Weinrib
% Originated 10/19/2018
% File Number: #0014

%% File explanation
% This file takes the noised up ENF recordings and processes them.
% That is to say, that it loads the both the noisy enf signals and the
% purely noise signals from file#0011, uses the function
% preprocess_noisy_rec (which is file number #0012) to bandpass filter and
% scale each signal, then uses the function ENF_STFT (which is file number
% #0013) to take the STFT of the signal and return both the estimated ENF
% (calculated as the frequency (i.e. frequency bin) with the maximum power
% for any given time window), as well as the entirety of the STFT as well.
% This file then places the STFT and estimation information in cell arrays
% titled 'processed_ENF_cellArray' and processed_noise_cellArray
% and saves them in something called 'Extracted_ENF_noisy_frame_16s/signal_list'


%% %% 
% This file 

%% 
clc, clear all, close all
%% Load Noisy ENF recordings
% Noised up actual Signals are in cell array noisy_sig_list

% Load purely noise recordings
% purely noise recordings are in cell array noise_pure_list
load('Noisy_Signals_#0011.mat');


%% First thing we need to do is preprocess our signals:
freq_dev = 5; % 5 hz frequency deviation allowed, so from 45-55 hz we will bandpass
% This function will return processed_ENF and processed_noise
tic
[processed_ENF, processed_noise] = preprocess_noisy_rec(ENF_plus_noise_Array,pure_Noise_Array,freq_dev,mean_freq_list)
% The returns of this function are the scaled and filtered time-domain
% signals of: 1) the Pure ENF + Noise, and 2) Pure noise
toc

%now we can clear our other variables to open up some space?
clear ENF_plus_noise_Array pure_Noise_Array
%%

% Define parameters for the STFT

% Parameters for all the recordings
frame_len = 16;
fs = 1000;
z_pad_fact = 8;
noverlap = 1024;
plotOn = false;
save_dir = 'Extracted_ENF_noisy_frame_16s/';

% parameters for the current recording
%signal = noisy_sig_list{2};
%signal_name = 'noisy_sig_list_1';

num_ENF_recs = length(processed_ENF);
num_noise_recs = length(processed_noise);

%% Now it's time to do the STFT on those files above, (and futurely,
%%% also some more feature extraction perhaps)


processed_ENF_cellArray = cell(1,num_ENF_recs);
for i = 1:num_ENF_recs
    signal = processed_ENF{i};
    signal_name = ['processed_ENF_' num2str(i)]
    % Call the STFT function
    [s,f,t,est_ENF] = ENF_STFT(signal,frame_len,fs,z_pad_fact,noverlap,plotOn,signal_name,save_dir);
    processed_ENF_cellArray{i} = {s,f,t,est_ENF};
end 
clear processed_ENF;
%% now for purely noise signals
% Note: jrw: go to GPU for this! 
processed_noise_cellArray = {};
for i = 1:num_noise_recs
    signal = processed_noise{i};
    signal_name = ['processed_noise_' num2str(i)]
    % Call the STFT function
    [s,f,t,est_ENF] = ENF_STFT(signal,frame_len,fs,z_pad_fact,noverlap,plotOn,signal_name,save_dir);
    processed_noise_cellArray{i} = {s,f,t,est_ENF};
end 
clear processed_noise;
%% Save my arrays
save('Extracted_ENF_noisy_frame_16s/signal_list','processed_ENF_cellArray','processed_noise_cellArray');
% We save it all in 'signal_list.mat'
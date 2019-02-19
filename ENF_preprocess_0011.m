%% Jonathan Weinrib
% Date of origination: 10/17/18
% Date of continuation, and date of origination with this name: 12/26/18
% File Number: #0011
%% Summary of this File: What I do
% This file creates a noised up version of ENF recordings.
% Basically, we previously read the files containing the clean ENF power
% signal (in .wav format), and converted them to .mat
% However, to do ML on these signals, we can't have a "clean" ENF signal,
% we need a noised up version. That's what this file does.

% To do this we need to know the power level of the original signal, so we
% can add in noise at a desired power level relative to the original signal
% level. So, the "clean", scaled, BP filtered ENF power signals are created
% and saved in file #0010, and then loaded here.

% In this file:
% First, load the scaled, bandpass filtered, ENF ground truth signal

% Second

%% Initial Preparation (not necessary)

% Clear Workspace
clc, clear all, close all


%% Making the signals with noise: 
%% First look at file number #0010 to get the appropriate signals.
% All we actually need from that file is sig_list, which contains our
% scaled, BP filtered ENF signals.
%load('Preprocessed_Clean_recordings_and_extracted_Ground_Truth','sig_list','num_sigs','fs');
load('Saved_#0010_Data');
%% Making the noise:
% We go back to the original signals and add noise to them to make noisy
% training data
% We also need to make pure noise to train our system on anyway

% The key element to doing this correctly is to have the appropriate power
% levels, as well as to know what our SNR's are for evaluation.
% Basically, it's not as easy as it seems

%% Part 1: Adding Noise to the clean ENF signal
%% Step 1: Finding the power in the signals:
% The power in each recording we compute below is take to be the 
% total power in the scaled, filtered ENF signals.

tic

% Define arrays to store power data from the ground truth ENF recordings.
sigPow_list   = zeros(1,num_sigs);
sigPowdB_list = zeros(1,num_sigs);

%sigPow_list = [];
%sigPowdB_list = [];


% To help find indices by which to locate in ENF_Files_Array
ind_list = [1];
for ind = 2:length(num_recordings_per_grid)
    ind_list(ind) = sum(num_recordings_per_grid(1:ind-1)) + 1;
end
%%

mean_freq_list = zeros(1,num_sigs)

% Finding the (Linear) power of each signal I have. For this we use the
% bandpassed scaled version of the signal.
for sigNum = 1:num_sigs % Loop through each signal we have

    % Find the grid that connects to this recording
    grid_ind = max(find(sigNum>=ind_list));
    % Find the corresponding recording number as well
    grid_rec_num = sigNum - ind_list(grid_ind) + 1;

    % Determine the Grid Number by which to access ENF_Files_Array
    % Determine the Recording Number by which to access ENF_Files_Array
    % Set the current signal
    
    % Note for power evaluation/confirmation: the recordings in
    % ENF_files_array (And so curr_sig here), are just the power from the
    % enf pure recordings in the [45 55] hz range (or 55 65)
    curr_sig = ENF_Files_Array{grid_ind,grid_rec_num};
    mean_freq_list(sigNum) = ENF_Filter_Freq_Array{grid_ind,grid_rec_num};

    % Find the total Power in our signal
    % in this case, bandpower just returns the average signal across the
    % whole signal. since we are baseband and have already filtered the rest
    % of the signal and noise out, this is totally okay!
    sig_power_lin = bandpower(curr_sig);

    % Save that power to a list of signal powers
    sigPow_list(sigNum) = sig_power_lin;

    % Find signal power in dB
    sig_power_db = pow2db(sig_power_lin);

    % Save that dB power to a list of signal powers
    sigPowdB_list(sigNum) = sig_power_db;
  
end % end saving signal powers

%% Step 2: Creating the Noised up ENF Signals
% Now that we have the power levels of each signal, we need to add noise in
% to each of those signals with a desired, variable, SNR. 
% We will use these recordings to make half of our dataset.


% Initialize cell to hold my different noisy signals
%noisy_sig_list = {};
ENF_plus_noise_Array = cell(1,num_sigs);

%Initialize arrays to hold the power level in the noised up signals
%noisySigPow_list = [];
noisySigPow_list = zeros(1,num_sigs);
%noisySigPowActual_list = [];
noisySigPowActual_list = zeros(1,num_sigs);
%noisySigPowdB_list = [];
noisySigPowdB_list = zeros(1,num_sigs);



% Loop through each signal in sig_list and create a noised up version of it
for sigNum = 1:num_sigs
    
    
    % Find the grid that connects to this recording
    grid_ind = max(find(sigNum>=ind_list));
    % Find the corresponding recording number as well
    grid_rec_num = sigNum - ind_list(grid_ind) + 1;
    
    
    % Set current signal
    curr_sig = ENF_Files_Array{grid_ind,grid_rec_num};
    
    % Add white gaussian noise to that signal
    % still need to figure out how the 'SNR', 'measured', and 'db' fields
    % affect things
    % create a noised up signal
    %NOTE: (JRW) later we will need to change this number to be more related to 
    % the actual SNR we want?
    % Note for power eval: noisy sig = white noise across all frequencies,
    % plus ENF in the [45 55] hz range
    noisy_sig = awgn(curr_sig,-26,'measured','db');
    % Add that noised up signal to our signal list
    %sigPowdB_list{sigNum} = noisy_sig;
    ENF_plus_noise_Array{sigNum} = noisy_sig;
    
    
    % noisySignal Powers
    % We are calculating this so that we can effectively utilize the powers
    % and SNR's correctly
    
    % Now we calculate the power in the noised up signal, over all
    % frequencies
    noisySig_power_lin = bandpower(noisy_sig);
    
    % Saving the actual noisy signal power, i.e. the signal power in the
    % desired bandwidth
    % NOTE: need to update, because right now i am filtering, but my stuff
    % is already baseband
    %Note: okay, b/c it's only used make noise? idk.
    
    mean_freq = ENF_Filter_Freq_Array{grid_ind,grid_rec_num};
    f1 = mean_freq - freq_deviation; f2 = mean_freq + freq_deviation;
    % QUESTION FOR NOAH::: do i need this?!
    %w1 = f1/fn; w2 = f2/fn; %cutoff frequencies (nyquist)
    
    noisySig_power_Actual = bandpower(noisy_sig,fs,[f1 f2]);
    % This saves the power of ENF + WGN in the band [f1 f2]. 
    % This is important, because when we have training data of just noise,
    % that's the power level we would want it to be at, e.g. equal to the
    % power of ENF + WGN in that band.
    
    % Save that power to a list of noisySignal powers
    noisySigPow_list(sigNum) = noisySig_power_lin;
    
    noisySigPowActual_list(sigNum) = noisySig_power_Actual;
    
 
    % Now we calculate the noisySignal power in dB
    noisySig_power_dB = pow2db(noisySig_power_lin);

    % Save that power to a list of noisySignal dB powers
    noisySigPowdB_list(sigNum) = noisySig_power_dB;
    
end % end created noisy signals and calculating the powers

%% Step 3: Calculating the SNR in each of the noisy Signals we have created

% Initialize array to store SNR list
%noisySig_SNR_list = [];
noisySig_SNR_list = zeros(1,num_sigs);
%SNR_Calculated_List = [];
SNR_Calculated_List = zeros(1,num_sigs);

% Initialize array to store the power of the noise channel
%noisy_channel_pow_list = [];
noisy_channel_pow_list = zeros(1,num_sigs);
%power_of_noise_alone_list = [];
power_of_noise_alone_list = zeros(1,num_sigs);

% Loop through each signal in sig_list
for sigNum = 1:num_sigs
    
    % Calculate the power in the noisy signal
    % you might think it's power_of_noisy_signal = noisySigPow_list(sigNum);
    % However, this is what happens when we spread the power of the
    % original ENF signal throughout the entire spectrum. However, in the desired range,
    % say [45 55] or something, the component from the noise will be much much smaller.
    % And, since we will be filtering the signal to a much smaller range during 
    % the preprocessing and estimation stages, the observed SNR will be
    % very different from that

    % So, to correctly calculate the power in the noisy signal, we have to
    % evaluate it in the upper and lower frequency limit range, f1 and f2,
    % where we have that f1 = mean_freq - freq_deviation; 
    %f2 = mean_freq + freq_deviation; given to us from our bandpass 
    % filtering stage.
    
    % NOTE JRW!: I have to change what mean_freq is apparently...
    power_of_noisy_signal = noisySigPowActual_list(sigNum);
    

    % maybe change this to be the bandpassed version of this? or rather,
    % check that the bandpassed version is the same first. TO DO JONATHAN
    power_of_original_signal = sigPow_list(sigNum);
    power_of_noise_alone = power_of_noisy_signal - power_of_original_signal;
    %powerdB_of_noise_alone = pow2db(power_of_noise_alone);
    % so power_of_noise_alone is the power of the noise in the [45 55] hz
    % range.

    % save the calculated power of the noise alone for later checking
    power_of_noise_alone_list(sigNum) = power_of_noise_alone;
    
 
    
    
    % We still calculate the power in the entire noisy signal so that we
    % may manually create a similar channel
    power_of_total_noisy_signal = noisySigPow_list(sigNum);
    power_of_noise_channel_alone = power_of_total_noisy_signal - power_of_original_signal;
    % so power_of_noise_channel_alone is the power of the noise channel,
    % across all frequencies
    
    
    % Save the power of the noise alone to a list so that we can make a
    % purely wgn channel with similar noise
    noisy_channel_pow_list(sigNum) = power_of_noise_channel_alone;
    
    
    
    %To calculate the SNR, we compare the power of the original signal
    % with the power of the noise component of the combined noisy signal.
    SNR_calculated = 10*log10(power_of_original_signal/power_of_noise_alone);
    
    % Save calculated SNR's to a list
    SNR_Calculated_List(sigNum) = SNR_calculated;
    
end

%Display the calculated SNRs
SNR_Calculated_List;
% for starters, it seems that using a value of -26 dB desired SNR when
% creating the noise get's us about -10 SNR in reality, which is a good
% start

%% Step 4: Creating recordings of just noise
% Now, we create recordings of just noise, and we try to make them similar
% to the actual recordings we have
% that is, we make a wgn sample to approximate the channel we have. I can
% either do that by taking away the WGN from the one I created (JONATHAN,
% this is probably a good thign to try at some point), or by creating
% another wgn channel with similar power level and all.

%pure_Noise_Array = {};
pure_Noise_Array = cell(1,num_sigs);

%power_of_pure_noise_inband = [];
power_of_pure_noise_inband = zeros(1,num_sigs);

%SNR_Calculated_List_pureNoise = [];
SNR_Calculated_List_pureNoise = zeros(1,num_sigs);



% Loop through each signal in sig_list
for sigNum = 1:num_sigs
    
    
    % Find the grid that connects to this recording
    grid_ind = max(find(sigNum>=ind_list));
    % Find the corresponding recording number as well
    grid_rec_num = sigNum - ind_list(grid_ind) + 1;    
    
    % Determine how long we want to make our purely noise recordings
    sig_length = length(ENF_Files_Array{grid_ind,grid_rec_num}); 
    
    desired_power = noisy_channel_pow_list(sigNum);
    % i make this the whole noisy channel, b/c when i use the wgn function
    % below, i can't just put it into a frequency range, i need to put it
    % across all frequencies, and so i calculated the power that the wgn
    % channel above created
    
    % Create white gaussian noise with length equal to the length of it's
    % corresponding recording, and with noise in the desried bandwidth
    % equal to the noise that we have in our ENF + noise signal
    
    % Create the signal
    noise_alone_signal = wgn(sig_length,1,desired_power,'linear');
    
    
    % Save the signal to a list. This constitutes the second half of 
    % the dataset.
    pure_Noise_Array{sigNum} = noise_alone_signal;
    % so noisy_pure_list is an array of just noise channels that
    % approximates the ENF+noise recordings we have. That is, it is equal
    % to just the noise component (statistically. it's not the exact same)
    % of the ENF+noise recordings, and has the same time length as them and
    
    % Note: the power of these signals i jsut created should :
    % this last bit is for testing, but basically i'm looking for
    % SNR_Calculated_List_pureNoise to be quite similar to
    % SNR_Calculated_List.
    % this would validate that i've done my creation of noise channel
    % correctly
    power_of_pure_noise_inband(sigNum) = bandpower(noise_alone_signal,fs,[f1 f2]);
    power_of_original_signal = sigPow_list(sigNum);
    SNR_Calculated_List_pureNoise(sigNum) = 10*log10(power_of_original_signal/power_of_pure_noise_inband(sigNum));
    
end

SNR_Calculated_List_pureNoise
SNR_Calculated_List
%power_of_pure_noise_inband
%power_of_noise_alone_list
% and it works!
% So i now successfully have my ENF + WGN time domain signal and my WGN
% time domain signal.
% Next step - feature extraction!!!

toc

%%
%For the moment, just gonna do it this way, and i can make it nicer later

%% Save signals 
save('Noisy_Signals_#0011','ENF_plus_noise_Array','pure_Noise_Array','mean_freq_list');

%% Next Step: STFT of signal:





%% Next Step Feature extraction
% Note: Skip for now
%% Next Step: Neural NETS!!!
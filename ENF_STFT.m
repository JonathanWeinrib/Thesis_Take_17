%%%% Function to preprocess Signal and do ENF estimation %%%%%
%%%% Jonathan Weinrib 
%%%% Date of Origination: 10/17/2018
% File Number: #0013

function [s,f,t,est_enf] = ENF_STFT(signal,frame_len,fs,z_pad_fact,noverlap,plotOn,signal_name,save_dir)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%Inputs
% signal = raw signal desired for the STFT
% frame_len = length of time frame used for window in seconds. can be 4,8,or 16 (default)
% fs = sampling frequency, Hz
% z_pad_fact = pad zeros by a factor of z_pad_fact
%noverlap = number of overlap points between shifting windows of the STFT. %default 1024?
% plotOn = do we want to plot the periodogram. Default (False)
% signal_name = string of the signal's name
%save_dir = location to put saved data. e.g. 'Extracted_ENF_frame_4s/est_enf_'

%RETURNS



% the final directory it will be saved to is for example:
%['Extracted_ENF_frame_4s/est_enf_' signal_name];


%%%% OUTPUT:
% We output the periodogram of the input signal as 
% s
% f
% t

% First step of STFT: define parameters

window_len = frame_len*fs; % We set the window size

% Now we set the number of points of the fft, nfft
% the number of points equals frame_len_time * fs, unless we pad zeros
% we are going to 
nfft = window_len * z_pad_fact;

%however, if grid signal is to long and the STFT is taking too long, we can split it up into
%sig_split_perc when we run the STFT
sig_split_perc = 1;


%now we get an estimation of the ENF signal by taking the STFT of our
%signal
[est_ENF,f,t] = spectrogram(signal(1:end/sig_split_perc),window_len,noverlap,nfft,fs);


% Scaling to ensure we have the appropriate power windowing
K = sum(hamming(window_len, 'periodic'))/window_len;
s = abs(est_ENF)/window_len/K;


% correction of the DC & Nyquist component
if rem(nfft, 2)                     % odd nfft excludes Nyquist point
    s(2:end, :) = s(2:end, :).*2;
else                                % even nfft includes Nyquist point
    s(2:end-1, :) = s(2:end-1, :).*2;
end

% convert amplitude spectrum to dB (min = -120 dB)
s = 20*log10(s + 1e-6); % jonathan: why 1e-6?

% plot the spectrogram
if plotOn
    figure;
    imagesc(t, f, s)
    set(gca,'YDir','normal')
    set(gca, 'FontName', 'Times New Roman', 'FontSize', 14)
    xlabel('Time, s')
    ylabel('Frequency, Hz')
    title('Amplitude spectrogram of the signal')
end


%estimate "enf" as the point witht the maximum power
num_times = length(t);
for i = 1:num_times

    [power_value, ind] = max(s(:,i));
    maxPower_freq = f(ind);
    est_enf(i) = maxPower_freq;
end
        
%         if plot_on
%             figure
%             plot(1:num_times,est_enf);
%             xlabel('Time (s)')
%             ylabel('Amplitude(dB)')
%             max_freq_value = max(est_enf)
%             min_freq_value = min(est_enf)
%             mean_freq_value = mean(est_enf)
%             median_freq_value = median(est_enf)
%             mode_freq_value = mode(est_enf)
%             freq_value_range = range(est_enf)
%         end


%% Save our values
save_directory = [save_dir signal_name];
save_ENF_loc = string(save_directory);
save(save_directory, 's', 'f', 't', 'est_enf');

end % end ENF_STFT
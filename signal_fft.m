function [t] = signal_fft(signal, fs)
% SIGNAL_FFT plots the input signal as a function of time and the Fast
% Fourier Transform of the same signal using fs as sampling frequency


ts = 1/fs;  % Sampling period
len = length(signal); % Length of the signal
t = 0:ts:(len-1)/fs; % Time vector

% Computing FFT signal
fft_1 = fft(signal - mean(signal)); %fft of the signal
fft_2 = fft_1(1:len/2); % one sided fft 
fft_2 = abs(fft_2); % obtain the magnitude of the fft and normalize the signal
frequencies = (0:len/2-1)*(fs/len); % scale the frequency axis


% Plotting Signal and FFT
figure("Name","Signal and FFT"); 
tiledlayout(2,1);

nexttile; plot(t, signal, 'r'); axis tight; grid on;
title("Original Signal"); xlabel("Time(s)"); ylabel("Amplitude"); legend("Signal");

nexttile; plot(frequencies, fft_2); axis tight; grid on;
title("One-Sided Amplitude FFT Spectrum"); xlabel("f(Hz)"); ylabel("Magnitude");
end
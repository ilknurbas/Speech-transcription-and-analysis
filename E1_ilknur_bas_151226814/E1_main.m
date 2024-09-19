clear all

curdir = fileparts(which('E1_main.m'));
addpath([curdir '/aux_scripts/']);


%% Task 6) Waveform reading and resampling
% Read waveform 251-136532-0016.flac
% and resample to 16 kHz (if not already; always ensure correct sample rate!)
%
% Required functions:   audioread()
%                       resample() 

[x,fs] = audioread([curdir '/data/251-136532-0016.flac']);
x;
fs; % 16000

% Read textGrid annotation of the signal
% (MODIFY THE PATH TO POINT TO YOUR ANNOTATION FILE FROM TASK 1)
% TG = readTextGrid([curdir '/251-136532-0016.TextGrid']);
TG = readTextGrid([curdir '/251-136532-0016.TextGrid']);

% Resample to 16 kHz (if not 16 kHz originally)

if(fs ~= 16000)
   x = resample(x, 16000, fs);
   fs = 16000;
end

x;
fs;  % 16000

%% Step 6) Windowing, waveform and spectrum plotting 
%--------------------------------------------------------------------------
% Extract a 20 ms hamming window from middle part of sound 'IH' of word "public", 
% and calculate its logarithmic magnitude spectrum.
% Create a plot with two panels, where the top panel shows the windowed 
% waveform and the bottom panel shows the corresponding log-spectrum.
%--------------------------------------------------------------------------
%
% Helpful functions:    fft(), hamming(), log10()
%                       figure(), plot(), subplot(), ylabel(), xlabel()
% 
%
% 
% Definitions:
%     D2a:  Signal magnitude in dB is defined as 
%           10*log10(linear_power_spectrum) 
%           20*log10(linear_amplitude_spectrum).
%
% Hints:
%    H2a:   you can find the sound timestamps using the reference annotation 
%           251-136532-0016_reference.TextGrid, read to TG struct earlier.
%
%    H2b:   Voiced speech spectrum typically has approx. 60 dB dynamic
%           range. 
%
%    H2c:  set(0,'defaulttextfontsize',16);  and
%          set(0,'defaultaxesfontsize',16);  can be used to set default
%          font sizes for figure axes and fonts.
%    H2d:  One-sided magnitude spectrum calculated from an even number of 
%          N samples has N/2+1 unique values from 0 to pi, and the rest are
%          symmetric to those. 
%    H2e:  Given a 16-kHz sampling rate, frequency band of the one-sided
%          magnitude spectrum should run from 0 to 8 kHz according to the
%          Nyquist criterion.



% Find onset and offset time of IH based on signal annotation
i = find(strcmp(TG.phones.val,'IH'),1);
t_start = TG.phones.start_at(i); % starting time of the sound
t_end = TG.phones.end_at(i);  % ending time of the sound
i;
t_start; %  0.9814;
t_end; %1.0398

% Use 20 ms Hamming window to extract a 20 ms segment from the signal.
% x is the sampled data hence change start and end time of the sound to
% samples
t_start_sampled = round(fs * t_start);
t_end_sampled = round(fs * t_end);
t_start_sampled; % 15703
t_end_sampled; % 16636
 
wl = 0.02 * fs; % window length in samples - 320
ww = hamming(wl); % define windowing function
wl;

% y =     % extract windowed signal y from full signal x
% found the starting point (in samples) of the windowed phone
% should be integer operands
window_start = round((t_start_sampled+t_end_sampled-wl)/2); 
window_start; % 16010
size(x); %148080           1
size(ww); % 320     1
y = x(window_start: (window_start + wl-1) ) .* ww;
size(y); % 320     1

% Calculate logarithmic magnitude spectrum (0-8000 Hz) to variable Y
% Remember to trim the mirrored (redundant) part from the spectrum 
% (hint H2d).

Y = abs(fft(y)); % untrimmed version comment it 
% Y = fft(y); % trimmed version uncomment it 
% Y = abs(Y(2:size(Y)/2+1));  % trimmed version uncomment it 
Y_db = 20*log10(Y); % db scale
Y_db;

figure(1);clf;
subplot(2,1,1);  
% Plot windowed signal waveform. Remember to add appropriate axis labels
% and units to the plot 
time = (0:length(y)-1) ./ fs;
length(fs);

plot(time, y)
title('windowed signal waveform')
xlabel("time(s)")
ylabel("amplitude")

subplot(2,1,2); 
% Plot log-magnitude spectrum. Remember to add appropriate axis labels and
% units to the plot. Note that speech dynamic range should be approx. 60
% dB.
length(Y_db);
length(fs);

freq_axis = linspace(0, fs/2, length(Y_db)); % untrimmed version comment it 
%freq_axis = linspace(0, fs/4, length(Y_db)); % trimmed version uncomment it 
plot(freq_axis, Y_db);
title('log-magnitude spectrum')
xlabel("frequency(Hz)")
ylabel("magnitude(dB)")





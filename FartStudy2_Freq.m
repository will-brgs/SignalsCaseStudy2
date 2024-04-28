%% Signals and Systems Case Study 2:
%% Introduction

% * Author:                   Mack LaRosa, Will Burgess, Leela Srinivas
% * Class:                    ESE 351
% * Date:                     Created 4/12/2024, Last Edited 4/28/2024

% main functional code block for case study. Includes design and up/down
% conversion for main objectives. Utilized poopfunc. Many functions of this
% file are then converted into ComSys for streamlined perforamnce analysis.
% Note that up/downscaling is not tested with the raised cosine pulse
% shape. Analysis of this shape is done with Performance_analysis_Rcos
%% Housekeeping
clear
close all
filepath = "C:\Users\Will\OneDrive - Washington University in St. Louis\. Signals & Systems\Case Study 2\Figures-body";
%% Generate Input Signal and Add Noise Factor, Bitrate = 1/Tp
Tp = 0.1; % Half pulse width
sample_period = Tp/40; % dt, pulse and recieve sample period
dt = sample_period;
sample_freq = 1/sample_period; % Frequency of pulse and recieve signal 

bit_rate = 1/(1 * Tp); %Fb, frequency of bits sent out
bit_period = 1/bit_rate; % Ts, Time between bits sent out
Ts = bit_period;
N = 20; % number of bits sent
sigma = 1;
%% Define Pulse Shapes
%w = linspace(-5, 5, (2*Ts)/dt+1); % frequency vector
t = -Ts:dt:Ts;
w = t.*1/dt;
% pulse_square_freq = ones(1,length(w));
% pulse_square_freq(1,1:round(1/4*length(pulse_square_freq))) = -1;
% pulse_square_freq(1,round(3/4*length(pulse_square_freq)):end) = -1;

numsymbols = 2;
pulse_rcos_time = rcosdesign(0.01,numsymbols,((length(t)-1)/numsymbols), 'normal');

pulse_sinc_time = sinc((2*t)/Ts);


pulse_rcos_freq = fftshift(fft(pulse_rcos_time));
pulse_sinc_freq = fftshift(fft(pulse_sinc_time));

fig1 = figure();
hold on
subplot(2,2,1), stem(t,(pulse_rcos_time), 'b')
xlabel('Time (s)'),ylabel('Amplitude'),title('Raised Cos Pulse in Time Domain')
subplot(2,2,2), plot(w,abs(pulse_rcos_freq),'r')
xlabel('Frequency (Hz)'),ylabel('Amplitude'),title('Raised Cosine Pulse in Frequency Domain')

subplot(2,2,3), stem(t,(pulse_sinc_time), 'b')
xlabel('Time (s)'),ylabel('Amplitude'),title('Sinc Pulse in Time Domain')
subplot(2,2,4), plot(w,abs(pulse_sinc_freq), 'r')
xlabel('Frequency (Hz)'),ylabel('Amplitude'),title('Sinc Pulse in Frequency Domain')


sgtitle('Pulse Shapes Utilized')
hold off
%% Use Function on Raised Cosine Pulse Shape
[~,~,xn, decoded, SNR] = poopFunc(abs(pulse_rcos_time), sigma, N);

fig2 = figure(); 
hold on
subplot(2,2,1),stem(xn,'b','filled')
xlabel('Index'),xlabel('Bit Index'), ylabel('Logic State')
subplot(2,2,2),stem(decoded, 'r','filled')
xlabel('Index'),xlabel('Bit Index'), ylabel('Logic State')
sgtitle('System Using Sinc Pulse')

% display performance results
error = (sum(xn ~= decoded))/(length(decoded)); 

disp('PERFORMANCE FOR SQUARE PULSE')
disp(['Bitrate: ' ,num2str(bit_rate), ' bits/second'])
disp(['Standard Deviation: ' , num2str(sigma)])
disp(['SNR: ' , num2str(SNR)])
disp(['Error: ' ,num2str(error*100),' percent'])
%% Use Function on Sinc Pulse
[~,~, xn, decoded, SNR] = poopFunc(abs(pulse_sinc_time), sigma, N);

subplot(2,2,3),stem(xn,'b','filled')
xlabel('Index'),ylabel('Amplitude'),title('Transmitted - Sinc Pulse')
subplot(2,2,4),stem(decoded, 'r','filled')
xlabel('Index'),ylabel('Amplitude'),title('Decoded - Sinc Pulse')
sgtitle('System Decoding Responses')
hold off

% display performance results
error = (sum(xn ~= decoded))/(length(decoded)); 

disp(' ');
disp('PERFORMANCE FOR SINC PULSE')
disp(['Bitrate: ' ,num2str(bit_rate), ' bits/second'])
disp(['Standard Deviation: ' , num2str(sigma)])
disp(['SNR: ' , num2str(SNR)])
disp(['Error: ' ,num2str(error*100),' percent'])

%% Up/down convert with 3 channel system - Sinc Pulse Shape
frequencies = [20,30,40];
[r,y, xn, ~, ~] = poopFunc(abs(pulse_sinc_time), sigma, N);

t_recieved = -Tp:dt:N * Ts + Tp -dt;
w_recieved = t_recieved .* 1/dt;

sinc_data_convolved = r;
band1_up = sinc_data_convolved .* cos(2*pi*frequencies(1)*t_recieved);
band2_up = sinc_data_convolved .* cos(2*pi*frequencies(2)*t_recieved);
band3_up = sinc_data_convolved .* cos(2*pi*frequencies(3)*t_recieved);

band1_up_fft = fft(band1_up);
band2_up_fft = fft(band2_up);
band3_up_fft = fft(band3_up);

fs = 1/dt; % sample frequency
Nfft = length(sinc_data_convolved); % length of fft
%f = 0:fs/Nfft:fs-fs/Nfft;
f = -(fs-fs/Nfft)/2:fs/Nfft:(fs-fs/Nfft)/2;
f = f * sample_period * 100;
%UPCONVERTING
fig3 = figure(); 
hold on
subplot(3,1,1),plot(f,abs(band1_up_fft))
xlabel('Frequency (Hz)'),ylabel('Amplitude'),title('Band 1, 20Hz')
subplot(3,1,2),plot(f,abs(band2_up_fft))
xlabel('Frequency (Hz)'),ylabel('Amplitude'),title('Band 2, 30Hz')
subplot(3,1,3),plot(f,abs(band3_up_fft))
xlabel('Frequency (Hz)'),ylabel('Amplitude'),title('Band 3, 40 Hz')
sgtitle('Three Upscaled bands - Sinc Pulse Shape')
hold off

fig4 = figure(); 
hold on
plot(f,abs(band1_up_fft),'r')
plot(f,abs(band2_up_fft),'m')
plot(f,abs(band3_up_fft),'b')
xlabel('Frequency (Hz)'),ylabel('Amplitude'),title('Merged Up-Converted Channels - Sinc Pulse Shape')
legend('20Hz Band','30Hz Band','40Hz Band')
hold off

upconverted_sinc = band1_up + band2_up + band3_up;

% DOWNCONVERTING
band1_down = upconverted_sinc .* cos(2*pi*frequencies(1)*t_recieved);
%band1_down = lowpass(band1_down, 6, sample_freq);
band1_down = conv(band1_down, pulse_sinc_time);
band1_down_fft = fft(band1_down);

band2_down = upconverted_sinc .* cos(2*pi*frequencies(2)*t_recieved);
%band2_down = lowpass(band2_down, 6, sample_freq);
band2_down = conv(band2_down, pulse_sinc_time);
band2_down_fft = fft(band2_down);

band3_down = upconverted_sinc .* cos(2*pi*frequencies(3)*t_recieved);
%band3_down = lowpass(band3_down, 6, sample_freq);
band3_down = conv(band3_down, pulse_sinc_time);
band3_down_fft = fft(band3_down);

fig5 = figure(); 
hold on
subplot(3,1,1),plot(real(band1_down_fft))
xlabel('Frequency (Hz)'),ylabel('Amplitude'),title('Band 1, 20Hz')
subplot(3,1,2),plot(real(band2_down_fft))
xlabel('Frequency (Hz)'),ylabel('Amplitude'),title('Band 2, 30Hz')
subplot(3,1,3),plot(real(band3_down_fft))
xlabel('Frequency (Hz)'),ylabel('Amplitude'),title('Band 3, 40 Hz')
sgtitle('Three Downconverted bands - Sinc Pulse Shape')
hold off

% Decode channnel 1
decoded_1 = zeros(1, N);
a = 0;
pulselen = length(pulse_sinc_time);
filterlen = length(band1_down);
factor = 1/(bit_rate * Tp); % find factor relating Ts and Tp, use that to modify pulselen
for i = pulselen + 1:(pulselen * factor + mod(factor, 2))/2:filterlen-pulselen * factor - 1
    a = a + 1;
    if(band1_down(i) > 0)
        decoded_1(a) = 1;
    else
       decoded_1(a) = -1;
    end
end
% Decode channel 2
decoded_2 = zeros(1, N);
a = 0;
pulselen = length(pulse_sinc_time);
filterlen = length(band2_down);
factor = 1/(bit_rate * Tp); % find factor relating Ts and Tp, use that to modify pulselen
for i = pulselen + 1:(pulselen * factor + mod(factor, 2))/2:filterlen-pulselen * factor - 1
    a = a + 1;
    if(band2_down(i) > 0)
        decoded_2(a) = 1;
    else
       decoded_2(a) = -1;
    end
end

% Decode channel 3
decoded_3 = zeros(1, N);
a = 0;
pulselen = length(pulse_sinc_time);
filterlen = length(band3_down);
factor = 1/(bit_rate * Tp); % find factor relating Ts and Tp, use that to modify pulselen
for i = pulselen + 1:(pulselen * factor + mod(factor, 2))/2:filterlen-pulselen * factor - 1
    a = a + 1;
    if(band3_down(i) > 0)
        decoded_3(a) = 1;
    else
       decoded_3(a) = -1;
    end
end

%Plot all 3 channels
fig6 = figure(); 

subplot(3,1,1), hold on
stem(xn, 'o', 'LineWidth', 1.5)
stem(decoded_1, 'x','LineWidth', 1)
xlabel('Bit Index'), ylabel('Logic State')
hold off, title('Channel 1'), legend('Transmitted Signal', 'Recieved Signal','location', 'east')

subplot(3,1,2), hold on
stem(xn, 'o','LineWidth', 1.5)
stem(decoded_2, 'x','LineWidth', 1)
xlabel('Bit Index'), ylabel('Logic State')
hold off, title('Channel 2'), legend('Transmitted Signal', 'Recieved Signal','location', 'east')

subplot(3,1,3),  hold on
stem(xn, 'o','LineWidth', 1.5)
stem(decoded_3, 'x', 'LineWidth', 1)
xlabel('Bit Index'), ylabel('Logic State')
title('Channel 3'), legend('Transmitted Signal', 'Recieved Signal','location', 'east')
sgtitle('Decoded Message Accuracy for all Three Chanels - Sinc Pulse Shape')
hold off

%% poopSending Test
disp("------------------------------------------------------------")

message = ('Will will wash your car!');
Sigma_vals = [0.7,1.2,2];
disp(["Original message: ", message])
for i = 1:3
disp("------------------------------------------------------------")
disp(['Transmision of sigma = ',num2str(Sigma_vals(i))]);
[~, ~, messageOut, SNR] = poopSend(pulse_sinc_time, Sigma_vals(i), message);
disp(['Calculated SNR: ',num2str(SNR)]);
end
%168 bit message

%% Save Figures
% exportgraphics(fig1, fullfile(filepath, 'Pulse_shapes.jpg'), 'resolution', 300);
% exportgraphics(fig2, fullfile(filepath, 'Decoding_responses.jpg'), 'resolution', 300);
% exportgraphics(fig3, fullfile(filepath, 'Upscaled_bands.jpg'), 'resolution', 300);
% exportgraphics(fig4, fullfile(filepath, 'Merged_bands.jpg'), 'resolution', 300);
% exportgraphics(fig5, fullfile(filepath, 'Downscaled_bands.jpg'), 'resolution', 300);
% exportgraphics(fig6, fullfile(filepath, 'Channel_accuracy.jpg'), 'resolution', 300);

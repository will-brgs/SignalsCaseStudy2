%% Housekeeping
close all
%% Generate Input Signal and Add Noise Factor, Bitrate = 1/Tp
Tp = 0.1; % Half pulse width
sample_period = Tp/50; % dt, pulse and recieve sample period
dt = sample_period;
sample_freq = 1/sample_period; % Frequency of pulse and recieve signal 

bit_rate = 1/(1 * Tp); %Fb, frequency of bits sent out
bit_period = 1/bit_rate; % Ts, Time between bits sent out
Ts = bit_period;

sigma = 0;
%% Define Pulse Shapes
%w = linspace(-5, 5, (2*Ts)/dt+1); % frequency vector
w = -1*Ts:dt:Ts;

% pulse_square_freq = ones(1,length(w));
% pulse_square_freq(1,1:round(1/4*length(pulse_square_freq))) = -1;
% pulse_square_freq(1,round(3/4*length(pulse_square_freq)):end) = -1;

pulse_rcos_freq = rcosdesign(.5, 3, length(w)+1);

% T = length(pulse_rcos_freq);
% pulse_rcos_freq = ones(1,length(w));
% 
% inside_cos_thing_left = ((pi*T)/beta) * (1:round(1/4*length(pulse_rcos_freq)) - (1-beta)/(2*T));
% pulse_rcos_freq(1, 1:round(1/4*length(pulse_rcos_freq))) = 0.5*(1+cos(inside_cos_thing_left));
% 
% inside_cos_thing_right = ((pi*T)/beta) * (1:round(1/4*length(pulse_rcos_freq)) - (1-beta)/(2*T));
% pulse_rcos_freq(1,round(3/4*length(pulse_square_freq)):end) = -1;
% 
% pulse_sinc_freq = sinc((20*w)/Ts);


pulse_rcos_time = ifftshift(ifft(pulse_rcos_freq));
pulse_sinc_time = ifftshift(ifft(pulse_sinc_freq));

figure, hold on
subplot(2,2,1), stem(pulse_rcos_freq)
xlabel('Frequency'),ylabel('Amplitude'),title('Raised Cosine Pulse in Frequency Domain')
subplot(2,2,2), plot(abs(pulse_rcos_time), 'r')
xlabel('Time'),ylabel('Amplitude'),title('Raised Cos Pulse in Time Domain')

subplot(2,2,3), stem(pulse_sinc_freq)
xlabel('Frequency'),ylabel('Amplitude'),title('Sinc Pulse in Frequency Domain')
subplot(2,2,4), plot(abs(pulse_sinc_time), 'r')
xlabel('Time'),ylabel('Amplitude'),title('Sinc Pulse in Time Domain')

sgtitle('Pulse Shapes Utilized')
hold off
%% Use Function on Square Pulse
[~, xn, decoded, SNR] = poopFunc(abs(pulse_rcos_time), sigma);

figure, hold on
subplot(2,1,1),stem(xn)
xlabel('Index'),ylabel('Amplitude'),title('Transmitted')
subplot(2,1,2),stem(decoded)
xlabel('Index'),ylabel('Amplitude'),title('Decoded')
sgtitle('System Using Sinc Pulse')
hold off

% display performance results
error = (sum(xn ~= decoded))/(length(decoded)); 

disp('PERFORMANCE FOR SQUARE PULSE')
disp(['Bitrate: ' ,num2str(bit_rate), ' bits/second'])
disp(['Standard Deviation: ' , num2str(sigma)])
disp(['SNR: ' , num2str(SNR)])
disp(['Error: ' ,num2str(error*100),' percent'])
%% Use Function on Sinc Pulse
[~, xn, decoded, SNR] = poopFunc(abs(pulse_sinc_time), sigma);

figure, hold on
subplot(2,1,1),stem(xn)
xlabel('Index'),ylabel('Amplitude'),title('Transmitted')
subplot(2,1,2),stem(decoded)
xlabel('Index'),ylabel('Amplitude'),title('Decoded')
sgtitle('System Using Square Pulse')
hold off

% display performance results
error = (sum(xn ~= decoded))/(length(decoded)); 

disp(' ');
disp('PERFORMANCE FOR SINC PULSE')
disp(['Bitrate: ' ,num2str(bit_rate), ' bits/second'])
disp(['Standard Deviation: ' , num2str(sigma)])
disp(['SNR: ' , num2str(SNR)])
disp(['Error: ' ,num2str(error*100),' percent'])

%% Up/down convert with 3 channel system
freq1 = 20;
freq2 = 30;
freq3 = 40;

band1 = pulse_sinc_time .* cos(2*pi*freq1*w);
band2 = pulse_sinc_time .* cos(2*pi*freq2*w);
band3 = pulse_sinc_time .* cos(2*pi*freq3*w);

fs = 1/dt; % sample frequency
Nfft = 1024; % length of fft
f = 0:fs/Nfft:fs-fs/Nfft;



% band1 = conv(pulse_sinc_freq, cos(2*pi*freq1*w)));
% band2 = conv(pulse_sinc_freq, cos(2*pi*freq2*w)));
% band3 = conv(pulse_sinc_freq, cos(2*pi*freq3*w)));

figure, hold on
subplot(3,1,1),plot(f,abs(fft(band1,Nfft)))
xlabel('Index'),ylabel('Amplitude'),title('Band 1, 20Hz')
subplot(3,1,2),plot(f,abs(fft(band2,Nfft)))
xlabel('Index'),ylabel('Amplitude'),title('Band 2, 30Hz')
subplot(3,1,3),plot(f,abs(fft(band3,Nfft)))
xlabel('Index'),ylabel('Amplitude'),title('Band 3, 40 Hz')
sgtitle('Three Upscaled bands')
hold off

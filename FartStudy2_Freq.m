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

sigma = 1;
%% Define Pulse Shapes
%w = linspace(-5, 5, (2*Ts)/dt+1); % frequency vector
t = -Ts:dt:Ts; % MUST BE 100 Long

% pulse_square_freq = ones(1,length(w));
% pulse_square_freq(1,1:round(1/4*length(pulse_square_freq))) = -1;
% pulse_square_freq(1,round(3/4*length(pulse_square_freq)):end) = -1;

numsymbols = 10;
pulse_rcos_time = rcosdesign(0.2,numsymbols,((length(t)-1)/numsymbols), 'normal');

pulse_sinc_time = sinc((10*t)/Ts);


pulse_rcos_freq = fftshift(fft(pulse_rcos_time));
pulse_sinc_freq = fftshift(fft(pulse_sinc_time));

figure, hold on
subplot(2,2,1), stem((pulse_rcos_time), 'b')
xlabel('Time'),ylabel('Amplitude'),title('Raised Cos Pulse in Time Domain')
subplot(2,2,2), plot(abs(pulse_rcos_freq),'r')
xlabel('Frequency'),ylabel('Amplitude'),title('Raised Cosine Pulse in Frequency Domain')

subplot(2,2,3), stem((pulse_sinc_time), 'b')
xlabel('Time'),ylabel('Amplitude'),title('Sinc Pulse in Time Domain')
subplot(2,2,4), plot(abs(pulse_sinc_freq), 'r')
xlabel('Frequency'),ylabel('Amplitude'),title('Sinc Pulse in Frequency Domain')


sgtitle('Pulse Shapes Utilized')
hold off
%% Use Function on Square Pulse
[~, xn, decoded, SNR] = poopFunc(abs(pulse_rcos_time), sigma);

figure, hold on
subplot(2,2,1),stem(xn,'b','filled')
xlabel('Index'),ylabel('Amplitude'),title('Transmitted - RCos Pulse')
subplot(2,2,2),stem(decoded, 'r','filled')
xlabel('Index'),ylabel('Amplitude'),title('Decoded - RCos Pulse')
sgtitle('System Using Sinc Pulse')

% display performance results
error = (sum(xn ~= decoded))/(length(decoded)); 

disp('PERFORMANCE FOR SQUARE PULSE')
disp(['Bitrate: ' ,num2str(bit_rate), ' bits/second'])
disp(['Standard Deviation: ' , num2str(sigma)])
disp(['SNR: ' , num2str(SNR)])
disp(['Error: ' ,num2str(error*100),' percent'])
%% Use Function on Sinc Pulse
[~, xn, decoded, SNR] = poopFunc(abs(pulse_sinc_time), sigma);

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

%% Up/down convert with 3 channel system
frequencies = [20,30,40];


band1 = pulse_sinc_time .* cos(2*pi*frequencies(1)*t);
band2 = pulse_sinc_time .* cos(2*pi*frequencies(2)*t);
band3 = pulse_sinc_time .* cos(2*pi*frequencies(3)*t);

fs = 1/dt; % sample frequency
Nfft = length(pulse_sinc_time); % length of fft
f = 0:fs/Nfft:fs-fs/Nfft;


% band1 = conv(pulse_sinc_freq, cos(2*pi*freq1*w)));
% band2 = conv(pulse_sinc_freq, cos(2*pi*freq2*w)));
% band3 = conv(pulse_sinc_freq, cos(2*pi*freq3*w)));

figure, hold on
subplot(3,1,1),plot(f,abs(fft(band1)))
xlabel('Index'),ylabel('Amplitude'),title('Band 1, 20Hz')
subplot(3,1,2),plot(f,abs(fft(band2)))
xlabel('Index'),ylabel('Amplitude'),title('Band 2, 30Hz')
subplot(3,1,3),plot(f,abs(fft(band3)))
xlabel('Index'),ylabel('Amplitude'),title('Band 3, 40 Hz')
sgtitle('Three Upscaled bands')
hold off

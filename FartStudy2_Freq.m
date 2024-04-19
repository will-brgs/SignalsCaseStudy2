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
w = -1*Ts:dt:1*Ts; % frequency vector

pulse_square_freq = ones(1,length(w));
pulse_square_freq(1,1:round(1/4*length(pulse_square_freq))) = -1;
pulse_square_freq(1,round(3/4*length(pulse_square_freq)):end) = -1;

pulse_sinc_freq = sinc((20*w)/Ts);


pulse_square_time = ifftshift(ifft(pulse_square_freq));
pulse_sinc_time = ifftshift(ifft(pulse_sinc_freq));

figure, hold on
subplot(2,2,1), stem(w,pulse_square_freq)
xlabel('Frequency'),ylabel('Amplitude'),title('Square Pulse in Frequency Domain')
subplot(2,2,2), plot(abs(pulse_square_time), 'r')
xlabel('Time'),ylabel('Amplitude'),title('Square Pulse in Time Domain')

subplot(2,2,3), stem(w,pulse_sinc_freq)
xlabel('Frequency'),ylabel('Amplitude'),title('Sinc Pulse in Frequency Domain')
subplot(2,2,4), plot(abs(pulse_sinc_time), 'r')
xlabel('Time'),ylabel('Amplitude'),title('Sinc Pulse in Time Domain')

sgtitle('Pulse Shapes Utilized')
hold off
%% Use Function on Square Pulse
[~, xn, decoded, SNR] = poopFunc(abs(pulse_square_time), sigma);

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

band1 = pulse_sinc_freq .* cos(1j*freq1*w);
band2 = pulse_sinc_freq .* cos(1j*freq2*w);
band3 = pulse_sinc_freq .* cos(1j*freq3*w);

band1 = conv(pulse_sinc_freq, cos(1j*freq1*w));
band2 = conv(pulse_sinc_freq, cos(1j*freq2*w));
band3 = conv(pulse_sinc_freq, cos(1j*freq3*w));

figure, hold on
subplot(3,1,1),plot(abs(band1))
xlabel('Index'),ylabel('Amplitude'),title('Band 1, 20Hz')
subplot(3,1,2),plot(abs(band2))
xlabel('Index'),ylabel('Amplitude'),title('Band 2, 30Hz')
subplot(3,1,3),plot(abs(band3))
xlabel('Index'),ylabel('Amplitude'),title('Band 3, 40 Hz')
sgtitle('Three Upscaled bands')
hold off
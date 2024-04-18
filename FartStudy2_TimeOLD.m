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
t = -1*Ts:dt:1*Ts; % time vector

pulse_square = ones(1,length(t));
pulse_square(1,1:round(1/4*length(pulse_square))) = 0;
pulse_square(1,round(3/4*length(pulse_square)):end) = 0;

pulse_sinc = sinc((2*t)/Ts);


ifft_pulse_square = ifft(pulse_square);
ifft_pulse_sinc = ifft(pulse_sinc);

figure, hold on
subplot(2,2,1), stem(t,pulse_square)
xlabel('Time'),ylabel('Amplitude'),title('Square Pulse Shape')
subplot(2,2,2), plot(abs(fftshift(fft_pulse_square)), 'r')
xlabel('Frequency(Hz)'),ylabel('Amplitude'),title('Square Pulse FFT')

subplot(2,2,3), stem(t,pulse_sinc)
xlabel('Time'),ylabel('Amplitude'),title('Sinc Pulse Shape')
subplot(2,2,4), plot(abs(fftshift(fft_pulse_sinc)), 'r')
xlabel('Frequency(Hz)'),ylabel('Amplitude'),title('Square Pulse FFT')

sgtitle('Pulse Shapes Utilized')
hold off
%% Use Function on Square Pulse
[~, xn, decoded, SNR] = poopFunc(pulse_square, sigma);

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
[~, xn, decoded, SNR] = poopFunc(pulse_sinc, sigma);

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
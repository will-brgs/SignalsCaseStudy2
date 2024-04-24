%% Housekeeping
clear
close all
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

% pulse_square_freq = ones(1,length(w));
% pulse_square_freq(1,1:round(1/4*length(pulse_square_freq))) = -1;
% pulse_square_freq(1,round(3/4*length(pulse_square_freq)):end) = -1;

numsymbols = 2;
pulse_rcos_time = rcosdesign(0.01,numsymbols,((length(t)-1)/numsymbols), 'normal');

pulse_sinc_time = sinc((2*t)/Ts);


pulse_rcos_freq = fftshift(fft(pulse_rcos_time));
pulse_sinc_freq = fftshift(fft(pulse_sinc_time));

figure, hold on
subplot(2,2,1), stem(t,(pulse_rcos_time), 'b')
xlabel('Time'),ylabel('Amplitude'),title('Raised Cos Pulse in Time Domain')
subplot(2,2,2), plot(t,abs(pulse_rcos_freq),'r')
xlabel('Frequency'),ylabel('Amplitude'),title('Raised Cosine Pulse in Frequency Domain')

subplot(2,2,3), stem(t,(pulse_sinc_time), 'b')
xlabel('Time'),ylabel('Amplitude'),title('Sinc Pulse in Time Domain')
subplot(2,2,4), plot(t,abs(pulse_sinc_freq), 'r')
xlabel('Frequency'),ylabel('Amplitude'),title('Sinc Pulse in Frequency Domain')


sgtitle('Pulse Shapes Utilized')
hold off
%% Use Function on Raised Cosine Pulse Shape
[~,~,xn, decoded, SNR] = poopFunc(abs(pulse_rcos_time), sigma);

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
[~,~, xn, decoded, SNR] = poopFunc(abs(pulse_sinc_time), sigma);

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
[r,y, xn, ~, ~] = poopFunc(abs(pulse_sinc_time), sigma);

t_recieved = -Tp:dt:N * Ts + Tp -dt;

sinc_data_convolved = r;
band1_up = sinc_data_convolved .* cos(2*pi*frequencies(1)*t_recieved);
band2_up = sinc_data_convolved .* cos(2*pi*frequencies(2)*t_recieved);
band3_up = sinc_data_convolved .* cos(2*pi*frequencies(3)*t_recieved);

band1_up_fft = fft(band1_up);
band2_up_fft = fft(band2_up);
band3_up_fft = fft(band3_up);

fs = 1/dt; % sample frequency
Nfft = length(sinc_data_convolved); % length of fft
f = 0:fs/Nfft:fs-fs/Nfft;

%UPCONVERTING
figure, hold on
subplot(3,1,1),plot(f,abs(band1_up_fft))
xlabel('Index'),ylabel('Amplitude'),title('Band 1, 20Hz')
subplot(3,1,2),plot(f,abs(band2_up_fft))
xlabel('Index'),ylabel('Amplitude'),title('Band 2, 30Hz')
subplot(3,1,3),plot(f,abs(band3_up_fft))
xlabel('Index'),ylabel('Amplitude'),title('Band 3, 40 Hz')
sgtitle('Three Upscaled bands - Sinc Pulse Shape')
hold off

figure, hold on
plot(f,abs(band1_up_fft),'r')
plot(f,abs(band2_up_fft),'m')
plot(f,abs(band3_up_fft),'b')
xlabel('Index'),ylabel('Amplitude'),title('Merged Up-Converted Channels - Sinc Pulse Shape')
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


figure, hold on
subplot(3,1,1),plot(real(band1_down_fft))
xlabel('Index'),ylabel('Amplitude'),title('Band 1, 20Hz')
subplot(3,1,2),plot(real(band2_down_fft))
xlabel('Index'),ylabel('Amplitude'),title('Band 2, 30Hz')
subplot(3,1,3),plot(real(band3_down_fft))
xlabel('Index'),ylabel('Amplitude'),title('Band 3, 40 Hz')
sgtitle('Three Downconverted bands - Sinc Pulse Shape')
hold off

% Downconverted Figure
% figure, hold on
% plot(f,abs((band1_down)),'r')
% plot(f,abs((band2_down)),'m')
% plot(f,abs((band3_down)),'b')
% xlabel('Index'),ylabel('Amplitude'),title('Merged Down-Converted Channels - Sinc Pulse Shape')
% legend('20Hz Band','30Hz Band','40Hz Band')
% hold off
% Decode 

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
figure
subplot(3,1,1), hold on
stem(xn, 'o', 'LineWidth', 1.5)
stem(decoded_1, 'x','LineWidth', 1)
hold off, title('Channel 1'), legend('Transmitted Signal', 'Recieved Signal','location', 'east')

subplot(3,1,2), hold on
stem(xn, 'o','LineWidth', 1.5)
stem(decoded_2, 'x','LineWidth', 1)
hold off, title('Channel 2'), legend('Transmitted Signal', 'Recieved Signal','location', 'east')

subplot(3,1,3),  hold on
stem(xn, 'o','LineWidth', 1.5)
stem(decoded_3, 'x', 'LineWidth', 1)
title('Channel 3'), legend('Transmitted Signal', 'Recieved Signal','location', 'east')
sgtitle('Decoded Message Accuracy for all Three Chanels')
hold off

%% poopSending Test
% [~, ~, messageOut, ~] = poopSend(pulse_sinc_time, 0, 'Abcde');
% binary = str2num(reshape(dec2bin('Abcde')',1,[])');
% 

% function [r,y, messageOut, SNR] = poopSend(pulse_shape, sigma, message)
pulse_shape = pulse_sinc_time;
sigma = 0;
message = 'Abcde';
binary = str2num(reshape(dec2bin(message)',1,[])');

Tp = 0.1; % Half pulse width
sample_period = Tp/40; % dt, pulse and recieve sample period
sample_freq = 1/sample_period; % Frequency of pulse and recieve signal 

bit_rate = 1/(1 * Tp); %Fb, frequency of bits sent out
bit_period = 1/bit_rate; % Ts, Time between bits sent out

pulse = pulse_shape;

N = length(binary);

%maxTime = N * bit_period;
binary(binary == 0) = -1

a = 0;
imp_train = zeros(1,N * bit_period * sample_freq);
for k = 1:length(imp_train)
    if mod(k - 1, sample_freq * bit_period) == 0
    a = a + 1;
    imp_train(k) = binary(a);
    else
    imp_train(k) = 0;    
    end
end

%sampleTimes = 0:sample_period:(N*bit_period)-sample_period;

y = conv(imp_train,pulse);

noise = sigma * max(y) * randn(1,length(y));
r = y + (noise * sigma);

filtered = conv(r,pulse);
decoded = zeros(1, N);

a = 0;
pulselen = length(pulse);
filterlen = length(filtered);

factor = 1/(bit_rate * Tp); % find factor relating Ts and Tp, use that to modify pulselen
% please name this something other than factor

for i = pulselen + 1:(pulselen * factor + mod(factor, 2))/2:filterlen-pulselen * factor - 1
    a = a + 1;
    if(filtered(i) > 0)
        decoded(a) = 1;
    else
       decoded(a) = 0;
    end
end

messageOut = char(bin2dec(num2str(reshape(decoded,7,[])')))';
disp(["Original message: ", message])
disp(["Decoded message: ", messageOut])
SNR = (sum(y.^2))/(sum(noise.^2));
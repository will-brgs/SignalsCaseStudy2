%% Housekeeping
close all
clear
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
pulse_rcos_time = rcosdesign(0.2,numsymbols,((length(t)-1)/numsymbols), 'sqrt');

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
[~,y, xn, ~, ~] = poopFunc(abs(pulse_sinc_time), sigma);

t_recieved = -Tp:dt:N * Ts + Tp -dt;

sinc_data_convolved = y;
band1 = sinc_data_convolved .* cos(2*pi*frequencies(1)*t_recieved);
band2 = sinc_data_convolved .* cos(2*pi*frequencies(2)*t_recieved);
band3 = sinc_data_convolved .* cos(2*pi*frequencies(3)*t_recieved);

fs = 1/dt; % sample frequency
Nfft = length(sinc_data_convolved); % length of fft
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
sgtitle('Three Upscaled bands - Sinc Pulse Shape')
hold off

figure, hold on
plot(f,abs(fft(band1)),'r')
plot(f,abs(fft(band2)),'m')
plot(f,abs(fft(band3)),'b')
xlabel('Index'),ylabel('Amplitude'),title('Merged Up-Converted Channels - Sinc Pulse Shape')
legend('20Hz Band','30Hz Band','40Hz Band')
hold off

upconverted_sinc = band1 + band2 + band3;

band1_down = upconverted_sinc .* sin(2*pi*frequencies(1)*t_recieved);
band1_down = lowpass(band1_down, 6, sample_freq);
%band1_down = conv(band1_down, pulse_sinc_time);

band2_down = upconverted_sinc .* sin(2*pi*frequencies(2)*t_recieved);
band2_down = lowpass(band2_down, 6, sample_freq);
%band2_down = conv(band2_down, pulse_sinc_time);


band3_down = upconverted_sinc .* sin(2*pi*frequencies(3)*t_recieved);
band3_down = lowpass(band3_down, 6, sample_freq);
%band3_down = conv(band3_down, pulse_sinc_time);

figure, hold on
subplot(3,1,1),plot(f,abs(fft(band1_down)))
xlabel('Index'),ylabel('Amplitude'),title('Band 1, 20Hz')
subplot(3,1,2),plot(f,abs(fft(band2_down)))
xlabel('Index'),ylabel('Amplitude'),title('Band 2, 30Hz')
subplot(3,1,3),plot(f,abs(fft(band3_down)))
xlabel('Index'),ylabel('Amplitude'),title('Band 3, 40 Hz')
sgtitle('Three Downconverted bands - Sinc Pulse Shape')
hold off
% Decode :3
decoded = zeros(1, N);

a = 0;
pulselen = length(pulse_sinc_time);
filterlen = length(band1_down);

factor = 1/(bit_rate * Tp); % find factor relating Ts and Tp, use that to modify pulselen

for i = pulselen + 1:(pulselen * factor + mod(factor, 2))/2:filterlen-pulselen * factor - 1
    a = a + 1;
    if(band1_down(i) > 0)
        decoded(a) = 1;
    else
       decoded(a) = -1;
    end
end

% figure, hold on
% plot(f,abs((band1_down)),'r')
% plot(f,abs((band2_down)),'m')
% plot(f,abs((band3_down)),'b')
% xlabel('Index'),ylabel('Amplitude'),title('Merged Down-Converted Channels - Sinc Pulse Shape')
% legend('20Hz Band','30Hz Band','40Hz Band')
% hold off
%% Up/down convert with 3 channel system -  Raised Cosine Shape
% frequencies = [20,30,40];
% [~,y, ~, ~, ~] = poopFunc(abs(pulse_rcos_time), sigma);
% 
% t_recieved = -Tp:dt:N * Ts + Tp -dt;
% 
% sinc_data_convolved = y;
% band1 = sinc_data_convolved .* cos(2*pi*frequencies(1)*t_recieved);
% band2 = sinc_data_convolved .* cos(2*pi*frequencies(2)*t_recieved);
% band3 = sinc_data_convolved .* cos(2*pi*frequencies(3)*t_recieved);
% 
% 
% figure, hold on
% subplot(3,1,1),plot(f,abs(fft(band1)))
% xlabel('Index'),ylabel('Amplitude'),title('Band 1, 20Hz')
% subplot(3,1,2),plot(f,abs(fft(band2)))
% xlabel('Index'),ylabel('Amplitude'),title('Band 2, 30Hz')
% subplot(3,1,3),plot(f,abs(fft(band3)))
% xlabel('Index'),ylabel('Amplitude'),title('Band 3, 40 Hz')
% sgtitle('Three Upscaled bands - RCos Pulse Shape')
% hold off
% 
% figure, hold on
% plot(f,abs(fft(band1)),'r')
% plot(f,abs(fft(band2)),'m')
% plot(f,abs(fft(band3)),'b')
% xlabel('Index'),ylabel('Amplitude'),title('Merged Up-Converted Channels - RCos Pulse Shape')
% legend('20Hz Band','30Hz Band','40Hz Band')
% hold off

%%
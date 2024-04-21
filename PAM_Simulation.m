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
%% Run Communication System For Different Sigma Values - Sinc Pulse Shape

frequencies = [20,30,40];
simlength = 100;
sigma = zeros(simlength,1);
SNR = zeros(simlength,1);
error_avg = zeros(simlength,1);
for i = 1:simlength
sigma(i) = 0.5 + (0.05 * i);
[SNR(i),error_1,error_2,error_3] = ComSys(pulse_sinc_time,frequencies,sigma(i));

% Calculate average error for each simulation in percent
error_avg(i) = (error_1 + error_2 + error_3)/3;
end
% p = polyfit(SNR, error_avg, 10);
% xfit = min(SNR):0.1:max(SNR);
model = @(b, SNR) b(1) * exp(b(2) * SNR);
initial_guess = [1, 0.1];
parameters = lsqcurvefit(model, initial_guess, SNR, error_avg);
xfit = linspace(min(SNR), max(SNR), 100);
yfit = model(parameters, xfit);

figure, hold on
scatter(SNR,error_avg, 'filled'), xlabel('SNR'),ylabel('Average Error Rate (%)')
plot((xfit),yfit, 'r', 'linewidth', 2)
title('Simulation of Various Sigma Values on Three Channels')
legend('Simulation Datapoints', 'Regression Line')
hold off

%% Run Communication System For Different Channel Bandwidth Values - Sinc Pulse Shape

frequency_ratios = [2,3,4];
simlength = 100;
sigma = 1;
bandwidth = zeros(simlength,1);
frequencies = zeros(length(frequency_ratios));
error_avg = zeros(simlength,1);

for i = 1:simlength
frequencies = frequency_ratios .* (100 - (0.95 * i));
bandwidth(i) = frequencies(2) - frequencies(1);

[~,error_1,error_2,error_3] = ComSys(pulse_sinc_time,frequencies,sigma);

% Calculate average error for each simulation in percent
error_avg(i) = (error_1 + error_2 + error_3)/3;
end
% p = polyfit(SNR, error_avg, 10);
% xfit = min(SNR):0.1:max(SNR);
model = @(b, SNR) b(1) * exp(b(2) * SNR);
initial_guess = [1, 0.1];
parameters = lsqcurvefit(model, initial_guess, bandwidth, error_avg);
xfit = linspace(min(bandwidth), max(bandwidth), 100);
yfit = model(parameters, xfit);

figure, hold on
scatter(bandwidth,error_avg, 'filled'), xlabel('Bandwidth (Hz)'),ylabel('Average Error Rate (%)')
plot((xfit),yfit, 'r', 'linewidth', 2)
title('Simulation of Various Sigma Values on Three Channels')
legend('Simulation Datapoints', 'Regression Line')
hold off
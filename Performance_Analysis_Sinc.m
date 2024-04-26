%% Introduction : 
% Provides all perfornace analysis seen for the sinc pulse shape. Generates
% both pulse shapes for easy comparrison. Performance analysis of SNR,
% bandwidth, and temporal pulse width compared. Also compares three channel
% performance to observe any potiential disparities. 
%% Housekeeping
clear
close all
options = optimoptions('lsqcurvefit', 'Display', 'off');
filepath = "C:\Users\Will\OneDrive - Washington University in St. Louis\. Signals & Systems\Case Study 2\Figures";
%% Generate Input Signal and Add Noise Factor, Bitrate = 1/Tp
Tp = 0.1; % Half pulse width
sample_period = Tp/40; % dt, pulse and recieve sample period
dt = sample_period;
sample_freq = 1/sample_period; % Frequency of pulse and recieve signal 

bit_rate = 1/(1 * Tp); %Fb, frequency of bits sent out
bit_period = 1/bit_rate; % Ts, Time between bits sent out
Ts = bit_period;
Ts = Tp;
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

fig1 = figure();
hold on
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
%% Run Communication System For Different Sigma Values, combined analysis - Sinc Pulse Shape

frequencies = [20,30,40];
simlength = 100;
internal_avg_length = 10;
sigma = zeros(simlength,1);
SNR = zeros(simlength,1);
error_avg = zeros(simlength,1);
N = 20;
for i = 1:simlength
internal_avg = zeros(internal_avg_length,1);
for j = 1:internal_avg_length
sigma(i) = 0.5 + (0.04 * i);
[SNR(i),error_1,error_2,error_3] = ComSys(pulse_sinc_time,frequencies,sigma(i),N);

% Calculate average error for each simulation in percent
internal_avg(j) = (error_1 + error_2 + error_3)/3;
end
error_avg(i) = sum(internal_avg) / length(internal_avg);
end

%generate regression line
model = @(b, SNR) b(1) * exp(b(2) * SNR);
initial_guess = [1, 0.1];
parameters = lsqcurvefit(model, initial_guess, SNR, error_avg);
xfit = linspace(min(SNR), max(SNR), 100);
yfit = model(parameters, xfit);

fig2 = figure();
hold on
scatter(SNR,error_avg, 'filled'), xlabel('SNR'),ylabel('Average Error Rate (%)')
grid on
plot((xfit),yfit, 'r', 'linewidth', 2)
title('Simulation of Various Sigma Values With Averaging of All Channelss')
legend('Simulation Datapoints', 'Regression Line')
hold off, grid off

%% Run Communication System For Different Sigma Values, channel analysis - Sinc Pulse Shape

frequencies = [20,30,40];
simlength = 100;
internal_avg_length = 10;
sigma = zeros(simlength,1);
SNR = zeros(simlength,1);
error_avg_1 = zeros(simlength,1);
error_avg_2 = zeros(simlength,1);
error_avg_3 = zeros(simlength,1);

N = 20;
for i = 1:simlength
internal_avg_1 = zeros(internal_avg_length,1);
internal_avg_2 = zeros(internal_avg_length,1);
internal_avg_3 = zeros(internal_avg_length,1);

for j = 1:internal_avg_length
sigma(i) = 0.5 + (0.04 * i);
[SNR(i),error_1,error_2,error_3] = ComSys(pulse_sinc_time,frequencies,sigma(i),N);
internal_avg_1(j) = error_1;
internal_avg_2(j) = error_2;
internal_avg_3(j) = error_3;
end
error_avg_1(i) = sum(internal_avg_1) / length(internal_avg_1);
error_avg_2(i) = sum(internal_avg_2) / length(internal_avg_2);
error_avg_3(i) = sum(internal_avg_3) / length(internal_avg_3);
end

ylimit = max(max(max(error_avg_1),max(error_avg_2)),max(error_avg_3))+1;

fig3 = figure();
%Channel 1
subplot(1,3,1)
scatter(SNR,error_avg_1, 'filled'), hold on
%regression line
model = @(b, SNR) b(1) * exp(b(2) * SNR);
initial_guess = [1, 0.1];
parameters = lsqcurvefit(model, initial_guess, SNR, error_avg_1);
xfit = linspace(min(SNR), max(SNR), 100);
yfit = model(parameters, xfit);
plot((xfit),yfit, 'r', 'linewidth', 2)
title('Channel 1')
xlabel('SNR'),ylabel('Average Error Rate (%)')
ylim([0,ylimit]);
legend('Simulation Datapoints', 'Regression Line')
hold off
%Channel 2
subplot(1,3,2)
scatter(SNR,error_avg_2, 'filled'), hold on
%regression line
model = @(b, SNR) b(1) * exp(b(2) * SNR);
initial_guess = [1, 0.1];
parameters = lsqcurvefit(model, initial_guess, SNR, error_avg_2);
xfit = linspace(min(SNR), max(SNR), 100);
yfit = model(parameters, xfit);
plot((xfit),yfit, 'r', 'linewidth', 2)
title('Channel 2')
xlabel('SNR'),ylabel('Average Error Rate (%)')
ylim([0,ylimit]);
legend('Simulation Datapoints', 'Regression Line')
hold off
%Channel 3
subplot(1,3,3)
scatter(SNR,error_avg_3, 'filled'), hold on
%regression line
model = @(b, SNR) b(1) * exp(b(2) * SNR);
initial_guess = [1, 0.1];
parameters = lsqcurvefit(model, initial_guess, SNR, error_avg_3);
xfit = linspace(min(SNR), max(SNR), 100);
yfit = model(parameters, xfit);
plot((xfit),yfit, 'r', 'linewidth', 2)
title('Channel 3')
xlabel('SNR'),ylabel('Average Error Rate (%)')
ylim([0,ylimit]);
legend('Simulation Datapoints', 'Regression Line')

sgtitle('Simulation of Various Sigma Values on Three Channels')
hold off

%% Run Communication System For Different Channel Bandwidth Values - Sinc Pulse Shape
frequency_ratios = [2,3,4];
simlength = 100;
sigma = 0.5;
internal_avg_length = 10;
bandwidth = zeros(simlength,1);
frequencies = zeros(length(frequency_ratios));
error_avg = zeros(simlength,1);
N = 20;

for i = 1:simlength
internal_avg = zeros(internal_avg_length,1);
for j = 1:internal_avg_length
frequencies = frequency_ratios .* (0.02* i);
bandwidth(i) = frequencies(2) - frequencies(1);

[~,error_1,error_2,error_3] = ComSys(pulse_sinc_time,frequencies,sigma,N);

internal_avg(j) = (error_1 + error_2 + error_3)/3;
end
% Calculate average error for each simulation in percent
error_avg(i) = sum(internal_avg) / length(internal_avg);
end
%generate regression line
p = polyfit(bandwidth, error_avg, 4);
xfit = min(bandwidth):0.1:max(bandwidth);
yfit = polyval(p, xfit);

fig4 = figure();
hold on
scatter(bandwidth,error_avg, 'filled'), xlabel('Bandwidth (Hz)'),ylabel('Average Error Rate (%)')
grid on
plot((xfit),yfit, 'r', 'linewidth', 2)
title('Simulation of Various Bandwidth Values on Three Channels')
legend('Simulation Datapoints', 'Regression Line')
hold off, grid off

%% Run Communication System For Different Temporal Pulse Durrations- Sinc Pulse Shape
frequencies = [20,30,40];
simlength = 100;
sigma = 0.5;
internal_avg_length = 10;
error_avg = zeros(simlength,1);
Ts_vals = zeros(simlength,1);
N = 20;
for i = 1:simlength
internal_avg = zeros(internal_avg_length,1);
%  Generate new pulse shape with new Ts

Ts = 0.1 + (i * 5e-3);
Ts_vals(i) = Ts;
pulse_sinc_time_Ts = sinc((2*t)/Ts);
for j = 1:internal_avg_length

[~,error_1,error_2,error_3] = ComSys(pulse_sinc_time_Ts,frequencies,sigma,N);

internal_avg(j) = (error_1 + error_2 + error_3)/3;
end
% Calculate average error for each simulation in percent
error_avg(i) = sum(internal_avg) / length(internal_avg);
end

%generate regression line
p = polyfit(Ts_vals, error_avg, 7);
xfit = linspace(min(Ts_vals),max(Ts_vals),simlength);
yfit = polyval(p, xfit);

fig5 = figure();
hold on
scatter(Ts_vals,error_avg, 'filled'), xlabel('Ts (s)'),ylabel('Average Error Rate (%)')
grid on
plot((xfit),yfit, 'r', 'linewidth', 2)
title('Simulation of Various Half-Pulse Temporal Widths on Three Channels')
legend('Simulation Datapoints', 'Regression Line', 'location', 'southeast')
hold off, grid off
%% Save Figures

exportgraphics(fig1, fullfile(filepath, 'analysis_pulseshapes.jpg'), 'resolution', 300);
exportgraphics(fig2, fullfile(filepath, 'analysis_sigma_merged.jpg'), 'resolution', 300);
exportgraphics(fig3, fullfile(filepath, 'analysis_sigma_independent.jpg'), 'resolution', 300);
exportgraphics(fig4, fullfile(filepath, 'analysis_bandwidth.jpg'), 'resolution', 300);
exportgraphics(fig5, fullfile(filepath, 'analysis_pulsewidth.jpg'), 'resolution', 300);

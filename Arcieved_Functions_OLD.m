%% Signals and Systems Matlab Homework #6
%% Introduction

%NOTE: THE FOLLOWING CODE IS NOT PART OF THE SUBMISSION BUT ACTS AS A
%HOLDING PLACE FOR PREVIOUSLY OMITTED CODE.
% * Author:                   Mack LaRosa, Will Burgess, Leela Srinivas
% * Class:                    ESE 351
% * Date:                     Created 4/12/2024, Last Edited 4/25/2024
%% Housekeeping
close all
clear
clc
%% Part  1: Create Binary PAM System
%% Generate Input Signal and Add Noise Factor, Bitrate = 1/Tp
Tp = 0.1; % Half pulse width
sample_period = Tp/50; % dt, pulse and recieve sample period
dt = sample_period;
sample_freq = 1/sample_period; % Frequency of pulse and recieve signal 

bit_rate = 1/(1 * Tp); %Fb, frequency of bits sent out
bit_period = 1/bit_rate; % Ts, Time between bits sent out
Ts = bit_period;

%  Generate Message
rect = ones(1,50);

pulse = 2 * conv(rect, rect);
%%
% replace pulse with sinc
t = -2*Ts:dt:2*Ts; % time vector
pulse = sinc(t/Ts);


pulse_fft = fftshift(pulse);

N = 20;

maxTime = N * bit_period;

xn = 2 * ((rand(1, N) > 0.5) - 0.5);
a = 0;
imp_train = zeros(1,N * bit_period * sample_freq);
for k = 1:length(imp_train)
    if mod(k - 1, sample_freq * bit_period) == 0
    a = a + 1;
    imp_train(k) = xn(a);
    else
    imp_train(k) = 0;    
    end
end
sampleTimes = 0:sample_period:(N*bit_period)-sample_period;


% Mack is a Bitchass Hoe: finish generation with signal


y = conv(imp_train,pulse);
% figure, subplot (2,1,1),plot(y)
% subplot(2,1,2),stem(xn)

sigma = 0;
noise = sigma * max(y) * randn(1,length(y));
r = y + (noise * sigma);

%figure,plot(r)
%% Part f:  Create various outputs

% i, p(t), p(w)
figure, hold on
subplot(2,1,1);
stem(pulse)
ylabel('Amplitude')
xlabel('Index')
title('Pulse P(t), Bitrate = 1/Tp')

subplot(2,1,2)
plot(pulse_fft)
ylabel('Amplitude')
xlabel('Frequency (Hz)')
title('Pulse Sprectum P(w), Bitrate = 1/Tp')
hold off

times = linspace(0, maxTime + 2 * Tp, length(y));

% ii, y(t)
figure, hold on
plot(times, y)
xlim([0 maxTime + 2 * Tp]);
ylabel('Amplitude')
xlabel('Time (s)')
title('Transmitted Signal y(t), Bitrate = 1/Tp')
hold off

% iii, r(t)
figure, hold on
plot(times,r)
xlim([0 maxTime + 2 * Tp]);
ylabel('Amplitude')
xlabel('Time (s)')
title('Recieved Signal r(t), Bitrate = 1/Tp')
hold off

% iv, Sent Messave vs. decoded message
times_sent = 0:bit_period:maxTime-bit_period;
figure, hold on
subplot(2,1,1)
stem(times_sent,xn)
ylabel('Amplitude')
xlabel('Time (s)')
title('Transmitted Signal xn, Bitrate = 1/Tp')
subplot(2,1,2)

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
       decoded(a) = -1;
    end
end

stem(times_sent,decoded);
ylabel('Amplitude')
xlabel('Time (s)')
title('Decoded Signal r(t), Bitrate = 1/Tp')
hold off

error = (sum(xn ~= decoded)/length(decoded)); 
SNR = (sum(y.^2))/(sum(noise.^2));

disp(['Bitrate: ' ,num2str(bit_rate), ' bits/second'])
disp(['Standard Deviation: ' , num2str(sigma)])
disp(['SNR: ' , num2str(SNR)])
disp(['Error: ' ,num2str(error),' percent'])

%% Generate Variable length Sinc Function

%% demo - Sinc pulse shape
Ts = .1; % symbol period (rate 1/Ts)
dt = .01; % sample period
t = -5*Ts:dt:5*Ts; % time vector
x = sinc(t/Ts); % define sinc, note Matlab convention sinc(x) = sin(pi*x)/(pi*x)
figure
subplot(2,1,1), plot(t,x)
xlabel('time (s)'), ylabel('x(t)'), title('Truncated sinc')
fs = 1/dt; % sample frequency
Nfft = 1024; % length of fft
f = [0:fs/Nfft:fs-fs/Nfft];
subplot(2,1,2), plot(f,abs(fft(x,Nfft)))
xlabel('frequency (Hz)'), ylabel('|X(j\omega)|')

pulse_sinc = x;

%%
[r, xn, decoded] = poopFunc(pulse, 1);
figure, plot(r)

figure, subplot(2,1,1),stem(xn)
subplot(2,1,2),stem(decoded)

%% Run Communication System For Different Bit Sent Values - Sinc Pulse Shape
frequencies = [20,30,40];
simlength = 100;
sigma = 0.5;
internal_avg_length = 10;
error_avg = zeros(simlength,1);
Ts = 10;

N_vals = zeros(simlength,1);

for i = 1:simlength
internal_avg = zeros(internal_avg_length,1);
%  Generate new pulse shape with new Ts

N = 9 + (i);
N_vals(i) = N;
for j = 1:internal_avg_length

[~,error_1,error_2,error_3] = ComSys(pulse_sinc_time,frequencies,sigma,N);

internal_avg(j) = (error_1 + error_2 + error_3)/3;
end
% Calculate average error for each simulation in percent
error_avg(i) = sum(internal_avg) / length(internal_avg);
end

%generate regression line
p = polyfit(N_vals, error_avg, 6);
xfit = linspace(min(N_vals),max(N_vals),simlength);
yfit = polyval(p, xfit);

figure, hold on
scatter(N_vals,error_avg, 'filled'), xlabel('Number of Bits Transmitted'),ylabel('Average Error Rate (%)')
grid on
plot((xfit),yfit, 'r', 'linewidth', 2)
title('Simulation of Various Numbers of Bits Transmitted Over Three Channels')
legend('Simulation Datapoints', 'Regression Line', 'location', 'southeast')
hold off, grid off

%%%

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
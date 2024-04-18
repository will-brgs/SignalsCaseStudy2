%% Signals and Systems Matlab Homework #6
%% Introduction
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
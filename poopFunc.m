
function [r,y, xn,  decoded, SNR] = poopFunc(pulse_shape, sigma, N)
Tp = 0.1; % Half pulse width
sample_period = Tp/40; % dt, pulse and recieve sample period
sample_freq = 1/sample_period; % Frequency of pulse and recieve signal 

bit_rate = 1/(1 * Tp); %Fb, frequency of bits sent out
bit_period = 1/bit_rate; % Ts, Time between bits sent out

pulse = pulse_shape;

%maxTime = N * bit_period;

xn = 2 * ((rand(1, N) > 0.5) - 0.5);
a = 0;
imp_train = zeros(1,round(N * bit_period * sample_freq));
for k = 1:length(imp_train)
    if mod(k - 1, sample_freq * bit_period) == 0
    a = a + 1;
    imp_train(k) = xn(a);
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

for i = pulselen + 1:(pulselen * factor + mod(factor, 2))/2-1:filterlen-pulselen * factor - 1
    a = a + 1;
    if(filtered(i) > 0)
        decoded(a) = 1;
    else
       decoded(a) = -1;
    end
end
SNR = (sum(y.^2))/(sum(noise.^2));
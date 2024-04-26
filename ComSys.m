function [SNR,error_1,error_2,error_3] = ComSys(pulse_shape_time, frequencies, sigma, N, Tp)
% Function Takes input arguments of factors we wish to analize for
% performance analysis. Acts as a way to work through the designed up and
% downscaling seen in fartstudy2_freq. Outputs the errors seen in each
% channel as well as total SNR. Note that the message is identical accross
% all 3 channels, but this does not imrpove performance whatsoever.

%   Detailed explanation goes here
sample_period = Tp/40; % dt, pulse and recieve sample period
dt = sample_period;

bit_rate = 1/(1 * Tp); %Fb, frequency of bits sent out
bit_period = 1/bit_rate; % Ts, Time between bits sent out
Ts = bit_period;

[r, y,xn, ~, SNR] = poopFunc(pulse_shape_time, sigma,N);

t_recieved = -Tp:dt:N * Ts + Tp -dt;

data_convolved = r;
band1_up = data_convolved .* cos(2*pi*frequencies(1)*t_recieved);
band2_up = data_convolved .* cos(2*pi*frequencies(2)*t_recieved);
band3_up = data_convolved .* cos(2*pi*frequencies(3)*t_recieved);

upconverted_sinc = band1_up + band2_up + band3_up;

% DOWNCONVERTING
band1_down = upconverted_sinc .* cos(2*pi*frequencies(1)*t_recieved);
band1_down = conv(band1_down, pulse_shape_time);

band2_down = upconverted_sinc .* cos(2*pi*frequencies(2)*t_recieved);
band2_down = conv(band2_down, pulse_shape_time);

band3_down = upconverted_sinc .* cos(2*pi*frequencies(3)*t_recieved);
band3_down = conv(band3_down, pulse_shape_time);


% Decode channnel 1
decoded_1 = zeros(1, N);
a = 0;
pulselen = length(pulse_shape_time);
filterlen = length(band1_down);
factor = 1/(bit_rate * Tp); % find factor relating Ts and Tp, use that to modify pulselen
for i = pulselen + 1:(pulselen * factor + mod(factor, 2))/2-1:filterlen-pulselen * factor - 1
    a = a + 1;
    if(band1_down(i) > 0)
        decoded_1(a) = 1;
    else
       decoded_1(a) = -1;
    end
end
error_1 = 100*(sum(xn ~= decoded_1))/(length(decoded_1)); 


% Decode channel 2
decoded_2 = zeros(1, N);
a = 0;
pulselen = length(pulse_shape_time);
filterlen = length(band2_down);
factor = 1/(bit_rate * Tp); % find factor relating Ts and Tp, use that to modify pulselen
for i = pulselen + 1:(pulselen * factor + mod(factor, 2))/2-1:filterlen-pulselen * factor - 1
    a = a + 1;
    if(band2_down(i) > 0)
        decoded_2(a) = 1;
    else
       decoded_2(a) = -1;
    end
end
error_2 = 100*(sum(xn ~= decoded_2))/(length(decoded_2)); 

% Decode channel 3
decoded_3 = zeros(1, N);
a = 0;
pulselen = length(pulse_shape_time);
filterlen = length(band3_down);
factor = 1/(bit_rate * Tp); % find factor relating Ts and Tp, use that to modify pulselen
for i = pulselen + 1:(pulselen * factor + mod(factor, 2))/2-1:filterlen-pulselen * factor - 1
    a = a + 1;
    if(band3_down(i) > 0)
        decoded_3(a) = 1;
    else
       decoded_3(a) = -1;
    end
end
error_3 = 100*(sum(xn ~= decoded_3))/(length(decoded_3)); 


end


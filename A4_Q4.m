clc; close all; clear all;

%Question 4
%Jinseng Vanderkloot
%101031534

C=0.25;
G1 = 1;
G2 = 1/2;
G3 = 1/8.5;
G4 = 1/0.1;
G0 = 1/1000;
L = 0.2;
alpha = 100;

%% C and G Matrix 
%F = V1, V2, V3, V4, V5, Vin, I4, IL 

CVec = [-C, C, 0, 0, 0, 0, 0, 0; 
      C, -C, 0, 0, 0, 0, 0, 0;
      0, 0, 0, 0, 0, 0, 0, 0;
      0, 0, 0, 0, 0, 0, 0, 0;
      0, 0, 0, 0, 0, 0, 0, 0;
      0, 0, 0, 0, 0, 0, 0, 0;
      0, 0, 0, 0, 0, 0, 0, 0;
      0, 0, 0, 0, 0, 0, 0, -L;];

GVec = [-G1, G1, 0, 0, 0, -1, 0, 0;
        G1, -G1-G2, 0, 0, 0, 0, 0, -1;
        0, 0, -G3, 0, 0, 0, 0, 1;
        0, 0, 0, -G4, G4, 0, 1, 0; 
        0, 0, 0, G4, -G4-G0, 0, 0, 0;
        1, 0, 0, 0, 0, 0, 0, 0; 
        0, 0, -alpha*G3, 1, 0, 0, 0, 0;
        0, 1, -1, 0, 0, 0, 0, 0;];

FVec = [0;0;0;0;0;1;0;0];

vin = linspace(-10,10);
outputV = zeros(size(vin));
freq = linspace(0.01, 100);

for a = 1:size(freq,2)
    w =(2*pi*freq(a));
    jw = 1i*w; % 1i makes value complex 1 -> 0+1j 
    V = (GVec+jw*CVec)\FVec; 
    outputV(a) = V(5);
end

%Plot Gain 
figure('Name','Circuit Frequency Response');
semilogx(freq, 20*log(abs(outputV)));
title('Frequency Response');
xlabel('Frequency  (Hz)');
ylabel('Gain (dB)');

% Question 4a) By inpection of the frequency response of the circuit, the
% result is a low pass filter with the cut off frequency being around 2Hz.

%Question 4b) The circuit has 40-46dB of gain for near DC values (any
%freqeuncy less than 2Hz). There is a slight slope for the gain so it is
%not ideal. 

%Question 4c) 
%   C(dV/dt) + GV = F
%   C[(Vn - Vn-1)/ delta_t] + G Vn = F(tn)
%   C (Vn/delta_t) - C (Vn-1 / delta_t) + G Vn = F(tn)
%   Vn = [C(Vn-1 / delta_t) + F(tn)] / (Vn/delta_t)

%% Question 4d - Step Responce 

t = 0:0.001:1;                          %Time 1000 steps 0 to 1 
F1 = [0;0;0;0;0;0;0;0];                 % initial values
vout = zeros(size(GVec,1),size(t,2));  % output vector
vin= zeros(1,size(t,2));               % input vector
vin(1,1) = 0;                           %initial input value 

%Loop through time 
for n=1:size(t,2)-1                      
    vout(:,n) = F1; 
    if n<30                        %Produce pulse
        vin(n+1)=0;
    else
        vin(n+1)=1;
    end
    F2 = FVec * vin(n+1);            % Capturing next input
    F1 = ((CVec./0.001) + GVec)\((CVec./0.001)*F1 + F2);
end
vout(:,n+1) = F1;

figure("Name","Step Response")
plot(t, vin,'r');
hold on
plot(t, vout(5,:),'b');
title('Step Response');
xlabel('Time (s)');
ylabel('Voltage (V)');
legend({'Input','Output'});

figure("Name","Step Fourier Response Spectrum")
semilogy(abs(fftshift(fft(vout(5,:)))))
xlabel('Frequency (Hz)');
ylabel('Amplitude');
xlim([250 750]);

t = 0:0.002:1;                          %Time 1000 steps 0 to 1 
F1 = [0;0;0;0;0;0;0;0];                 % initial values
vout = zeros(size(GVec,1),size(t,2));  % output vector
vin= zeros(1,size(t,2));               % input vector
vin(1,1) = 0;                           %initial input value 

%Loop through time 
for n=1:size(t,2)-1                      
    vout(:,n) = F1; 
    if n<30                        %Produce pulse
        vin(n+1)=0;
    else
        vin(n+1)=1;
    end
    F2 = FVec * vin(n+1);            % Capturing next input
    F1 = ((CVec./0.002) + GVec)\((CVec./0.002)*F1 + F2);
end
vout(:,n+1) = F1;

figure("Name","Step Response - Double Time Step")
plot(t, vin,'r');
hold on
plot(t, vout(5,:),'b');
title('Step Response - Double Time Step');
xlabel('Time (s)');
ylabel('Voltage (V)');
legend({'Input','Output'});

figure("Name","Step Fourier Response Spectrum - Larger Time Step")
semilogy(abs(fftshift(fft(vout(5,:)))))
xlabel('Frequency (Hz)');
ylabel('Amplitude');
xlim([0 500])

%% Question 4d - Sine wave @ Differnt frequency 

t = 0:0.001:1;
freq2 = 1/0.03; % Frequency 33Hz
F1 = [0;0;0;0;0;0;0;0];                 % initial values
vout2 = zeros(size(GVec,1),size(t,2)); % Vector to hold outputs
vin2 = zeros(1,size(t,2));          % Vector to hold inputs
vin2(1,1) = sin(2*pi*freq2*t(1));

%Loop through time 
for n=1:size(t,2)-1
    vout2(:,n) = F1;
    vin2(1,n+1) = sin(2*pi*freq2*t(n+1));
    F2 = FVec * vin2(1,n+1);
    F1 = ((CVec./0.001) + GVec)\((CVec./0.001)*F1 + F2);
end
vout2(:,n+1) = F1;

figure("Name","Sine Wave Response 33Hz")
plot(t, vin2,'r');
hold on
plot(t, vout2(5,:),'b');
title('Sine Wave Response 33Hz');
xlabel('Time (s)');
ylabel('Voltage (V)');
legend({'Input','Output'});

figure("Name","Sine Wave Fourier Response Spectrum 33Hz")
semilogy(abs(fftshift(fft(vout2(5,:)))))
xlabel('Frequency (Hz)');
ylabel('Amplitude');
xlim([250 750])

t = 0:0.002:1;
freq2 =1/0.03; % Frequency 33Hz
F1 = [0;0;0;0;0;0;0;0];                 % initial values
vout2 = zeros(size(GVec,1),size(t,2)); % Vector to hold outputs
vin2 = zeros(1,size(t,2));          % Vector to hold inputs
vin2(1,1) = sin(2*pi*freq2*t(1));

%Loop through time 
for n=1:size(t,2)-1
    vout2(:,n) = F1;
    vin2(1,n+1) = sin(2*pi*freq2*t(n+1));
    F2 = FVec * vin2(1,n+1);
    F1 = ((CVec./0.002) + GVec)\((CVec./0.002)*F1 + F2);
end
vout2(:,n+1) = F1;

figure("Name","Sine Wave Response 33Hz - Double Time Step")
plot(t, vin2,'r');
hold on
plot(t, vout2(5,:),'b');
title('Sine Wave Response 33Hz - Double Time Step');
xlabel('Time (s)');
ylabel('Voltage (V)');
legend({'Input','Output'});

figure("Name","Sine Wave Fourier Response Spectrum 33Hz")
semilogy(abs(fftshift(fft(vout2(5,:)))))
xlabel('Frequency (Hz)');
ylabel('Amplitude');
xlim([0 500])

t = 0:0.001:1;
freq2 = 1; % Frequency in 1Hz
F1 = [0;0;0;0;0;0;0;0];                 % initial values
vout2 = zeros(size(GVec,1),size(t,2)); % Vector to hold outputs
vin2 = zeros(1,size(t,2));          % Vector to hold inputs
vin2(1,1) = sin(2*pi*freq2*t(1));

%Loop through time 
for n=1:size(t,2)-1
    vout2(:,n) = F1;
    vin2(1,n+1) = sin(2*pi*freq2*t(n+1));
    F2 = FVec * vin2(1,n+1);
    F1 = ((CVec./0.001) + GVec)\((CVec./0.001)*F1 + F2);
end
vout2(:,n+1) = F1;

figure("Name","Sine Wave Response 1Hz")
plot(t, vin2,'r');
hold on
plot(t, vout2(5,:),'b');
title('Sine Wave Response');
xlabel('Time (s)');
ylabel('Voltage (V)');
legend({'Input','Output'});

figure("Name","Sine Wave Fourier Response Spectrum 1Hz")
semilogy(abs(fftshift(fft(vout2(5,:)))))
xlabel('Frequency (Hz)');
ylabel('Amplitude');
xlim([250 750])

freq2 = 1/0.001; % Frequency in 1000Hz
F1 = [0;0;0;0;0;0;0;0];                 % initial values
vout2 = zeros(size(GVec,1),size(t,2)); % Vector to hold outputs
vin2 = zeros(1,size(t,2));          % Vector to hold inputs
vin2(1,1) = sin(2*pi*freq2*t(1));

%Loop through time 
for n=1:size(t,2)-1
    vout2(:,n) = F1;
    vin2(1,n+1) = sin(2*pi*freq2*t(n+1));
    F2 = FVec * vin2(1,n+1);
    F1 = ((CVec./0.001) + GVec)\((CVec./0.001)*F1 + F2);
end
vout2(:,n+1) = F1;

figure("Name","Sine Wave Response 1000Hz")
plot(t, vin2,'r');
hold on
plot(t, vout2(5,:),'b');
title('Sine Wave Response');
xlabel('Time (s)');
ylabel('Voltage (V)');
legend({'Input','Output'});

figure("Name","Sine Wave Fourier Response Spectrum 1000Hz")
semilogy(abs(fftshift(fft(vout2(5,:)))))
xlabel('Frequency (Hz)');
ylabel('Amplitude');
xlim([250 750])

%% Question 4d - Gaussian Pulse Response
F1 = [0;0;0;0;0;0;0;0];                 % initial values
vout3 = zeros(size(GVec,1),size(t,2)); % Vector to hold outputs
vin3 = zeros(1,size(t,2));          % Vector to hold inputs

t = 0:0.001:1;

%https://www.gaussianwaves.com/2014/07/generating-basic-signals-gaussian-pulse-and-power-spectral-density-using-fft/
std = 0.03;
delay = 0.06^2;
vin3(1,1) = 1/(sqrt(2*pi*delay))*(exp(-t(1).^2/(2*delay)));

for n=1:size(t,2)-1
    vout3(:,n) = F1;
    vin3(1,n+1) = 1/(sqrt(2*pi*delay))*(exp(-t(n+1).^2/(2*delay)));
    F2 = FVec * vin3(1,n+1);
    F1 = ((CVec./0.001) + GVec)\((CVec./0.001)*F1 + F2);
end
vout3(:,n+1) = F1;

figure("Name","Gaussian Pulse Response")
plot(t, vin3,'r');
hold on
plot(t, vout3(5,:),'b');
title('Gaussian Pulse Response');
xlabel('Time (s)');
ylabel('Voltage (V)');
legend({'Input','Output'});

figure("Name","Gaussian Pulse Fourier Response Spectrum")
semilogy(abs(fftshift(fft(vout3(5,:)))))
title('Gaussian Pulse Fourier Response');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
xlim([250 750])


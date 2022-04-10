clc; close all; clear all;

%Question 5
%Jinseng Vanderkloot 
%101031534


C=0.25;
Cn = 0.00001;
G1 = 1;
G2 = 1/2;
G3 = 1/8.5;
G4 = 1/0.1;
G0 = 1/1000;
L = 0.2;
alpha = 100;

%F = V1, V2, V3, V4, V5, Vin, I4, IL 

% 0 = Iin - G1(V1-V2) - C*d(V1-V2)/dt             
% 0 = G1(V1-V2) + C*d(V1-V2)/dt - V2*G2 - IL
% 0 = IL - G3*V3 - V3*Cn + In 
% 0 = I4 - G4*(V4-V5)
% 0 = G4*(V4-V5) - Go*Vo 
% Vin = V1
% V4 = I3*Alpha = -Alpha*G3*V3
% dIL = dt(V2-V3)/L  

CVec = [-C, C, 0, 0, 0, 0, 0, 0; 
      C, -C, 0, 0, 0, 0, 0, 0;
      0, 0, Cn, 0, 0, 0, 0, 0;
      0, 0, 0, 0, 0, 0, 0, 0;
      0, 0, 0, 0, 0, 0, 0, 0;
      0, 0, 0, 0, 0, 0, 0, 0;
      0, 0, -alpha*Cn, 0, 0, 0, 0, 0;
      0, 0, 0, 0, 0, 0, 0, -L;];

GVec = [-G1, G1, 0, 0, 0, -1, 0, 0;
        G1, -G1-G2, 0, 0, 0, 0, 0, -1;
        0, 0, -G3, 0, 0, 0, 0, 1;
        0, 0, 0, -G4, G4, 0, 1, 0; 
        0, 0, 0, G4, -G4-G0, 0, 0, 0;
        1, 0, 0, 0, 0, 0, 0, 0; 
        0, 0, -alpha*G3, 1, 0, 0, 0, 0;
        0, 1, -1, 0, 0, 0, 0, 0;];

FVec = [0;0;1;0;0;1;0;0]; %Include new current source for R3

t = 0:0.001:1; 
F1 = [0;0;0;0;0;0;0;0]; 
In = 0.001*randn(1,size(t,2));
vout = zeros(size(GVec,1),size(t,2));
vin= zeros(1,size(t,2));               
vin(1,1) = In(1);      

std = 0.03;
delay = 0.06^2;
vin(1,1) = 1/(sqrt(2*pi*delay))*(exp(-t(1).^2/(2*delay)));


for n=1:size(t,2)-1
    vout(:,n) = F1;
    vin(1,n+1) = 1/(sqrt(2*pi*delay))*(exp(-t(n+1).^2/(2*delay)));
    F2 = FVec * In(n+1); %Add noise to node 3
    F2(6) = vin(n+1); %Add input 
    F1 = ((CVec./0.001) + GVec)\((CVec./0.001)*F1 + F2);
end
vout(:,n+1) = F1;

figure("Name","Current Noise Histogram")
histogram(In);
xlabel('Bin (#)');
ylabel('Current (A)');
title('Noise Histogram');


figure("Name","Gaussian Pulse Response")
plot(t, vin,'r');
hold on
plot(t, vout(5,:),'b');
title('Gaussian Pulse Response');
xlabel('Time (s)');
ylabel('Voltage (V)');
legend({'Input','Output'});

figure("Name","Gaussian Pulse Fourier Response Spectrum")
semilogy(abs(fftshift(fft(vout(5,:)))))
title('Gaussian Pulse Fourier Response');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
xlim([250 750])

%% Changing Cn 

t = 0:0.001:1; 
F1 = [0;0;0;0;0;0;0;0]; 
In = 0.001*randn(1,size(t,2));
vout = zeros(size(GVec,1),size(t,2));
vin= zeros(1,size(t,2));               
vin(1,1) = In(1);      

std = 0.03;
delay = 0.06^2;
vin(1,1) = 1/(sqrt(2*pi*delay))*(exp(-t(1).^2/(2*delay)));

a = linspace(0.000005, 0.00005,4);

for x = 1:size(a,2)
    Cn = a(x);

    CVec = [-C, C, 0, 0, 0, 0, 0, 0; 
      C, -C, 0, 0, 0, 0, 0, 0;
      0, 0, Cn, 0, 0, 0, 0, 0;
      0, 0, 0, 0, 0, 0, 0, 0;
      0, 0, 0, 0, 0, 0, 0, 0;
      0, 0, 0, 0, 0, 0, 0, 0;
      0, 0, -alpha*Cn, 0, 0, 0, 0, 0;
      0, 0, 0, 0, 0, 0, 0, -L;];

    for n=1:size(t,2)-1
        vout(:,n) = F1;
        vin(1,n+1) = 1/(sqrt(2*pi*delay))*(exp(-t(n+1).^2/(2*delay)));
        F2 = FVec * In(n+1); %Add noise to node 3
        F2(6) = vin(n+1); %Add input
        F1 = ((CVec./0.001) + GVec)\((CVec./0.001)*F1 + F2);
    end
    vout(:,n+1) = F1;

    figure("Name","Gaussian Pulse Response ")
    plot(t, vin,'r');
    hold on
    plot(t, vout(5,:),'b');
    title('Gaussian Pulse Response, Capacitance(F) = ',Cn);
    xlabel('Time (s)');
    ylabel('Voltage (V)');
    legend({'Input','Output'});

    figure("Name","Gaussian Pulse Fourier Response Spectrum")
    semilogy(abs(fftshift(fft(vout(5,:)))))
    title('Gaussian Pulse Fourier Response');
    xlabel('Frequency (Hz)');
    ylabel('Amplitude');
    xlim([250 750]);
end

% Increasing the value of the capacitor chnages the effect of the noise on
% the circuit, the value cannot be increase to much or elese the circuit
% fails. 

%% Changing Time Step 

C=0.25;
Cn = 0.00001;
G1 = 1;
G2 = 1/2;
G3 = 1/8.5;
G4 = 1/0.1;
G0 = 1/1000;
L = 0.2;
alpha = 100;

CVec = [-C, C, 0, 0, 0, 0, 0, 0; 
      C, -C, 0, 0, 0, 0, 0, 0;
      0, 0, Cn, 0, 0, 0, 0, 0;
      0, 0, 0, 0, 0, 0, 0, 0;
      0, 0, 0, 0, 0, 0, 0, 0;
      0, 0, 0, 0, 0, 0, 0, 0;
      0, 0, -alpha*Cn, 0, 0, 0, 0, 0;
      0, 0, 0, 0, 0, 0, 0, -L;];

GVec = [-G1, G1, 0, 0, 0, -1, 0, 0;
        G1, -G1-G2, 0, 0, 0, 0, 0, -1;
        0, 0, -G3, 0, 0, 0, 0, 1;
        0, 0, 0, -G4, G4, 0, 1, 0; 
        0, 0, 0, G4, -G4-G0, 0, 0, 0;
        1, 0, 0, 0, 0, 0, 0, 0; 
        0, 0, -alpha*G3, 1, 0, 0, 0, 0;
        0, 1, -1, 0, 0, 0, 0, 0;];

FVec = [0;0;1;0;0;1;0;0]; %Include new current source for R3

t = 0:0.002:1; 
F1 = [0;0;0;0;0;0;0;0]; 
In = 0.002*randn(1,size(t,2));
vout = zeros(size(GVec,1),size(t,2));
vin= zeros(1,size(t,2));               
vin(1,1) = In(1);      

std = 0.03;
delay = 0.06^2;
vin(1,1) = 1/(sqrt(2*pi*delay))*(exp(-t(1).^2/(2*delay)));


for n=1:size(t,2)-1
    vout(:,n) = F1;
    vin(1,n+1) = 1/(sqrt(2*pi*delay))*(exp(-t(n+1).^2/(2*delay)));
    F2 = FVec * In(n+1); %Add noise to node 3
    F2(6) = vin(n+1); %Add input 
    F1 = ((CVec./0.002) + GVec)\((CVec./0.002)*F1 + F2);
end
vout(:,n+1) = F1;

figure("Name","Gaussian Pulse Response - Time Step = 0.002")
plot(t, vin,'r');
hold on
plot(t, vout(5,:),'b');
title('Gaussian Pulse Response- Time Step = 0.002');
xlabel('Time (s)');
ylabel('Voltage (V)');
legend({'Input','Output'});

figure("Name","Gaussian Pulse Fourier Response Spectrum")
semilogy(abs(fftshift(fft(vout(5,:)))))
title('Gaussian Pulse Fourier Response- Time Step = 0.002');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
xlim([0 500]);

t = 0:0.0005:1; 
F1 = [0;0;0;0;0;0;0;0]; 
In = 0.0005*randn(1,size(t,2));
vout = zeros(size(GVec,1),size(t,2));
vin= zeros(1,size(t,2));               
vin(1,1) = In(1);      

std = 0.03;
delay = 0.06^2;
vin(1,1) = 1/(sqrt(2*pi*delay))*(exp(-t(1).^2/(2*delay)));


for n=1:size(t,2)-1
    vout(:,n) = F1;
    vin(1,n+1) = 1/(sqrt(2*pi*delay))*(exp(-t(n+1).^2/(2*delay)));
    F2 = FVec * In(n+1); %Add noise to node 3
    F2(6) = vin(n+1); %Add input 
    F1 = ((CVec./0.0005) + GVec)\((CVec./0.0005)*F1 + F2);
end
vout(:,n+1) = F1;

figure("Name","Gaussian Pulse Response - Time Step = 0.0005")
plot(t, vin,'r');
hold on
plot(t, vout(5,:),'b');
title('Gaussian Pulse Response- Time Step = 0.0005');
xlabel('Time (s)');
ylabel('Voltage (V)');
legend({'Input','Output'});

figure("Name","Gaussian Pulse Fourier Response Spectrum")
semilogy(abs(fftshift(fft(vout(5,:)))))
title('Gaussian Pulse Fourier Response- Time Step = 0.0005');
xlabel('Frequency (Hz)');
ylabel('Amplitude');
xlim([500 1500]);

clc; close all; clear all;

%Jinseng Vanderkloot 
% 101031534

%This is PA7 with DC and AC sweep of circuit for A4 

% 1A Differential Equations

% 0 = Iin - G1(V1-V2) - C*d(V1-V2)/dt             
% 0 = G1(V1-V2) + C*d(V1-V2)/dt - V2*G2 - IL
% 0 = IL - G3*V3
% 0 = I4 - G4*(V4-V5)
% 0 = G4*(V4-V5) - Go*Vo 
% Vin = V1
% V4 = I3*Alpha = -Alpha*G3*V3
% dIL = dt(V2-V3)/L  


% 1B Frequency Domain 
        %Values to put into matrix 

% 0 = Iin - G1(V2-V1) - jwC*(V2-V1)
        % V1*(-G1)   Iin   V2*(G1)
% 0 = G1(V2-V1) + jwC*(V2-V1) - V2*G2 - IL
        % V1(G1)   V2(-G1 + -G2)   -IL 
% 0 = IL - G3*V3
        % IL    V3(-G3)
% 0 = I4 - G4*(V4-V5)
        % I4   V4(-G4)  V5 (G4)  
% 0 = G4*(V4-V5) - Go*V5
        %V4(G4)    V5 (-G4-G0) 
% Vin = V1 
        %Vin(t)     V1 
% V4 = I3*Alpha = -Alpha*G3*V3
        %V4    V3 (-Alpha*G3) 
% IL = (1/jwL)*(V2-V3)
        %1/L   V2    -V3 

C=0.25;
G1 = 1;
G2 = 1/2;
G3 = 1/8.5;
G4 = 1/0.1;
G0 = 1/1000;
L = 0.2;
alpha = 100;

%F = V1, V2, V3, V4, V5, Vin, I4, IL 

CVec = [-C, C, 0, 0, 0, 0, 0, 0; 
      C, -C, 0, 0, 0, 0, 0, 0;
      0, 0, 0, 0, 0, 0, 0, 0;
      0, 0, 0, 0, 0, 0, 0, 0;
      0, 0, 0, 0, 0, 0, 0, 0;
      0, 0, 0, 0, 0, 0, 0, 0;
      0, 0, 0, 0, 0, 0, 0, 0;
      0, 0, 0, 0, 0, 0, 0, -L;]; %Not 1/L 

GVec = [-G1, G1, 0, 0, 0, -1, 0, 0;
        G1, -G1-G2, 0, 0, 0, 0, 0, -1;
        0, 0, -G3, 0, 0, 0, 0, 1;
        0, 0, 0, -G4, G4, 0, 1, 0; 
        0, 0, 0, G4, -G4-G0, 0, 0, 0;
        1, 0, 0, 0, 0, 0, 0, 0; 
        0, 0, -alpha*G3, 1, 0, 0, 0, 0;
        0, 1, -1, 0, 0, 0, 0, 0;];

FVec = [0;0;0;0;0;1;0;0];

%DC Sweep
vin = linspace(-10,10);
threeV = zeros(size(vin));
outputV = zeros(size(vin));
for a = 1:size(vin,2)
    V = (GVec+CVec)\(FVec.*vin(a));
    threeV(a) = V(3);
    outputV(a) = V(5);
end 

figure('name',"DC Sweep")
plot(vin,threeV,'g');
hold on
plot(vin,outputV,'r');
hold off
title('DC Sweep');
xlabel('Input Voltage (V)');
ylabel('Output Voltage (V)');
legend('V3','Vout');

% AC Sweep
freq = linspace(0.01, 100);
threeV = zeros(size(vin));
outputV = zeros(size(freq));
for a = 1:size(freq,2)
    w =(2*pi*freq(a));
    jw = 1i*w; % 1i makes value complex 1 -> 0+1j 
    V = (GVec+jw*CVec)\FVec; 
    outputV(a) = V(5);
    threeV(a) = V(3);
end

figure("Name","AC Sweep")
plot(2*pi*freq,real(outputV),'g'); %Real() gets rid of imag value warning 
xlim([0 100])
hold on
plot(2*pi*freq,real(threeV),'r');
hold off
xlabel('Frequency  (Hz)')
ylabel('Output (V)')
legend('Vout','V3');


%Capacitor normal Distribution 
C = 0.25 + 0.05 * randn(1000,1); %Random C values 
outputV = zeros(size(C));
for a = 1:size(C,1)
    CVec(1:2,1:2)=C(a); %CVec does not auto update for changing C so need this statement
    V = (GVec+(1i*pi)*CVec)\FVec; 
    outputV(a) = V(5);
end 

figure("Name","Capacitor Distribution")
histogram(C);
xlabel('Capacitance (F)');
ylabel('Number (#)');

figure("Name","Output Distribution");
histogram(real(outputV));
xlabel('Vout (V)');
ylabel('Number (#)');


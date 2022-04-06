clc
clear all
close all 

%Modify Assignment 3 Questions 2 to calculate the resistance for variable
%voltage across area 


%Initial
m0 = 9.10938215e-31;            % electron mass
mn = 0.26*m0;                   % Effective mass
Temp = 300;                     % Inital Temp (K)
kb = 1.3806504e-23;             % Boltzmann constant
tmn = 0.2e-12;                  % Mean time between collision 
q = 1.602e-19;                  % Charge of an Electron (Assignment 3 for force) 

% Region Area 
wArea = 200e-9;
lArea = 100e-9;

%Thermal Velocity (Question 1.A) 
vt=sqrt((2*kb*Temp)/mn);        % Sim in 2D so (2*kb*Temp), 3D is (3*kb*Temp)

%Electron motion
numElec = 1000;                  %Number of simulated Electrons.
numEPlot = 30;                  %Number of plotted Electrons
dt = (lArea*wArea)/2;           %Typically 1/100 of region size
stepsTot = 200;                 %Total amount of steps (1000 was a long simulation)
tTot= stepsTot*dt;              %Total Simulation time
x = zeros(1,numElec);           %Inital X matrix
y = zeros(1,numElec);           %Inital y matrix
vx = zeros(1,numElec);          %Inital velocity x matrix
vy = zeros(1,numElec);          %Inital velocity y matrix
vtot = zeros(1,numElec);        %Inital velocity matrix
xfield = zeros(1,numElec);      %Inital Position in E field x
yfield = zeros(1,numElec);      %Inital Position in E field y
vxForce = zeros(1,numElec);      %Inital Force due to position in E field x
vyForce = zeros(1,numElec);      %Inital Force due to position in E field y

%Probability of Scatter
scatOn = 1;                     %Turn Scatter on (1) or off(0)
Pscatter = 1-exp(-dt/tmn);      %Scatter Equation
tScatter = zeros(1,numElec);

%Bottle Neck Boundary reduced to prevent bleeding
% X limit (no particle between X1 and X2 - Cosider Y limits)
boxX1=80e-9;
boxX2=120e-9;
% Y Limit (no particles between 0 and Y1 and Y2 to Y limit)
boxY1 = 40e-9;
boxY2 = 60e-9;

%Electric Field
x0 = 0.5; %voltage at right side of area increased to 0.5V to see huge effect
x1 = 0; %Voltage at left side of area
nx = wArea*1e9; % # of colums (changes to have same ratio as A1 area 2/1)
ny = lArea*1e9; % # of rows
xBox = 40; %Width of box (changes to same box ratio as A1 area l= 40 and w=40)
yBox = 40; %Hight of box
boxCond = 0.01;
% V=A2_Function(nx, ny, xBox, yBox, boxCond, x0, x1);
% Vmap = reshape(V, [ny, nx]);
% [Ex,Ey] = gradient(-Vmap);

% figure('name', 'Voltage Surface Plot')
% surf(Vmap'),view(120,20);
% 
% figure('name', 'Electric Field Vector Plot');
% quiver(Ex,Ey);


%Electron Graph
for cnt = 1:numElec
    x(cnt)=rand()*wArea;
    y(cnt)=rand()*lArea;
    %If the electrons are place in the box, re-roll position
    while (x(cnt)>=boxX1 && x(cnt)<=boxX2 && (y(cnt)<=boxY1 || y(cnt>=boxY2))) %Relocate them if in boundary
        x(cnt)=rand()*wArea;
        y(cnt)=rand()*lArea;
    end
    vx(cnt)=(vt/sqrt(2))*randn();  % velocity * Gaussian dist
    vy(cnt)=(vt/sqrt(2))*randn();  % velocity * Gaussian dist
    vtot(cnt)= sqrt (vx(cnt)^2)+(vy(cnt)^2);
    colors= rand(numElec,3);
end

%Boundary Energy/Velocity loss coefficient (reduction in velocity =
%reduction in temprature)
vloss = 0.95;

vol = 0.1 : 0.5 : 10.1;
vStep = 1;

%% S3 Main Loop
for voltage = vol

    V=A2_Function(nx, ny, xBox, yBox, boxCond, voltage, x1);
    Vmap = reshape(V, [ny, nx]);
    [Ex,Ey] = gradient(-Vmap);
    t=0;
    intCNT = 2;
    while t < tTot
        t = t + dt;

        %Store old position
        oldx = x;
        oldy = y;

        %Update to new position
        x(1:numElec) = x(1:numElec) + (vx(1:numElec).*dt);
        y(1:numElec) = y(1:numElec) + (vy(1:numElec).*dt);

        vtot(1:numElec)= sqrt ((vx(1:numElec).^2)+(vy(1:numElec).^2));

        for check = 1:numElec
            %Scatter
            if scatOn==1
                if Pscatter > rand()
                    vx(check)=(vt/sqrt(2))*randn();
                    vy(check)=(vt/sqrt(2))*randn();
                    tScatter(check)= 0; %If collision, time goes to 0
                else
                    tScatter(check)= tScatter(check) + dt; %track time increaing while no collision
                end
            end

            % Bring X and Y of particle location to values able to be tracked
            % in electric field mesh to apply field forces
            xfield(check) = round(x(check)*1e9);
            yfield(check) = round(y(check)*1e9);

            %Make sure field effects still apply to particles which go across
            %area
            if xfield(check)>200
                xfield(check) = 200;
            end
            if xfield(check)<1
                xfield(check) = 1;
            end
            if yfield(check)>100
                yfield(check) = 100;
            end
            if yfield(check)<1
                yfield(check) = 1;
            end
            %Apply field force to particles
            vxForce(check)= ((q *Ex(yfield(check),xfield(check))/(wArea/nx))/mn)*dt; % V=(q*E)/mass F=q*E
            vyForce(check)= ((q *Ey(yfield(check),xfield(check))/(lArea/ny))/mn)*dt; % V=(q*E)/mass F=q*E
            vx(check) = vx(check)+ vxForce(check);
            vy(check) = vy(check)+ vyForce(check);

            %Apply Boundary Conditions
            %If bottom contact, bounce off in opposite direction
            if (y(check)<=0)
                y(check) = 0;
                vy(check) = -vy(check);
            end
            %If top contact, bounce off in opposite directio
            if (y(check)>=lArea)
                y(check) = lArea;
                vy(check) = -vy(check);
            end
            %If left side of box, come out right side
            if(x(check)<=0)
                x(check) = x(check) + wArea;
            end
            %If right side of box, come out left side
            if(x(check)>=wArea)
                x(check) = x(check) - wArea;
            end

            %Apply bottle neck conditions
            %If contact on left walls of boundary (not in Gap)
            if (oldx(check)<boxX1 && x(check)>=boxX1 && (y(check)<= boxY1 || y(check)>= boxY2) && oldx(check)>((1/5)*wArea) && oldx(check)<((4/5)*wArea))
                x(check)=boxX1;
                vx(check) = -(vx(check)*vloss);
            end
            %If contact on right walls of boundary (not in Gap)
            if (oldx(check)>boxX2 && x(check)<=boxX2 && (y(check)<= boxY1 || y(check)>= boxY2) && oldx(check)>((1/5)*wArea)&& oldx(check)<((4/5)*wArea))
                x(check)=boxX2;
                vx(check) = -(vx(check)*vloss);
            end
            %If contact with bottom boundary in gap
            if (x(check)>boxX1 && x(check)< boxX2 && oldy(check)>boxY1 && y(check)<= boxY1)
                y(check)= boxY1;
                vy(check) = -(vy(check)*vloss);
            end
            %If contact with top boundary in gap
            if (x(check)>boxX1 && x(check)< boxX2 && oldy(check)<boxY2 && y(check)>=boxY2)
                y(check)=boxY2;
                vy(check) = -(vy(check)*vloss);
            end
        end
        %Current over time with Field - Mistake multiplying current by L and W
        %- Only multiply with W cause density is A/m * m
        Ix(:,intCNT)= q*10e19*mean(vx)* wArea; %Current = e * n * vd * Area (Particle density)

        intCNT = intCNT +1;
    end
    AvgIx (:,vStep) = mean (Ix);
    vStep = vStep +1;
end
figure(1); hold on;
p = polyfit(vol,AvgIx,1);
line = polyval(p,vol);
plot(vol,AvgIx,"r");
plot(vol,line,'b');
title('Current in Area for Changing voltage'),xlabel('Voltage(V)', 'FontSize', 10), ylabel('Current(A)', 'FontSize', 10);
legend ('Current', 'Linear Fit');
fprintf('Linear equation is I = %e * V + %e \r', p(1), p(2));
R3 = 1/(p(1));
fprintf('The resistance for R3 = %e Ohms \r',R3);
hold on;
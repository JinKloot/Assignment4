%Make this into a function like in the intro to the lab which inputs size
% of area, size of boxes (placed into the middle x and bottom and top).
% and box conduction values

%Jinseng Vanderkloot 
%101031534

function [V] = A2_Function(nx, ny, xBox, yBox,boxCond,x0,x1)
%Inputs: 
%Area x dimension, Area y dimension, box x dimension in middle of area, 
%Box y dimension from bottom to high and from top down, box conductivity,
%x0 = volatge at left side, x1 = volatge at right side. 

global Carea %NEEDS TO BE GLOBAL 

% Add bottleneck
Carea = ones(nx,ny); %set conduction area to 1
% In area, place boxes with new conduction (faster than for loop) 
Carea(nx/2 - xBox/2:nx/2 + yBox/2,1:yBox) = boxCond; %Bottom Box
Carea(nx/2 - xBox/2:nx/2 + yBox/2,ny-yBox:ny) = boxCond; %Top Box

G = sparse(nx*ny,ny*nx);
F = zeros(nx*ny,1);

for i = 1:nx
    for j = 1:ny
        n = j + (i-1) * ny;     % middle
        nxm = j + (i-2) * ny;   % right
        nxp = j + i * ny;       % left
        nym = j-1 + (i-1) * ny; % top
        nyp = j+1 + (i-1) * ny; % down
        if i == 1 %Left Boundary V=Vo
            G(n,n) = 1;
            F(n,1) = x0;
        elseif i == nx %Right Boundary V=Vo
            G(n,n) = 1;
            F(n,1) = x1;
        elseif j == 1 %Bottom Boundary (Free)
            bxm = (Carea(i,j) + Carea(i-1,j)) / 2;
            bxp = (Carea(i,j) + Carea(i+1,j)) / 2;
            byp = (Carea(i,j) + Carea(i,j+1)) / 2;

            G(n,n) = -(bxm + bxp + byp);
            G(n,nxm) = bxm;
            G(n,nxp) = bxp;
            G(n,nyp) = byp;
        elseif j == ny %Top Boundary (Free)
            bxm = (Carea(i,j) + Carea(i-1,j)) / 2;
            bxp = (Carea(i,j) + Carea(i+1,j)) / 2;
            bym = (Carea(i,j) + Carea(i,j-1)) / 2;

            G(n,n) = -(bxm + bxp + bym);
            G(n,nxm) = bxm;
            G(n,nxp) = bxp;
            G(n,nym) = bym;
        else %Middle
            bxm = (Carea(i,j) + Carea(i-1,j)) / 2;
            bxp = (Carea(i,j) + Carea(i+1,j)) / 2;
            byp = (Carea(i,j) + Carea(i,j+1)) / 2;
            bym = (Carea(i,j) + Carea(i,j-1)) / 2;

            G(n,n) = -(bxm + bxp + bym + byp);
            G(n,nxm) = bxm;
            G(n,nxp) = bxp;
            G(n,nym) = bym;
            G(n,nyp) = byp;
        end
    end
end

V = G\F;

end
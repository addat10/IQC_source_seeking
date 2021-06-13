% Plot the shape of the OS interaction field
clear
clc
addpath('./flocking_scripts/.')
%% Parameters of the flocking field
epsilon         = 0.1; % Used to define the sigma_norm
da              = 7;   % Equilibrium distance(or desired distance) to neighbours
ra              = 1.2 * da; % Sensing radius (no interaction with agents outside ra disc)
h               = 0.9; % Bump function goes to zero after h
%% Plot the scalar field
x=-5:0.1:5;
y=-5:0.1:5;
V=zeros(size(x,2),size(y,2));
for i=1:size(x,2)
    for j=1:size(y,2)
        q=[x(1,i);y(1,j)];
        V(i,j)= psi_alpha(sigma_norm(q,epsilon),ra,da,h); 
    end
end

%% Plot
[X,Y]=meshgrid(x,y);
figure()
surf(X,Y,V)
figure()
contour(X,Y,V,10)

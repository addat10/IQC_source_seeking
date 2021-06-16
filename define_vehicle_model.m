function [G]=define_vehicle_model()
% Generic second order vehicle model
A=[0 1; 0 -1]; 
B=[0;-1];
C=[1 0];
D=0;
G=ss(A,B,C,D);
end
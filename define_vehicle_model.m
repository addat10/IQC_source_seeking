function [G]=define_vehicle_model()
% This function defines a closed-loop vehicle model that has flocking force
% as input and position and velocity as output
% General LTI vehicle models should be in the form 1/s*G_vel where G_vel is
% a velocity tracking controller.

%% Generic second order vehicle: mass with friction in 1D
c_damp=1;
A=[0 1; 0 -c_damp]; 
B=[0;1];
C=[1 0];
D=0;
G=ss(A,B,C,D);
end
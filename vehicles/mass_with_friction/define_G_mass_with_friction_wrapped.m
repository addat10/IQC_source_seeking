function [G]=define_G_mass_with_friction_wrapped(dim,c_damp,mass,step_size)
% This function defines a closed-loop vehicle model that has flocking force
% as input and position and velocity as output
% General LTI vehicle models should be in the form 1/s*G_vel where G_vel is
% a velocity tracking controller.

%% Generic second order vehicle: mass with friction in 1D
inv_mass=1/mass;
A=[0 1; 0 -c_damp*inv_mass]; 
B=[0;1*inv_mass*step_size];
C=[1 0];
D=0;
% Higher dimensional vehicles
A_hat=kron(A,eye(dim));
B_hat=kron(B,eye(dim));
C_hat=kron(C,eye(dim));
D_hat=kron(D,eye(dim));

G=ss(A_hat,B_hat,C_hat,D_hat);
end
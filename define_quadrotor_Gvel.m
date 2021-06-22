function [G]=define_quadrotor_Gvel(dim)
% This function defines a closed-loop quadrotor model that has flocking
% force as input and it preserves the integral action

%% Define open loop linearized quadrotor model
% m=0.640; % Mass of the quadrotor
% g=9.81; % gravitational acceleration
% A = [ 0 1 0 0 0 0 0 0  0 0 0 0 ;  % x
%       0 0 0 0 0 0 0 0 -g 0 0 0 ;  % xdot
%       0 0 0 1 0 0 0 0  0 0 0 0 ;  % y
%       0 0 0 0 0 0 0 0  0 0 g 0 ;  % ydot
%       0 0 0 0 0 1 0 0  0 0 0 0 ;  % z
%       0 0 0 0 0 0 0 0  0 0 0 0 ;  % zdot
%       0 0 0 0 0 0 0 1  0 0 0 0 ;  % psi
%       0 0 0 0 0 0 0 0  0 0 0 0 ;  % psidot
%       0 0 0 0 0 0 0 0  0 1 0 0 ;  % theta
%       0 0 0 0 0 0 0 0  0 0 0 0 ;  % thetadot
%       0 0 0 0 0 0 0 0  0 0 0 1 ;  % phi
%       0 0 0 0 0 0 0 0  0 0 0 0 ]; % phidot
% B = [ 0   0 0 0 ;
%       0   0 0 0 ;
%       0   0 0 0 ;
%       0   0 0 0 ;
%       0   0 0 0 ;
%       1/m 0 0 0 ;
%       0   0 0 0 ;
%       0   1 0 0 ;
%       0   0 0 0 ;
%       0   0 1 0 ;
%       0   0 0 0 ;
%       0   0 0 1 ];

% Define double integrator to generate references
if dim~=2
    error('Currently just works with dimension=2')
end
c_damp=1;
A=kron([0 1; 0 -c_damp],eye(dim)); 
B=kron([0;1],eye(dim));
C=eye(2*dim);
D=zeros(2*dim,dim);
G_double_int=ss(A,B,C,D);
load('G_quad_cl');

G=G_quad_cl*G_double_int;
%G=ss(A_hat,B_hat,C_hat,D_hat);
end
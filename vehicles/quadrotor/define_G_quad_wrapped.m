function [G]=define_G_quad_wrapped(dim)
% This function defines a closed-loop quadrotor model that has flocking
% force as input and it preserves the integral action

% Define double integrator to generate references
if dim~=2
    error('Currently just works with dimension=2')
end
c_damp=5;
inv_mass=(1/0.640);
step_size=0.2;
A=kron([0 1; 0 -c_damp*inv_mass],eye(dim)); 
B=kron([0;1*inv_mass*step_size],eye(dim));
C=eye(2*dim);
D=zeros(2*dim,dim);
G_double_int=ss(A,B,C,D);

G_quad_cl=get_quad_G_cl();
G=G_quad_cl*G_double_int;
% Compare the wrapped plant with the generic vehicle model
sigma(G_double_int)
hold on
sigma(G,'r')
end
function [Psi_GI,M]=define_ZF_multiplier(m,L,G_veh)
% This function returns the augmented model matrices for Psi*[G;I] which
% are required for analysis with LMIs. In case of static multipliers, Psi
% is just an identity map of aproriate dimension. 
% Note that this function returns a single fixed multiplier and not a class
% of multipliers.

[A,B,C,D]=ssdata(G_veh);

% Negative feedback(if gradient is the delta block, we want grad descent)
B=-B;

% Define Psi
[ny,nu]=size(D);
psi=ss([],[],[],eye(nu+ny));

% Build Psi*[G;I] with Dynamic Multiplier
Psi_GI=psi*ss(A,B,[C;0,0],[D;1]);


% Static multiplier
M=[-m*L,    (m+L)/2;...
   (m+L)/2, -1];
end
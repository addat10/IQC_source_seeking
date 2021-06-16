function [Psi_GI,M]=define_ZF_multiplier(m,L,G)
[A,B,C,D]=ssdata(G);
% Build G with Dynamic Multiplier
Psi_GI=ss(A,B,[C;0,0],[D;1]);

% Static multiplier
M=[-m*L,    (m+L)/2;...
   (m+L)/2, -1];
end
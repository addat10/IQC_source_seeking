%% Design of a Formation Controller
% -------------------------------------------------------------------------
% script   : exercise_MAS_qcopter
% -------------------------------------------------------------------------
% Author   : Marcus Bartels
% Version  : December 10th, 2015
% Copyright: MB, 2015
% -------------------------------------------------------------------------
%
% 1. Define the agent model given in the appendix of the Lecture Notes
% 2. Formulate agent + uncertain topology in LFT form
% 3. Generate a generalized plant
% 4. Synthesize a robust formation controller
%
% IMPORTANT: This is the exercise version. 
%   Lines marked with <FIXME> are to be filled by you!
%   Unless you do so, the script will not be executable.
%
% -------------------------------------------------------------------------

%% state space model of the quadrocopter
g = 9.81;   % Gravity constant

m = 0.640;  %  Mass of the Quadrocopter
  
A = [ 0  1  0  0  0  0  0  0  0  0  0  0   ;
      0  0  0  0  0  0  0  0 -g  0  0  0   ;
      0  0  0  1  0  0  0  0  0  0  0  0   ;
      0  0  0  0  0  0  0  0  0  0  g  0   ;
      0  0  0  0  0  1  0  0  0  0  0  0   ;
      0  0  0  0  0  0  0  0  0  0  0  0   ;
      0  0  0  0  0  0  0  1  0  0  0  0   ;
      0  0  0  0  0  0  0  0  0  0  0  0   ;
      0  0  0  0  0  0  0  0  0  1  0  0   ;
      0  0  0  0  0  0  0  0  0  0  0  0   ;
      0  0  0  0  0  0  0  0  0  0  0  1   ;
      0  0  0  0  0  0  0  0  0  0  0  0 ] ;
   
B = [ 0  0  0  0  0 1/m 0  0  0  0  0  0   ;
      0  0  0  0  0  0  0  1  0  0  0  0   ;
      0  0  0  0  0  0  0  0  0  1  0  0   ;
      0  0  0  0  0  0  0  0  0  0  0  1 ]';

nx = size(A,1); % number of states
  
C = zeros(3,nx);
C(1,1) = 1; % x
C(2,3) = 1; % y
C(3,5) = 1; % z

nu = size(B,2); % number of inputs
ny = size(C,1); % number of outputs

D = zeros(ny,nu);

P = ss(A,B,C,D);
Paug = ss(A,B,[C; eye(nx)], [D; zeros(nx,nu)]);
% Augmented plant with state vector as additional output 

%% setup of the generalized plant
% Factorize C
[U,S,V]=svd(C);
Dd = U*sqrt(S(:,1:size(U,2)));
Cd = sqrt(S)*V';
nd = size(Cd,1); % size of Delta

% Generate an LFT model of one agent
Blft = [zeros(nx,nd), B];
Clft = [Cd; C; eye(nx)];
Dlft = [zeros(nd,nd),   zeros(nd,nu);
            Dd,         zeros(ny,nu);
        zeros(nx,nd),   zeros(nx,nu)];
Plft = ss(A,Blft,Clft,Dlft);    

% Define the shaping filters
ws=10;
Ms=0.001;

wk=1000;
Mk=10;
ck=1000;

s=tf('s');
w_S = (ws/Ms)*1/(s+ws);
Ws = diag(ones(1,ny))*w_S;

w_KS = (ck/Mk)*(s+wk)/(s+ck*wk);
Wk = diag(ones(1,nu))*w_KS;

% Compose the generalized plant
systemnames = 'Plft Ws Wk';
inputvar    = sprintf('[wd(%d); r(%d); u(%d)]', nd,ny,nu);
input_to_Plft = '[wd; u]';
input_to_Ws = sprintf('[r-Plft(%d:%d)]', nd+1,nd+ny);
input_to_Wk = '[u]';
outputvar   = sprintf('[Plft(1:%d); Ws; Wk; r-Plft(%d:%d); -Plft(%d:%d)]', nd,nd+1,nd+ny,nd+ny+1,nd+ny+nx);
cleanupsysic = 'yes';
GP = sysic; % compose the defined setup

NMEAS = ny+nx;  % number of measured outputs
NCON = nu;      % number of control inputs

%% Controller Synthesis
[K,CL,gam1,INFO] = hinfsyn(GP,NMEAS,NCON,'METHOD','LMI');

%% Simulation Setup
% Define the topology
L = [1,0,-1,0;...
     0,1,-1,0;...
     -0.5,-0.5,1,0;...
     -0.5,0,-0.5,1];
Lq = kron(L,eye(ny));

% Define the reference
ref = [2 4 6 8]';
iny = [1; 0; 0];
r = kron(eye(4),iny)*ref;

%%
G_cl=feedback(Paug*K,eye(15));
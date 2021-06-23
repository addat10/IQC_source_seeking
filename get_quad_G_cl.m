function G_quad_wrapped=get_quad_G_cl()
% This function defines the quadrotor model, designs a tracking controller
% and returns wrapped up closed loop of the quadrotor that takes as input
% position and velocities and outputs quadrotor positions.
% -------------------------------------------------------------------------
    %% Define the model of the quadrocopter
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
    
    % Note the transpose at the end of B
    B = [ 0  0  0  0  0 1/m 0  0  0  0  0  0   ;
          0  0  0  0  0  0  0  1  0  0  0  0   ;
          0  0  0  0  0  0  0  0  0  1  0  0   ;
          0  0  0  0  0  0  0  0  0  0  0  1 ]';
    C = eye(12);  
    D = zeros(12,4);
    P = ss(A,B,C,D);
    ref_dim=6; % Position and velocity references
    %% Design a tracking controller
    % Hinf design
    %G_quad_wrapped = hinf_design_2DOF_four_block(P,ref_dim);
    
    G_quad_wrapped = lqr_design_with_static_FF(P,ref_dim);
      
end
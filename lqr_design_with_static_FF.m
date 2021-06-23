function G_quad_wrapped = lqr_design_with_static_FF(P,ref_dim)
% This function designs an LQR controller

    [A,B,C,D]=ssdata(P);   
    %LQR design    
    [nx,nu]=size(B);    
    Q=eye(nx);
    R=0.01*eye(nu);
    N=zeros(nx,nu);    
    [K,S,lambdas] = lqr(P,Q,R);
    A_cl=(A-B*K);
    K_ff=(-C*inv(A_cl)*B)\eye(nx);    
    G_cl=ss(A_cl,B*K_ff,C,0);
    %% Wrap the closed loop such that it takes as input 
    % [x_des,y_des,xdot_des,ydot_des]
    % z_des and zdot_des are fixed to zero at the moment
    B_1=[ 1, 0, 0, 0;
          0, 0, 1, 0;
          0, 1, 0, 0;
          0, 0, 0, 1;
          0, 0, 0, 0;
          0, 0, 0, 0];
    B_in=[B_1;zeros(6,4)];
    C_out=zeros(2,12);
    C_out(1,1)=1;C_out(2,3)=1;
    G_quad_wrapped=C_out*G_cl*B_in;  
end
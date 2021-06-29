design_controller_for_quadrotor;
S1=[1 0 0 0;
    0 1 0 0;
    0 0 0 0;
    1 0 0 0;
    0 0 1 0;
    0 1 0 0;
    0 0 0 1;
    0 0 0 0;
    0 0 0 0];
S=[S1;zeros(6,4)];

B_in=S;
C_out=[eye(2),zeros(2,13)];
G_quad_cl=C_out*G_cl*B_in;

save('G_quad_cl','G_quad_cl')
function [status,X]=verify_exp_stab_ZF_FBCC(G_veh,alpha,sec_1,sec_2,cond_tol,tol)
% This function runs the analysis LMI with cvx and returns the status and
% the storage function matrix P. 

    [A_G,B_G,C_G,D_G]=ssdata(G_veh);
    
    % Negative feedback(if gradient is the delta block, we want grad descent)
    B_G=-B_G;
    
    nx=size(A_G,1);
    [ny,nu]=size(D_G);    
    if nu~=ny
        error('Delta block must be square')
    end
    n_delta=nu;   
    
    % Define the filter H used in ZF multiplier
    % Impulse response must be non-negative at all times
    % Must be stable
    % Impulse response must have an L1 norm of less than 1
    a=-5; % Pole of the filter H
    A_h=a;
    B_h=1;
    C_h=-a; % Makes sure that the L1 norm of the impulse response is 1
    
    % Define the augmented ny must be equal to nu
    A_PSI=(A_h-2*alpha)*eye(n_delta); % For the alpha IQC condition
    B_PSI=B_h*[sec_2*eye(n_delta),-1*eye(n_delta)];
    C_PSI=[-C_h*eye(n_delta);zeros(n_delta);zeros(2*n_delta,n_delta)];
    D_PSI=[ sec_2*eye(n_delta),-1*eye(n_delta);...
            -sec_1*eye(n_delta), 1*eye(n_delta);...
            eye(2*n_delta)];
    
    PSI=ss(A_PSI,B_PSI,C_PSI,D_PSI);      
    
    % Build Psi*[G;I] with Dynamic Multiplier
    PSI_GI=PSI*ss(A_G,B_G,[C_G;zeros(n_delta,nx)],[D_G;eye(n_delta)]);
    
    
    [status,X]=verify_exp_stab_ZF_FBCC_LMI(PSI_GI,alpha,sec_1,sec_2,cond_tol,tol);
end
function [status,X]=verify_exp_stab_ZF_FBCC_LMI(Psi_GI,alpha,sec_1,sec_2,cond_tol,tol)
% This function runs the exp-stab analysis KYP Lemma LMI
    status=false;
    [A,B,C,D]=ssdata(Psi_GI);
    n=size(A,1);
    dim=size(B,2);
    
    M=[ zeros(dim),eye(dim);...
        eye(dim),zeros(dim)];
    
    cvx_begin sdp 
    cvx_precision high    
    
    variable X(n,n) symmetric
    variable lambda_1
    
    variable Q(dim,dim) symmetric
    variable S(dim,dim)
    variable R(dim,dim) symmetric

    % Multiplier class
    % Create the multiplier class
    d=[ sec_1,sec_1;
        sec_1,sec_2;
        sec_2,sec_1;
        sec_2,sec_2];
    E=@(delta)[eye(dim);diag(delta)];
    
    C1=E(d(1,:))'*[Q,S;S',R]*E(d(1,:));
    C2=E(d(2,:))'*[Q,S;S',R]*E(d(2,:));
    C3=E(d(3,:))'*[Q,S;S',R]*E(d(3,:));
    C4=E(d(4,:))'*[Q,S;S',R]*E(d(4,:));

    % FDI via KYP Lemma
    L1=[A'*X + X*A + 2*alpha*X,    X*B;
        B'*X,                      zeros(dim)];
    L2=[C';D']*[lambda_1*M,zeros(2*dim);zeros(2*dim),[Q,S;S',R]]*[C,D];   
    LMI=L1+L2;

    minimize 1; 
    subject to:
    R<= 0;
    C1>=0;
    C2>=0;
    C3>=0;
    C4>=0;
    lambda_1>=0;
    X >= tol*eye(n);
    X <= cond_tol*tol*eye(n);
    
    LMI<=0;    
    cvx_end 
    if strcmp('Solved',cvx_status)
        status=true;
    end
end
function [status,P]=verify_exp_stab_ZF_NC_delete(G_veh,alpha,sec_1,sec_2,cond_tol,tol)
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
    a=-2; % Pole of the filter H
    A_h=a;
    B_h=1;
    C_h=-a; % Makes sure that the L1 norm of the impulse response is 1
    
    
    % Define the augmented ny must be equal to nu
    A_PSI=(A_h)*eye(n_delta); % For the alpha IQC condition
    B_PSI=B_h*[-sec_1*eye(n_delta),1*eye(n_delta)];
    C_PSI=[-C_h*eye(n_delta);zeros(n_delta);zeros(2*n_delta,n_delta)];
    D_PSI=[ -sec_1*eye(n_delta),1*eye(n_delta);...
            sec_2*eye(n_delta), -1*eye(n_delta);...
            sec_2*eye(n_delta),-1*eye(n_delta);...
            -sec_1*eye(n_delta), 1*eye(n_delta)];
    
    PSI=ss(A_PSI,B_PSI,C_PSI,D_PSI);
    M=[ zeros(n_delta),eye(n_delta);...
        eye(n_delta),zeros(n_delta)];  
    
    % Build Psi*[G;I] with Dynamic Multiplier
    PSI_GI=PSI*ss(A_G,B_G,[C_G;zeros(n_delta,nx)],[D_G;eye(n_delta)]);
    
    % Run the exponenetial analysis LMI
    [status,P]=verify_exp_stab_alpha_ZF_LMI(PSI_GI,M,alpha,cond_tol,tol);
    
end
function [status,P]=verify_exp_stab_alpha_ZF_LMI(PSI_GI,M,alpha,cond_tol,tol)
 % This function runs the exp-stab analysis KYP Lemma LMI
    status=false;
    % Get state-space matrics
    [A,B,C,D]=ssdata(PSI_GI);    
    [n,nu]=size(B);

    cvx_begin sdp 
    cvx_precision high
    variable lambda_1
    variable lambda_2
    variable P(n,n) symmetric 

    L1=[A'*P + P*A + 2*alpha*P,    P*B;
        B'*P,                      zeros(nu)];
    L2=[C';D']*[lambda_1*M,zeros(2*nu);zeros(2*nu),lambda_2*M]*[C,D];

    LMI=L1+L2;

    minimize 1; 
    subject to:
    P >= tol*eye(n);
    P <= cond_tol*tol*eye(n);
    LMI<=-tol*eye(n+nu);
    lambda_1>=tol;
    lambda_2>=tol;
    cvx_end 
    if strcmp('Solved',cvx_status)
        status=true;
    end
end
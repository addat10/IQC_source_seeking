function [status,P]=verify_exp_stab_ZF_basis_3(G,alpha,sec_1,sec_2,rho,n_psi,cond_tol,tol)
% This function runs the analysis LMI with cvx and returns the status and
% the storage function matrix P. The multipliers used are as described in
% Veenman and Scherer
    
    [A_G,B_G,C_G,D_G]=ssdata(G);
    
    % Negative feedback(if gradient is the delta block, we want grad descent)
    B_G=-B_G;
    
    nx=size(A_G,1);
    [ny,nu]=size(D_G);    
    if nu~=ny
        error('Delta block must be square')
    end
    n_delta=nu; 
    
    % Define the multiplier psi as in Scherer and Veenman  
    if n_psi>=1 % Dynamic Multiplier
        A_psi=rho*diag(ones(1,n_psi))+diag(ones(1,n_psi-1),-1);
        B_psi=[1;zeros(n_psi-1,1)];
        C_psi=[zeros(1,n_psi);eye(n_psi)];
        D_psi=[1;zeros(n_psi,1)];
        % Repeat along all channels in the delta block
        psi=ss(kron(A_psi,eye(nu)),kron(B_psi,eye(nu)),kron(C_psi,eye(nu)),kron(D_psi,eye(nu)));
    else
        D_psi=[1;zeros(n_psi,1)];
        psi=ss([],[],[],kron(D_psi,eye(nu)));
    end
    

    % Define the system psi_tilde for the positivity constraint
    psi_b_til=[];
    if n_psi>=1
         s=tf('s');        
        for i=1:n_psi
            basis=s^(i-1)/(s-rho)^(i-1);
            psi_b_til=[psi_b_til;basis];
        end    
    end
    
    % Build Psi*[G;I] with Dynamic Multiplier
    B_in=[  -sec_1*eye(n_delta), 1*eye(n_delta);...
            sec_2*eye(n_delta),-1*eye(n_delta)];
    GI=B_in*ss(A_G,B_G,[C_G;zeros(n_delta,nx)],[D_G;eye(n_delta)]);
    PSI_GI=[psi,zeros((1+n_psi)*nu,nu);zeros((1+n_psi)*nu,nu),psi]*GI;
    
    % Run the exponenetial analysis LMI
    [status,P]=verify_exp_stab_alpha_ZF_basis_LMI(PSI_GI,psi,psi_b_til,alpha,cond_tol,tol);
    
end
function [status,X]=verify_exp_stab_alpha_ZF_basis_LMI(PSI_GI,psi,psi_b_til,alpha,cond_tol,tol)
 % This function runs the exp-stab analysis KYP Lemma LMI
    status=false;
    
    % Get state-space matrices and dimensions
    [A,B,C,D]=ssdata(PSI_GI);
    [Av,Bv,~,~]=ssdata(minreal(psi));
    % Get sizes
    [n,nu]=size(B);
    n_psi=size(Av,1);
    n_psi_ch=n_psi/nu;
    if n_psi_ch>=1
        [Av_til,Bv_til,Cv_til,Dv_til]=ssdata(minreal(psi_b_til));    
        n_psi_til=size(Av_til,1);
    end
    % Build R_nu as in Scherer and Veenman
    Rv=diag(sqrt(factorial((0:n_psi_ch-1))));
    
    % Start describing the optimization problem
    cvx_begin sdp 
    cvx_precision high
    
    % Define the variables
    variable X(n,n) symmetric 
    variable P0
    if n_psi_ch>=1
        variable P1(1,n_psi_ch)
        variable P2(1,n_psi_ch)
        variable P3(1,n_psi_ch)
        variable P4(1,n_psi_ch)
        if nu==3
            % Additional variables for channel 2
            variable P0_2            
            variable P1_2(1,n_psi_ch)
            variable P2_2(1,n_psi_ch)
            variable P3_2(1,n_psi_ch)
            variable P4_2(1,n_psi_ch)
            
            % Additional variables for channel 3
            variable P0_3
            variable P1_3(1,n_psi_ch)
            variable P2_3(1,n_psi_ch)
            variable P3_3(1,n_psi_ch)
            variable P4_3(1,n_psi_ch)
        end        
    end
    
    % Define LMIs now
    
    % Common part of the main LMI
    L1=[A'*X + X*A + 2*alpha*X,    X*B;
        B'*X,                      zeros(nu)];
    
    % Multiplier part of the LMI
    if n_psi_ch==0
        P_12=P0;
    elseif n_psi_ch>=1
        P_12=[  P0, P4-P3;
                P2'-P1',    zeros(n_psi_ch)];
        if nu==3
        P_12=[diag([P0,P0_2,P0_3]),   diag([P4-P3,P4_2-P3_2,P4_3-P3_3]);
              diag([P2'-P1',P2_2'-P1_2',P2_3'-P1_3'])  ,    zeros(n_psi)];
        end
    end
    L2=[C';D']*[zeros(n_psi+nu),P_12;P_12',zeros(n_psi+nu)]*[C,D];

    LMI=L1+L2;
    
    % LMIs for positivity constraint
    if n_psi_ch==1
        % First order multipliers
        MAT=Rv*Dv_til;        
        pos_Lemma=@(P)(MAT'*diag(P)*MAT); 
    elseif n_psi_ch>=2
        % Higher order multipliers
        variable X1(n_psi_til,n_psi_til) symmetric
        variable X2(n_psi_til,n_psi_til) symmetric
        variable X3(n_psi_til,n_psi_til) symmetric
        variable X4(n_psi_til,n_psi_til) symmetric    
        MAT=[eye(n_psi_til),    zeros(n_psi_til,nu);...
            Av_til,             Bv_til;...
            Rv*Cv_til,          Rv*Dv_til];
        
        pos_Lemma=@(X,P)(MAT'*...
                        [zeros(n_psi_til),  X,                              zeros(n_psi_til,n_psi);...
                        X,                  zeros(n_psi_til),               zeros(n_psi_til,n_psi);...
                        zeros(n_psi,n_psi_til),   zeros(n_psi,n_psi_til),   diag(P)]*...
                        MAT);
    end
    
    minimize 1; 
    subject to:
    X >= tol*eye(n);
    X <= cond_tol*tol*eye(n);
    LMI<=-tol*eye(n+nu); 
    
    % L1 Norm constraint
    if n_psi_ch==0
        kron(eye(nu),P0)>=tol*eye(nu);
    elseif n_psi_ch>=1 
        kron(eye(nu),P0) + (P1+P2+P3+P4)*inv(Av)*Bv >=tol*eye(nu);
        if nu==3
            kron(eye(nu),P0_2) + (P1_2+P2_2+P3_2+P4_2)*inv(Av)*Bv >=tol*eye(nu);
            kron(eye(nu),P0_3) + (P1_3+P2_3+P3_3+P4_3)*inv(Av)*Bv >=tol*eye(nu);
        end        
    end
    
    % Positivity constraint
    if n_psi_ch==1        
        pos_Lemma(P1)>=tol*eye(n_psi_ch);
        
        pos_Lemma(P2)>=tol*eye(n_psi_ch);
        
        pos_Lemma(P3)>=tol*eye(n_psi_ch);
        
        pos_Lemma(P4)>=tol*eye(n_psi_ch);
        if nu==3
            pos_Lemma(P1_2)>=tol*eye(n_psi_ch);
            pos_Lemma(P2_2)>=tol*eye(n_psi_ch);
            pos_Lemma(P3_2)>=tol*eye(n_psi_ch);
            pos_Lemma(P4_2)>=tol*eye(n_psi_ch);
            
            pos_Lemma(P1_3)>=tol*eye(n_psi_ch);
            pos_Lemma(P2_3)>=tol*eye(n_psi_ch);
            pos_Lemma(P3_3)>=tol*eye(n_psi_ch);
            pos_Lemma(P4_3)>=tol*eye(n_psi_ch);
        end
    elseif n_psi_ch>=2
        X1>=tol*eye(n_psi_til)
        pos_Lemma(X1,P1)>=tol*eye(n_psi_til+nu);

        X2>=tol*eye(n_psi_til)
        %pos_Lemma(X2,P2)>=tol*eye(n_psi_til+nu);

        X3>=tol*eye(n_psi_til)
        pos_Lemma(X3,P3)>=tol*eye(n_psi_til+nu);

        X4>=tol*eye(n_psi_til)
        %pos_Lemma(X4,P4)>=tol*eye(n_psi_til+nu);
    end
    
    % Impose these constraints if uncertainty phi is not odd
%     if n_psi_ch>=1
%         P2==zeros(1,n_psi_ch);
%         P4==zeros(1,n_psi_ch);
%     end
    
    
    % Impose these constraints if only searching for causal multipliers
     %P1==P2;
    
    % Impose these constraints if only searching for anti-causal multipliers
     %P3==P4;
    
    cvx_end 
    if strcmp('Solved',cvx_status)
        status=true;
    end
end
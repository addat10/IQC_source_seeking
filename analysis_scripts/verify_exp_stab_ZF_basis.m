function [status,P]=verify_exp_stab_ZF_basis(G_veh,alpha,sec_1,sec_2,odd_flag,causal_flag,rho,n_psi,cond_tol,tol)
% This function runs the analysis LMI with cvx and returns the status and
% the storage function matrix P with Zames Falb Multipliers
    
    [A_G,B_G,C_G,D_G]=ssdata(G_veh);
    
    % Negative feedback(if gradient is the delta block, we want grad descent)
    B_G=-B_G;
    
    nx=size(A_G,1);
    [ny,nu]=size(D_G);    
    if nu~=ny
        error('Delta block must be square')
    end
    n_delta=nu; 
    
    % Define the multiplier psi as in Scherer and Veenman  
    if n_psi>=1
        A_psi_alpha=(rho-2*alpha)*diag(ones(1,n_psi))+diag(ones(1,n_psi-1),-1);
        B_psi_alpha=[1;zeros(n_psi-1,1)];
        C_psi_alpha=[zeros(1,n_psi);eye(n_psi)];
        D_psi_alpha=[1;zeros(n_psi,1)];
        psi_alpha=ss(kron(A_psi_alpha,eye(nu)),kron(B_psi_alpha,eye(nu)),kron(C_psi_alpha,eye(nu)),kron(D_psi_alpha,eye(nu)));
    else
        D_psi_alpha=[1;zeros(n_psi,1)];
        psi_alpha=ss([],[],[],kron(D_psi_alpha,eye(nu)));
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
    PSI_GI=[psi_alpha,zeros((1+n_psi)*nu,nu);zeros((1+n_psi)*nu,nu),psi_alpha]*GI;
    
    % Run the exponenetial analysis LMI
    [status,P]=verify_exp_stab_alpha_ZF_basis_LMI(PSI_GI,psi_alpha,psi_b_til,alpha,odd_flag,causal_flag,cond_tol,tol);
    
end
function [status,X]=verify_exp_stab_alpha_ZF_basis_LMI(PSI_GI,psi_alpha,psi_b_til,alpha,odd_flag,causal_flag,cond_tol,tol)
 % This function runs the exp-stab analysis KYP Lemma LMI
    status=false;
    % Get state-space matrices and dimensions
    [A,B,C,D]=ssdata(PSI_GI);
    [Av_alpha,Bv_alpha,~,~]=ssdata(minreal(psi_alpha));
    
    Av=Av_alpha+2*alpha*eye(size(Av_alpha,1));
    Bv=Bv_alpha;
    
    [n,nu]=size(B);
    n_psi=size(Av,1);
    n_psi_ch=n_psi/nu;
    % Get R_nu as in Scherer and Veenman
    Rv=diag(sqrt(factorial((0:n_psi_ch-1))));
    
    if n_psi_ch>=1
        [Av_til,Bv_til,Cv_til,Dv_til]=ssdata(minreal(psi_b_til));    
        n_psi_til=size(Av_til,1);
    end
    

    cvx_begin sdp 
    cvx_precision high
    
    variable X(n,n) symmetric 
    variable P0
    variable P1(1,n_psi_ch)
    variable P2(1,n_psi_ch)
    variable P3(1,n_psi_ch)
    variable P4(1,n_psi_ch)

    if n_psi_ch>=2       
        variable X1(n_psi_til,n_psi_til) symmetric
        variable X2(n_psi_til,n_psi_til) symmetric
        variable X3(n_psi_til,n_psi_til) symmetric
        variable X4(n_psi_til,n_psi_til) symmetric
    end
    
    P_12=[  P0, P4-P3;
            P2'-P1',    zeros(n_psi_ch)];

    L1=[A'*X + X*A + 2*alpha*X,    X*B;
        B'*X,                      zeros(nu)];
    L2=[C';D']*[zeros(n_psi+nu),kron(P_12,eye(nu));kron(P_12',eye(nu)),zeros(n_psi+nu)]*[C,D];

    LMI=L1+L2;
    
    if n_psi_ch>=2
        MAT=[eye(n_psi_til),    zeros(n_psi_til,1);...
            Av_til,             Bv_til;...
            Rv*Cv_til,          Rv*Dv_til];
        
        pos_Lemma=@(X,P)(MAT'*...
                        [zeros(n_psi_til),  X,                              zeros(n_psi_til,n_psi);...
                        X,                  zeros(n_psi_til),               zeros(n_psi_til,n_psi);...
                        zeros(n_psi,n_psi_til),   zeros(n_psi,n_psi_til),   diag(P)]*...
                        MAT);
    elseif n_psi_ch==1
        MAT=[Rv*Dv_til];        
        pos_Lemma=@(P)(MAT'*diag(P)*MAT);    
    end
    

    minimize 1; 
    subject to:
    X >= tol*eye(n);
    X <= cond_tol*tol*eye(n);
    LMI<=-tol*eye(n+nu); 
    
    % Norm constraint
    if n_psi_ch>=1        
        kron(eye(nu),P0) + kron((P1+P2+P3+P4),eye(nu))*inv(Av)*Bv >=tol*eye(nu);
    elseif n_psi_ch==0
        kron(eye(nu),P0)>=tol*eye(nu);
    end
    
    % Positivity constraint for the two cases of odd and non-odd uncertainty 
    % Case 1: If the uncertainty is odd
    if odd_flag==1 
        if n_psi_ch==1
            pos_Lemma(P1)>=tol*eye(n_psi_ch);
            pos_Lemma(P3)>=tol*eye(n_psi_ch);
            pos_Lemma(P2)>=tol*eye(n_psi);
            pos_Lemma(P4)>=tol*eye(n_psi);
        elseif n_psi_ch>=2
            X1>=tol*eye(n_psi_til)
            pos_Lemma(X1,P1)>=tol*eye(n_psi_til+nu);
            X3>=tol*eye(n_psi_til)
            pos_Lemma(X3,P3)>=tol*eye(n_psi_til+nu);
            X2>=tol*eye(n_psi_til)
            pos_Lemma(X2,P2)>=tol*eye(n_psi_til+nu);
            X4>=tol*eye(n_psi_til)
            pos_Lemma(X4,P4)>=tol*eye(n_psi_til+nu);
        end
        if causal_flag==1
        % Impose these constraints if only searching for causal multipliers
            P1==P2;
        end
        if causal_flag==-1
        % Impose these constraints if only searching for anti-causal multipliers
            P3==P4;
        end
    end
    %  Case 2: if the uncertainty is not odd
    if odd_flag==0        
        if causal_flag==1
            P2==zeros(1,n_psi_ch);
            P4==zeros(1,n_psi_ch);
            P1==P2;
            if n_psi_ch==1               
                pos_Lemma(P3)>=tol*eye(n_psi_ch);            
            elseif n_psi_ch>=2                
                X3>=tol*eye(n_psi_til)
                pos_Lemma(X3,P3)>=tol*eye(n_psi_til+nu);            
            end            
        end
        if causal_flag==-1            
            P2==zeros(1,n_psi_ch);
            P4==zeros(1,n_psi_ch);
            P3==P4;
            if n_psi_ch==1               
                pos_Lemma(P1)>=tol*eye(n_psi_ch);            
            elseif n_psi_ch>=2                
                X1>=tol*eye(n_psi_til)
                pos_Lemma(X1,P1)>=tol*eye(n_psi_til+nu);            
            end            
        end        
    end  
    
    
    cvx_end 
    if strcmp('Solved',cvx_status)
        status=true;
    end
end
function [alpha_best,P_ret]=bisection_exponent(G_veh,m,L,alpha_lims,cond_tol,cvx_tol,bisect_tol,multiplier_flag)
% This function runs a bisection algorithm for obtaining the best
% exponent alpha by calling the analysis routine for a fixed alpha
% multiple times.
    
    % Analyze for alpha=alpha_low
    switch multiplier_flag
        case 1  % Circle Criterion
            [status,P_ret]=verify_exp_stab_CC(G_veh,alpha_lims(1),m,L,cond_tol,cvx_tol);
        case 2  % Full block circle criterion
            [status,P_ret]=verify_exp_stab_FBCC(G_veh,alpha_lims(1),m,L,cond_tol,cvx_tol);
        case 3  % Zames Falb Multipliers with CC
            [status,P_ret]=verify_exp_stab_ZF(G_veh,alpha_lims(1),m,L,cond_tol,cvx_tol);
        case 4 % Zames Falb Multipliers with FBCC
            [status,P_ret]=verify_exp_stab_ZF_FBCC(G_veh,alpha_lims(1),m,L,cond_tol,cvx_tol);
        case 5 % Zames Falb multipliers with non-causal multilpiers
            [status,P_ret]=verify_exp_stab_ZF_NC(G_veh,alpha_lims(1),m,L,cond_tol,cvx_tol);
        case 6
            rho=-1;
            psi_order=1;
            odd_flag=0;
            causal_flag=0; % 1: causal, -1:anti-causal, 0:non-causal
            [status,P_ret]=verify_exp_stab_ZF_basis(G_veh,alpha_lims(1),m,L,odd_flag,causal_flag,rho,psi_order,cond_tol,cvx_tol);

    end

    if ~status
        %error('Infeasible for alpha_low. Choose a smaller alpha_low')
        alpha_best=-1;
        return
    end
    
    % Analyze for alpha=alpha_high
    switch multiplier_flag
        case 1  % Circle Criterion
            [status,P]=verify_exp_stab_CC(G_veh,alpha_lims(2),m,L,cond_tol,cvx_tol);
        case 2  % Full block circle criterion
            [status,P]=verify_exp_stab_FBCC(G_veh,alpha_lims(2),m,L,cond_tol,cvx_tol);
        case 3  % Zames Falb Multipliers with CC
            [status,P]=verify_exp_stab_ZF(G_veh,alpha_lims(2),m,L,cond_tol,cvx_tol);
        case 4 % Zames Falb Multipliers with FBCC
            [status,P]=verify_exp_stab_ZF_FBCC(G_veh,alpha_lims(2),m,L,cond_tol,cvx_tol);
        case 5
            [status,P]=verify_exp_stab_ZF_NC(G_veh,alpha_lims(2),m,L,cond_tol,cvx_tol);
        case 6
            [status,P]=verify_exp_stab_ZF_basis(G_veh,alpha_lims(2),m,L,odd_flag,causal_flag,rho,psi_order,cond_tol,cvx_tol);
    end    
    if status        
        alpha_best=alpha_lims(2); % Return alpha_high if feasible
        P_ret=P;
    else
    
    % Start the bisection
    while alpha_lims(2)-alpha_lims(1)>bisect_tol
        alpha_mid=mean(alpha_lims);           
        switch multiplier_flag
            case 1  % Circle Criterion
                [status,P]=verify_exp_stab_CC(G_veh,alpha_mid,m,L,cond_tol,cvx_tol);
            case 2  % Full block circle criterion
                [status,P]=verify_exp_stab_FBCC(G_veh,alpha_mid,m,L,cond_tol,cvx_tol);
            case 3  % Zames Falb Multipliers
                [status,P]=verify_exp_stab_ZF(G_veh,alpha_mid,m,L,cond_tol,cvx_tol);
            case 4 % Zames Falb Multipliers with FBCC
                [status,P]=verify_exp_stab_ZF_FBCC(G_veh,alpha_mid,m,L,cond_tol,cvx_tol);
            case 5
                [status,P]=verify_exp_stab_ZF_NC(G_veh,alpha_mid,m,L,cond_tol,cvx_tol);
            case 6
                [status,P]=verify_exp_stab_ZF_basis(G_veh,alpha_mid,m,L,odd_flag,causal_flag,rho,psi_order,cond_tol,cvx_tol);
        end            
        if status
            alpha_lims(1)=alpha_mid;
            P_ret=P;
        else
            alpha_lims(2)=alpha_mid;        
        end
    end
    alpha_best=alpha_lims(1);
    end
end
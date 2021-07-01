function [alpha_best,P]=bisection_exponent(G_veh,m,L,alpha_lims,cond_tol,cvx_tol,bisect_tol,multiplier_flag)
% This function runs a bisection algorithm for obtaining the best
% exponent alpha by calling the analysis routine for a fixed alpha
% multiple times.
    
    % Analyze for alpha=alpha_low
    switch multiplier_flag
        case 1  % Circle Criterion
            [status,~]=verify_exp_stab_CC(G_veh,alpha_lims(1),m,L,cond_tol,cvx_tol);
        case 2  % Full block circle criterion
            [status,~]=verify_exp_stab_FBCC(G_veh,alpha_lims(1),m,L,cond_tol,cvx_tol);
        case 3  % Zames Falb Multipliers
            [status,P]=verify_exp_stab_ZF(G_veh,alpha_lims(1),m,L,cond_tol,cvx_tol);
    end

    if ~status
        error('Infeasible for alpha_low. Choose a smaller alpha_low')
    end
    
    % Analyze for alpha=alpha_high
    switch multiplier_flag
        case 1  % Circle Criterion
            [status,P]=verify_exp_stab_CC(G_veh,alpha_lims(2),m,L,cond_tol,cvx_tol);
        case 2  % Full block circle criterion
            [status,P]=verify_exp_stab_FBCC(G_veh,alpha_lims(2),m,L,cond_tol,cvx_tol);
        case 3  % Zames Falb Multipliers
            [status,P]=verify_exp_stab_ZF(G_veh,alpha_lims(2),m,L,cond_tol,cvx_tol);
    end    
    if status        
        alpha_best=alpha_lims(2); % Return alpha_high if feasible
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
        end            
        if status
            alpha_lims(1)=alpha_mid;
        else
            alpha_lims(2)=alpha_mid;        
        end
    end
    alpha_best=alpha_lims(1);
    end
end
function [alpha_best,P]=bisection_exponent_CC(G_veh,m,L,alpha_lims,cond_tol,cvx_tol,bisect_tol,multiplier_flag)
% This function runs a bisection algorithm for obtaining the best
% exponent alpha by calling the analysis routine for a fixed alpha
% multiple times.
    
    % Analyze for alpha=alpha_low
    switch multiplier_flag
        case 1  [status,P]=verify_exp_stab_CC(G_veh,alpha_lims(1),m,L,cond_tol,cvx_tol);
        case 2
        case 3
    end

    if ~status
        error('Infeasible for alpha_low. Choose a smaller alpha_low')
    end
    
    % Analyze for alpha=alpha_high
    switch multiplier_flag
        case 1  [status,P]=verify_exp_stab_CC(G_veh,alpha_lims(2),m,L,cond_tol,cvx_tol);
        case 2
        case 3
    end
    
    if status        
        alpha_best=alpha_lims(2); % Return alpha_high if feasible
    else
        while alpha_lims(2)-alpha_lims(1)>bisect_tol
            alpha_mid=mean(alpha_lims);           
            switch multiplier_flag
                    case 1  [status,P]=verify_exp_stab_CC(G_veh,alpha_mid,m,L,cond_tol,cvx_tol);
                    case 2
                    case 3
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
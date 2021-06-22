function [alpha_best,P]=bisection_exponent(Psi_G_I,M,alpha_lims,cvx_tol,bisect_tol)
% This function runs a bisection algorithm for obtaining the best
% exponent alpha by calling the analysis routine for a fixed alpha
% multiple times.
    
    % Analyze for alpha=alpha_low
    [status,P]=verify_exp_stab(Psi_G_I,M,alpha_lims(1),cvx_tol);
    if ~status
        error('Infeasible for alpha_low. Choose a smaller alpha_low')
    end
    
    % Analyze for alpha=alpha_high
    [status,P]=verify_exp_stab(Psi_G_I,M,alpha_lims(2),cvx_tol);
    if status        
        alpha_best=alpha_lims(2); % Return alpha_high if feasible
    else
        while alpha_lims(2)-alpha_lims(1)>bisect_tol
            alpha_mid=mean(alpha_lims);
            [status,P]=verify_exp_stab(Psi_G_I,M,alpha_mid,cvx_tol);
            if status
                alpha_lims(1)=alpha_mid;
            else
                alpha_lims(2)=alpha_mid;        
            end
        end
        alpha_best=alpha_lims(1);
    end
end
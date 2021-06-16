function [alpha_best]=bisection_exponent(Psi_G_I,M,alpha_lims,cvx_tol,bisect_tol)
    if ~verify_exp_stab(Psi_G_I,M,alpha_lims(1),cvx_tol)
        error('Infeasible for alpha_low')
    end
    if verify_exp_stab(Psi_G_I,M,alpha_lims(2),cvx_tol)
        alpha_best=alpha_lims(2);
    else
        while alpha_lims(2)-alpha_lims(1)>bisect_tol
            alpha_mid=mean(alpha_lims);
            if verify_exp_stab(Psi_G_I,M,alpha_mid,cvx_tol)
                alpha_lims(1)=alpha_mid;
            else
                alpha_lims(2)=alpha_mid;        
            end
        end
        alpha_best=alpha_lims(1);
    end
end
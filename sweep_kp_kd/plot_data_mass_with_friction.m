%% Plot results parameter sweep
clear
clc
%% Performance curves
id=50; % represents the fixed value of a
figure()
plot_data_perf(id,'.\data_mass_with_friction\coarse_grid_mult_flag_1')
% hold on    
% plot_data_perf(id,'.\data_mass_with_friction\mult_flag_non_causal_6')
% plot_data_perf(id,'.\data_mass_with_friction\backup\mult_flag_causal_6')
% plot_data_perf(id,'.\data_mass_with_friction\backup\mult_flag_anti_causal_6')
legend('CC','ZFb','ZFb causal','ZFb anti-causal')
xlabel('kd')
ylabel('alpha')
ylim([0,2])
title('Convergence rates(exponents)')
%% Stability regions
% figure()
% plot_data_stab('.\data_odd_phi\mult_flag_1')
% hold on
% plot_data_stab('.\data_odd_phi_alpha\mult_flag_causal_6')
% plot_data_stab('.\data_odd_phi_alpha\mult_flag_anti_causal_6')
% plot_data_stab('.\data_odd_phi_alpha\mult_flag_non_causal_6')
% legend('CC','ZFb causal','ZFb anti causal','ZFb non causal')
% xlabel('a')
% ylabel('b')
% title('stability region')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function []=plot_data_perf(id,file)
    save_path=['.\',file];
    data=load(save_path,'alpha_best','kd');
    alpha_best=data.alpha_best;
    kd=data.kd;    
    plot(kd,alpha_best(id,:))    
end
function []=plot_data_stab(file)
    save_path=['.\',file];
    data=load(save_path,'alpha_best','kp','kd');
    alpha_best=data.alpha_best;
    kp=data.kp;
    kd=data.kd;    
       
    plot(a_vec,b_vec)    
end
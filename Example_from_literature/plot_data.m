%% Plot results for Example 7 20 from Scherer and Weilands LMI notes
clear
clc
%% Performance curves
id=1; % represents the fixed value of a
figure()
plot_data_perf('mult_flag_1')
hold on    
plot_data_perf('mult_flag_6')
legend('CC','ZFb')
xlabel('b')
ylabel('alpha')
ylim([0,0.4])
title('Convergence rates(exponents)')
%% Stability regions
figure()
plot_data_stab('mult_flag_1')
hold on
plot_data_stab('mult_flag_6')
legend('CC','ZFb')
xlabel('a')
ylabel('b')
title('stability region')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function []=plot_data_perf(file)
    save_path=['.\',file];
    data=load(save_path,'alpha_best','L');
    alpha_best=data.alpha_best;
    L=data.L;    
    plot(L,alpha_best)    
end
function []=plot_data_stab(file)
    save_path=['.\',file];
    data=load(save_path,'alpha_best','L');
    alpha_best=data.alpha_best;
    L=data.L; 
    id=max(find(alpha_best(1,:)>=0));
    b_vec=L(1,id);
    plot(1,b_vec,'*')    
end
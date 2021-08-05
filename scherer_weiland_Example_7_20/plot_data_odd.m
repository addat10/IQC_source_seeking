%% Plot results for Example 7 20 from Scherer and Weilands LMI notes
clear
clc
%% Performance curves
id=8; % represents the fixed value of a
figure()
plot_scherer_weiland_data_perf(id,'.\data_odd_phi\mult_flag_1')
hold on    
plot_scherer_weiland_data_perf(id,'.\data_odd_phi\mult_flag_causal_6')
plot_scherer_weiland_data_perf(id,'.\data_odd_phi\mult_flag_anti_causal_6')
plot_scherer_weiland_data_perf(id,'.\data_odd_phi\mult_flag_non_causal_6')
legend('CC','ZFb causal','ZFb anti causal','ZFb non causal')
xlabel('b')
ylabel('alpha')
ylim([0,0.4])
title('Convergence rates(exponents) for Example 7.20')
%% Stability regions
figure()
plot_scherer_weiland_data_stab('.\data_odd_phi\mult_flag_1')
hold on
plot_scherer_weiland_data_stab('.\data_odd_phi\mult_flag_causal_6')
plot_scherer_weiland_data_stab('.\data_odd_phi\mult_flag_anti_causal_6')
plot_scherer_weiland_data_stab('.\data_odd_phi\mult_flag_non_causal_6')
legend('CC','ZFb causal','ZFb anti causal','ZFb non causal')
xlabel('a')
ylabel('b')
title('stability region')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function []=plot_scherer_weiland_data_perf(id,file)
    save_path=['.\',file];
    data=load(save_path,'alpha_best','L');
    alpha_best=data.alpha_best;
    L=data.L;    
    plot(L,alpha_best(id,:))    
end
function []=plot_scherer_weiland_data_stab(file)
    save_path=['.\',file];
    data=load(save_path,'alpha_best','L');
    alpha_best=data.alpha_best;
    L=data.L;  
    a_vec=0.4:0.2:1.8; 
    b_vec=zeros(1,size(a_vec,2));
    for i=1:size(a_vec,2)
        id=max(find(alpha_best(i,:)>=0));
        b_vec(1,i)=L(1,id);
    end    
    plot(a_vec,b_vec)    
end
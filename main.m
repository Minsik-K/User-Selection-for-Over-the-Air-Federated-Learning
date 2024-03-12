%----------------------------------------------------------
% Implementation of the paper [1]
% 
% [1] Minsik Kim, Lee Swindlehurst and Daeyoung Park,
%    " Beamforming Vector Design and Device Selection in Over-the-Air Federated Learning".
%
% This script reproduces the Figs. 2, 3, 4. for proposed algorithm in [1]
%
%---------------------------------------------------------- 

%% Fig. 2. The number of selected devices versus the MSE requirement
clear all; clc
% close all
% randn('seed', 42);

disp("start time : " + string(datetime('now')));
disp(" ----------------------------- ");

N = 6;  % N: Number of BS antennas
K = 60; % K: Number of devices

% MSE threshold gamma (dB)
MSE_threshold_dB =  -12:4:20;

devices_CVX = zeros(length(MSE_threshold_dB),1);
devices_Subgrad = zeros(length(MSE_threshold_dB),1);
devices_RB = zeros(length(MSE_threshold_dB),1);

simulation_Iter = 100;
% parfor j = 1 : length(MSE_threshold_dB)
for j = 1 : length(MSE_threshold_dB)

    MSE_threshold = 10^(MSE_threshold_dB(j)/10);
    
    disp("Current MSE: " + MSE_threshold_dB(j))
    disp("Current time: " + datestr(datetime('now')))
    disp(" ----------------------------- ");
        
    for iter=1:simulation_Iter
        H = (randn(N,K)+1i*randn(N,K)) / sqrt(2);
        
% % imbalaced channel gain
%         temp = rand(1,60)*-3;
%         temp = sqrt(10.^(temp./10));
%         H = H.*temp;

        tic % CVX
        [m_cvx, num_of_users_CVX, ~ ] = func_device_selection(H, MSE_threshold,1);            
        toc;
        
        tic % Subgradient
        [m_subgrad, num_of_users_Subgrad, ~ ] = func_device_selection(H, MSE_threshold,0);    
        toc;
                
        tic % RB
        [m_rb, num_of_users_RB, ~ ] = func_device_selection(H, MSE_threshold,4); 
        toc;
        
        % Check the performance
%         [MSE_threshold_dB(j) dev_CVX(iter) dev_Subgrad(iter) dev_DC(iter) dev_RB(iter) dev_SDR(iter)]

        devices_CVX(j) = devices_CVX(j) + num_of_users_CVX / simulation_Iter;
        devices_Subgrad(j) = devices_Subgrad(j) + num_of_users_Subgrad / simulation_Iter;
        devices_RB(j) = devices_RB(j) + num_of_users_RB / simulation_Iter;
    end
end
filename = "results\Selected_device_vs_MSE_threshold_"+"N("+string(N)+")K("+string(K)+")";
filename = filename + string(datetime('now','Format','yy_MM_dd_(HH_mm)'));
save(filename)

figure;
plot(MSE_threshold_dB, devices_CVX, 'o-'); hold on; grid;
plot(MSE_threshold_dB, devices_Subgrad, 's-');
plot(MSE_threshold_dB, devices_RB, '*-'); % Random beamforming
plot(MSE_threshold_dB, 60.*exp(-1./(10.^(MSE_threshold_dB/10))), '-.'); % Random beamforming

xlabel('MSE requirement (dB)','Interpreter','latex');
ylabel('Number of selected devices','Interpreter','latex');
legend('CVX','Subgradient','Random Beamforming','Random Beamforming(Analysis)', ...
    'location','northwest','Interpreter','latex')

figure_name = filename;
savefig(figure_name)

disp("finish time : " + string(datetime('now')));



%% Fig. 3. The number of selected devices versus the number of devices
clear all; clc
% close all
% randn('seed',85);

disp("start time : " + string(datetime('now')));
disp(" ----------------------------- ");

N = 6; % N: Number of BS antennas;
Ks = 30:10:100; % K: Number of devices

% MSE threshold gamma (dB)
MSE_threshold_dB = 4;
MSE_threshold = 10^(MSE_threshold_dB/10);

devices_CVX = zeros(length(Ks),1);
devices_Subgrad = zeros(length(Ks),1);
devices_RB = zeros(length(Ks),1);

simulation_Iter = 100;
% parfor j=1:length(Ks)
for j=1:length(Ks)
    K=Ks(j);
    
    disp("Current K: " + Ks(j))
    disp("Current time: " + datestr(datetime('now')))
    disp(" ----------------------------- ");

    for iter=1:simulation_Iter
        H = (randn(N,K)+1i*randn(N,K)) / sqrt(2);
        
        tic % CVX
        [m_cvx, num_of_users_CVX, ~ ] = func_device_selection(H, MSE_threshold,1);     
        toc;
        
        tic % Subgradient
        [m_subgrad, num_of_users_Subgrad, ~ ] = func_device_selection(H, MSE_threshold,0); 
        toc;
         
        tic % RB
        [m_rb, num_of_users_RB, ~ ] = func_device_selection(H, MSE_threshold,4); 
        toc;
        
        devices_CVX(j) = devices_CVX(j) + num_of_users_CVX / simulation_Iter;
        devices_Subgrad(j) = devices_Subgrad(j) + num_of_users_Subgrad / simulation_Iter;
        devices_RB(j) = devices_RB(j) + num_of_users_RB / simulation_Iter;
    end
end
filename = "results\Selected_device_vs_Total_number_of_Dev_"+"N("+string(N)+")MSE("+string(MSE_threshold_dB)+")";
filename = filename + string(datetime('now','Format','yy_MM_dd_(HH_mm)'));
save(filename)

figure;
plot(Ks, Ks, '--'); hold on; grid;          % Oracle
plot(Ks, devices_CVX, 'o-');
plot(Ks, devices_Subgrad, 's-');
plot(Ks, Ks.*exp(-1/MSE_threshold), '-.')   % Random beamforming(Analysis)
plot(Ks, devices_RB, '*-');
xlabel('Number of devices','Interpreter','latex');
ylabel('Number of selected devices','Interpreter','latex');
legend('Oracle-Aircomp','CVX','Subgrad','Random beamforming(Analysis)', ...
    'Random beamforming','location','southeast','Interpreter','latex')

figure_name = filename;
savefig(figure_name)

disp("finish time : " + string(datetime('now')));


%% Fig. 4. The relationship between the scaled MSE upper-bound and the MSE requirement.
clear all; clc
% close all

% randn('seed',85);

disp("start time : " + string(datetime('now')));
disp(" ----------------------------- ");

N = 6;  % N: number of receive antennas
K = 60; % K: number of devices

% MSE threshold gamma (dB)
MSE_threshold_dB = -8:2:10;

sub_devices_inv2=zeros(length(MSE_threshold_dB),1);
cvx_devices_inv2=zeros(length(MSE_threshold_dB),1);
rb_devices_inv2=zeros(length(MSE_threshold_dB),1);
    
simulation_Iter = 100;
for j=1:length(MSE_threshold_dB)
    MSE_threshold = 10^(MSE_threshold_dB(j)/10);
    
        disp("Current MSE: " + MSE_threshold_dB(j))
        disp("Current time: " + datestr(datetime('now')))
        disp(" ----------------------------- ");
    
    parfor iter=1:simulation_Iter
        H = (randn(N,K)+1i*randn(N,K)) / sqrt(2);
        
        % subgradient
        [~, num_of_users, ~] = func_device_selection(H, MSE_threshold, 0);
        sub_dev(iter) = num_of_users+10^-20;
        
        % CVX
        [~, num_of_users, ~] = func_device_selection(H, MSE_threshold, 1);
        cvx_dev(iter) = num_of_users+10^-20;
        
        % Random beamforming
        [~, num_of_users, ~] = func_device_selection(H, MSE_threshold, 4);
        rb_dev(iter) = num_of_users+10^-20;
    end
    
    sub_devices_inv2(j)=mean(1./sub_dev.^2);
    cvx_devices_inv2(j)=mean(1./cvx_dev.^2);
    rb_devices_inv2(j)=mean(1./rb_dev.^2);
end
filename = "results\MSE_upper-bound_vs_the_MSE_threshold_"+"N("+string(N)+")K("+string(K)+")";
filename = filename + string(datetime('now','Format','yy_MM_dd_(HH_mm)'));
save(filename)

figure;
plot(MSE_threshold_dB, MSE_threshold_dB + 10*log10(rb_devices_inv2'), '<-'); hold on; grid;
plot(MSE_threshold_dB, MSE_threshold_dB + 10*log10(sub_devices_inv2'), '^-')
plot(MSE_threshold_dB, MSE_threshold_dB + 10*log10(cvx_devices_inv2'), 'o-')
xlabel('MSE threshold ($\gamma$) (dB)','Interpreter','latex');
ylabel('$\gamma$ E$[1/S^2]$ (dB)','Interpreter','latex');
axis([-8 inf -inf -10])
legend('Random beamforming','Subgrad','CVX','location','northeast','Interpreter','latex')

figure_name = filename;
savefig(figure_name)

disp("finish time : " + string(datetime('now')));
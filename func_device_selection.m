function [m_final, num_of_users, active_user]= func_device_selection(H_original, MSE_threshold, mode, sigma, m_old)
%----------------------------------------------------------
%  
% mode 0 : Subgradient method   % min ||x||_1 subject to tau - |g_i^H * m|^2 <= x_i,   ||m||^2 <= 1
% mode 1 : CVX                  % min ||x||_1 subject to tau - |g_i^H * m|^2 <= x_i,   ||m||^2 <= 1
% mode 2 : CVX                  % min ||x||_1 subject to tau - |g_i^H * m|^2 <= x_i,   ||m||^2 <= 1, x>=0
% mode 3 : Heuristic method     % max sum_i log(1+ theta*|y_1|^2/tau)  sujectt to y=G^H*m, ||m||^2 <= 1
% mode 4 : Random beamforming   % lower bound on the performance of beamforming algorithms
%  
%---------------------------------------------------------- 

if nargin<3
    mode = 0;
end

if nargin<4
    sigma = 0;
end

K = size(H_original, 2); % K: Number of Users
N = size(H_original, 1); % K: Number of Ant.

% imperfect channel
G = sqrt(1-sigma^2)*H_original + sigma*(randn(N,K)+1i*randn(N,K)) / sqrt(2);

if nargin<4
    [m, num_of_users, active_user] = func_beam_design(G, MSE_threshold, mode);
else
    [m, num_of_users, active_user] = func_beam_design(G, MSE_threshold, mode, m_old);
end

m_final=m;

% Check active users
if norm(m_final)<0.1
    num_of_users=0;
    active_user=[];
else
    active=(abs(H_original'*m_final).^2>=1/MSE_threshold);
    active_user=(1:K);
    active_user=active_user(active);
    num_of_users=sum(active);
end



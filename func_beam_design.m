function [m, num_of_users, active_user] = func_beam_design(G, MSE_threshold, mode, m_l)
%----------------------------------------------------------
% Algorithm for following problem 
%
%   min ||x||_1
%   subject to tau - |g_i^H * m|^2 <= x_i
%              ||m||^2 <= 1
%
% mode 0 : Subgradient method   % min ||x||_1 subject to tau - |g_i^H * m|^2 <= x_i,   ||m||^2 <= 1
% mode 1 : CVX                  % min ||x||_1 subject to tau - |g_i^H * m|^2 <= x_i,   ||m||^2 <= 1
% mode 2 : CVX                  % min ||x||_1 subject to tau - |g_i^H * m|^2 <= x_i,   ||m||^2 <= 1, x>=0
% mode 3 : Heuristic method     % max sum_i log(1+ theta*|y_1|^2/tau)  sujectt to y=G^H*m, ||m||^2 <= 1
% mode 4 : Random beamforming   % lower bound on the performance of beamforming algorithms
%  
%---------------------------------------------------------- 

% G: imperfect channel
% G = sqrt(1-sigma^2)*H_original + sigma*(randn(N,K)+1i*randn(N,K)) / sqrt(2);

[N,K] = size(G); % N: Number of antennas   K: Number of devices

tau = 1/MSE_threshold ;
theta = 10^9;

if nargin<3
    mode = 0;   % subgradient method (proposed low-complexity method)
end

if nargin<4
    m_l = ones(N,1)/sqrt(N);
    
    max_users=sum(abs(G'*m_l).^2>= tau);
    
    for j=1:K
        m = G(:,j)/norm(G(:,j));
        users=sum(abs(G'*m).^2>= tau);
        if users>max_users
            max_users=users;
            m_l=m;
        end
    end
end
m=m_l;

% Stop criteria (if all devices are already selected)
active=(abs(G'*m).^2>=tau );
num_of_users=sum(active);

if num_of_users==K
    active_user=(1:K)';
    active_user=active_user(active);
    return;
end

Iter = 100;     % Maximum iteration
Iter_inner = 100;
eps = 10^-3;

switch mode
    
    case 0  % subgradient method      
        max_num_of_users=-1;
        for iter=1:Iter
           % Initialization
            m_old = m_l;
            m = m_l;
            
            lambda = ones(K,1);
            lambda = lambda - 0.02*(tau-abs(G'*m_l).^2>0);
            
            for inner=1:1000
                
                % Update beamforming vector m
                m = G*(diag(lambda)*(G'*m_l));
                m = m/norm(m);   % projection m onto set { m | ||m||^2<=1}
                
                % |G'm|^2 >= |G'm^(l)|^2 + 2*Re{(G'm^(l))^T*(G'(m-m^(l))}   in this code
                %          = 2*Re{(G'm^(l))^T*(G'm)} - |G'm^(l)|^2          in the paper
                
                % constraint: tau - |G'm|^2 < x
                
                % Update MSE violations
                x = max(tau - abs(G'*m_l).^2 - 2*real( conj((G'*m_l)) .* (G'*(m-m_l))), 0);
                
                % Update dual variable based on projected subgradient method
                xi = 0.1/sqrt(inner);  % Step size of lambda
                lambda = max(lambda - xi*(x - tau + abs(G'*m_l).^2 + 2*real( conj((G'*m_l)) .* (G'*(m-m_l)))), 10^-10 );
                
                % Check the duality gap
                duality_gap = abs(lambda'*(tau + abs(G'*m_l).^2 - 2*real( conj((G'*m_l)) .* (G'*(m-m_l))) -x));
                
                if norm(m-m_old)<eps
                    break;
                else
                    m_old = m;
                end
                
                if inner==Iter_inner
                    duality_gap;
                end
            end
            
            num_of_users=sum(abs(G'*m).^2>=tau);           
            if num_of_users>max_num_of_users
                m_final = m;
                max_num_of_users = num_of_users;
            end            
            
            % Stop criteria for MM algorithm
            if norm(m-m_l)<eps || max_num_of_users==K
                m=m_final;
                break;
            else
                m_l = m;
            end
        end   
        m = m_final;
        
    case 1  % CVX     min ||x||_1 subject to tau - |g_i^H * m|^2 <= x_i,   ||m||^2 <= 1     
        max_num_of_users=-1;
        for iter=1:Iter
            
            cvx_begin quiet
            variable m(N) complex
            variable x(K)
            dual variables lambda{K}
            minimize norm(x,1)
            for i=1:K
                lambda{i}:    x(i) >= tau - abs(G(:,i)'*m_l)^2 - 2*real(m_l'*G(:,i)*G(:,i)'*(m-m_l));
            end
            norm(m) <= 1;
            cvx_end
            
            num_of_users=sum(abs(G'*m).^2>=tau);
            if num_of_users>max_num_of_users
                m_final = m;
                max_num_of_users = num_of_users;
            end
            
            % Stop criteria
            if norm(m-m_l)<eps || num_of_users==K
                break;
            else
                m_l = m;
            end
            
        end
        m = m_final;
        
    case 2  % CVX   min ||x||_1 subject to tau - |g_i^H * m|^2 <= x_i,   ||m||^2 <= 1, x>=0       
        max_num_of_users=-1;
        for iter=1:Iter
            
            cvx_begin quiet
            variable m(N) complex
            variable x(K)
            dual variables a{K};
            dual variables b{K};
            minimize sum(x)
            for i=1:K
                a{i}: x(i) >= tau - abs(G(:,i)'*m_l)^2 - 2*real(m_l'*G(:,i)*G(:,i)'*(m-m_l));
                b{i}: x(i) >= 0;
            end
            norm(m) <= 1;
            cvx_end
            
            num_of_users=sum(abs(G'*m).^2>=tau);
            if num_of_users>max_num_of_users
                m_final = m;
                max_num_of_users = num_of_users;
            end
            
            % Stop criteria
            if norm(m-m_l)<m_eps || num_of_users==K
                break;
            else
                m_l = m;
            end
            
        end  
        m = m_final;
        
        
    case 3  % Heuristic method
        max_num_of_users=-1;
        for iter=1:Iter
            y = G'*m_l;
            
            lambda = 1./(theta*abs(y).^2+tau) .* y;  % best
%             lambda = exp(-theta*(abs(y)-sqrt(tau))) .* y;
%             lambda = exp(-theta*(abs(y).^2-tau)) .* y;
%             lambda = 1./(theta*abs(y)+sqrt(tau))./abs(y).* y;
%             
%             index = abs(y)<sqrt(tau);
%             lambda = (abs(y).^2-tau).*(-100*index + 1*(1-index));
            
            m = G*lambda;
            m = m/norm(m);
            
            num_of_users=sum(abs(G'*m).^2>=tau);
            if num_of_users>max_num_of_users
                m_final = m;
                max_num_of_users = num_of_users;
            end
            
            % Stop criteria
            if norm(m-m_l)<m_eps || num_of_users==K
                break;
            else
                m_l = m;
            end
        end
        m = m_final;
        
    case 4  % Random beamforming       
        m = randn(N,1)+randn(N,1)*1i;
        m = m/norm(m);
end

y = G'*m;
active=(abs(y).^2>=1/MSE_threshold );
num_of_users=sum(active);

active_user=(1:K)';
active_user=active_user(active);


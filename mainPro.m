function [] = mainPro(ID)
% main function for Bayesian estimates
%% basic setup
% tot = 3e4; burn = 25e3; 
tot = 6e3; burn = 5e3; 
% tot = 200; burn = 100; 
objrate = 0.44;
reportrate = 1; verbose = 1;

%% decipher job identifier
nChain = 201;  % set number of MCMC runs
ch = str2double(num2str(ID));  ch0 = ch; 
nvar = ceil(ch/nChain); ch = ch - (nvar-1)*nChain;

%% import data, initialize parameters
% [x, E] = readData();  
if ch == 1; load(strcat('Drq1m1p1.mat')); x = x'; E.lbs = E.lbs'; E.ubs = E.ubs'; ch=ch-1; 
else ch=ch-1; load(strcat('Ds',num2str(ch),'q1m1p1.mat')); x = x'; E.lbs = E.lbs'; E.ubs = E.ubs'; 
end
fprintf('model = %d, replication = %d:\n', nvar, ch)
x0 = x;
nx = numel(x); %np = max(E.flags); 
% load('out.mat') % use the initial value from mainproFreq.m
% load(strcat('inits',num2str(nvar),'.mat')) %'init_para','Hinds','updateFlags'

% load(strcat('inits',num2str(nvar),'.mat')) 
% x = init_para(:,3)'; E.lbs = init_para(:,4)'; E.ubs = init_para(:,5)';
np = max(E.flags); %numel(updateFlags); 

%% boundary adjustment, make initial value not hit the boundary
% E.lbs(isinf(E.lbs)) = -30; E.ubs(isinf(E.ubs)) = 30; 
for k = 1:nx; if x(k)==E.lbs(k); x(k)=x(k)+1e-4; end; if x(k)==E.ubs(k); x(k)=x(k)-1e-4; end; end

% x(E.flags==9) = 0.1; 

rng('default'); rng(ch0*8);
% x = x.*unifrnd(.9,1.1,[1,length(x)]) + 1e-4;
loglike = sum(getLogLik(x, E, 1:3, nvar));

%% adaptive MCMC setup
batchLen = 50; 
accepts = zeros(1,nx); batchNum = 0; batchTot = tot/batchLen; rates = zeros(batchTot, nx);
tunings = -5*ones(1,nx); %initial tuning
% accepts = zeros(1,np); batchNum = 0; batchTot = tot/batchLen; rates = zeros(batchTot, np);
% tunings =-3*ones(1,np); %initial tuning
% tunings([4, 6]) = -6; 
% tunings(7:8) = -[6.5, 5]; 
% tunings([2,3]) = -5; 
% tunings(1) = -1; 

%% run MCMC
matpara = zeros(tot-burn, nx); Ls =zeros(1, tot-burn); 
tic
for iter = 1:tot
    if verbose == 1; fprintf('%4d[%3.3f]', [iter,x(E.flags==5)]);  if(~mod(iter,5)); fprintf('\n'); end; end
    %if verbose == 1; fprintf('%4d[%3.3f]', [iter,loglike]);  if(~mod(iter,5)); fprintf('\n'); end; end
    
%     % component updates
%     for k = 1:nx
%         if E.flags(k)~=7
%             x1 = x; x1(k) = invLogit(normrnd(logit(x(k), E.lbs(k), E.ubs(k)), exp(tunings(k))), E.lbs(k), E.ubs(k));
%             loglike1 = sum(getLogLik(x1, E, 1:3, nvar));
%             MH = loglike1 - loglike + log(x1(k) - E.lbs(k)) + log(E.ubs(k) - x1(k)) -  log(x(k) - E.lbs(k)) - log(E.ubs(k) - x(k));
%             u1 = log(rand(1)); if u1 <= MH;  x = x1;  loglike = loglike1; accepts(k) = accepts(k)+1;  end
%         end
%     end
    
    for k0 = 1:np
        if any([1:6,8] == k0) % do not update nu0
            k = find(E.flags == k0);  %updateFlags(k0)
            %         if k0==6 && nvar == 2 %updateFlags(k0)==6
            %             k = Hinds;
            %         end
            x1 = x; x1(k) = invLogit(normrnd(logit(x(k), E.lbs(k), E.ubs(k)), exp(tunings(k0))), E.lbs(k), E.ubs(k));
            loglike1 = sum(getLogLik(x1, E, 1:3, nvar));
            MH = loglike1 - loglike + sum(log(x1(k) - E.lbs(k)) + log(E.ubs(k) - x1(k)) -  log(x(k) - E.lbs(k)) - log(E.ubs(k) - x(k)));
            u1 = log(rand(1));
            if u1 <= MH
                x = x1;  loglike = loglike1; accepts(k0) = accepts(k0)+1;
            end
        end
    end
    
    if iter > burn; matpara(iter-burn,:) = x; Ls(iter-burn) = loglike; end
    if ~mod(iter,batchLen)
        % if reportrate == 1;disp(num2str(accumarray(E.flags', accepts/batchLen,[numel(E.dims),1],@(x) mean(x))', 2)); end
        if reportrate == 1 ; disp(num2str(accepts/batchLen)); end %&& iter > burn
        batchNum = batchNum+1; 
        accepts = accepts/batchLen;
        rates(batchNum,:) = accepts;
        tunings = tunings + sign((accepts>objrate)-0.5).*min(0.01,1/sqrt(batchNum));
        accepts = zeros(1,nx);
    end
end
runtime = toc/3600;
fprintf('\n%d iterations are done with elapsed time %.2f hours.\n', tot, runtime)

plot(matpara(:,E.flags==5))
mean(exp(matpara(:,E.flags==5)))

save(strcat('out',num2str(nvar),'_',num2str(ch),'.mat'),'runtime','matpara','Ls','rates','x0')
end

function [y] = logit(x, a, b)
y = log((x-a)./(b-x)); 
end

function [x] = invLogit(y, a, b)
x = a + (b-a)./(1+exp(-y)); 
end


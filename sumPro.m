function [] = sumPro()
%% summay procedure for mainPro.m output
global E

% nsim = 1e2; 
% for nvar = 1:2
%     for i = 0:nsim
%         load(strcat('out',num2str(nvar),'_',num2str(i),'.mat'))
%         if i == 0; mat = nan(numel(x0),4,1+nsim); end
%         mat(:,:,i+1) = [x0; mean(matpara,1); quantile(matpara,[.025,.975],1)]';
%     end
%     save(strcat('paras',num2str(nvar),'.mat'),'mat')
% end


dirs = './mat/'; 
figdirs = './mat/figures/';
filename = strcat(dirs, 'fig.ps');
% [~, E] = readData();
load(strcat('Drq2m1p1.mat'));
% E.J2 = 2;
JK2=E.J2*E.K2; JK = E.J*E.K0; AK=(E.A-1)*E.K0; KK = E.K*E.K1;
nam = [ strcat(repmat({'\lambda (j='},[1,JK2]),cellfun(@num2str,num2cell(kron(ones(1,E.K2),1:E.J2)),'UniformOutput', false), ...
    repmat({', k='},[1,JK2]), cellfun(@num2str,num2cell(kron(1:E.K2, ones(1,E.J2))),'UniformOutput', false), repmat({')'},[1,JK2])), ...
    'q','\tau','\eta','M','\nu_0', ...
    strcat(repmat({'N_{j1k} (k='},[1,JK]),cellfun(@num2str,num2cell(kron(ones(1,E.J),1:E.K0)),'UniformOutput', false), ...
    repmat({', j='},[1,JK]), cellfun(@num2str,num2cell(kron(1:E.J, ones(1,E.K0))),'UniformOutput', false), repmat({')'},[1,JK])), ...
    strcat(repmat({'N_{1ak} (k='},[1,AK]),cellfun(@num2str,num2cell(kron(ones(1,E.A-1),1:E.K0)),'UniformOutput', false), ...
    repmat({', a='},[1,AK]), cellfun(@num2str,num2cell(kron(1:(E.A-1), ones(1,E.K0))),'UniformOutput', false), repmat({')'},[1,AK])), ...
    strcat(repmat({'\pi (k='},[1,KK]),cellfun(@num2str,num2cell(kron(ones(1,E.K),1:E.K1)),'UniformOutput', false), ...
    repmat({', k1='},[1,KK]), cellfun(@num2str,num2cell(kron(1:E.K, ones(1,E.K1))),'UniformOutput', false), repmat({')'},[1,KK])), ...
    ];

DICs = zeros(2, 4); 
for nvar = 1:2
    k0 = 1;
    chs = 0; %1:3; 
    if nvar == 1; chs = [2,3,4,5]; end
    if nvar == 2; chs = [1,2,6,8]; end
    nch = length(chs);
    niter = 5e3; burn = 0; thin = 5;  %[8,3,6]
    nsample = (niter-burn)/thin;
    iteruse = burn + (1:nsample)*thin; tot = nch*nsample; Xmeans = 0; 
    for ch = 1:nch
        ch0 = chs(ch);
        load(strcat(dirs,'out',num2str(nvar),'_',num2str(ch0),'.mat'))
        [matpara, flags, xmean] = translatePara(matpara(iteruse,:)); % translatePara(matpara(iteruse,:), E);
        if k0 ==1 %for initializing...
            npara = size(matpara,2); MCSamples = nan(nch, npara, nsample); LsAll = nan(nsample, nch);
            matParas = nan(tot, npara);
        end
        MCSamples(ch, :, :) = matpara'; %#ok<*NODEF>
        LsAll(:, ch) = Ls(iteruse);
        matParas((ch-1)*nsample + (1:nsample), :) = matpara;
        Xmeans = Xmeans + xmean; 
        k0 = k0+1;
    end
    R = psrf(MCSamples); %display([R,mean(R(3:end)), std(R(3:end))]); %boxplot(R)
    Xmeans = Xmeans/nch; 
    
    sum(getLogLik(Xmeans, E, 1:3, nvar))
    
    plotit = 0;
    for mycondi = 1:plotit
        %subplot(3,4,1), plot(-2*LsAll); title('Deviance','FontSize',12)
        
        %         k0 = 1;
        %         for i = 1:min(npara, 10)
        %             if i == 8; k0 = k0+1; end
        %             subplot_tight(2,5,i,[0.1,0.04]), plot(reshape(MCSamples(:,k0,:),[nch,nsample])')
        %             title(nam{k0},'FontSize',12)
        %             if nvar==1; print('-painters', '-dpsc2', '-r300', filename); end
        %             k0 = k0 + 1;
        %         end
        %         tightfig;
        
        %         for i = 1:16
        %             subplot(4,4,i), plot(reshape(MCSamples(:,npara-i,:),[nch,nsample])')
        %             title(nam{npara-i},'FontSize',12)
        %         end
        
        % check the traceplot
        nbin = 10;
        for j = 1:ceil(npara/nbin)
            fprintf('%4d', j);  if(~mod(j,20)); fprintf('\n'); end;
            for i = 1:nbin %npara1 %size(matParas, 2)%
                if (j-1)*nbin+i <= npara
                    subplot_tight(2,5,i,[0.1,0.05])
                    k = (j-1)*nbin+i;
                    plot(reshape(MCSamples(:,k,:),[nch,nsample])')
                    set(gca, 'FontSize', 7)
                    title(nam{k},'FontSize',7)
                else
                    subplot_tight(2,5,i,[0.1,0.05]), plot(nan)
                end
                
                %if k0==1; print('-painters', '-dpsc2', '-r30', filename);
                %else print('-painters', '-dpsc2', '-r30', '-append', filename);  end
            end
            
            tightfig;
            saveas(gcf, strcat(figdirs,'mcmc',num2str(nvar),'_',num2str(j)), 'png')
        end
        
    end
    for mycondi = 1:(1-plotit)
        nams1 = {'pure diffusion model','pure overlap model'};
        subplot_tight(1,2,nvar,[0.1,0.08]); plot(-2*LsAll); title(strcat('Deviance for',nams1{nvar}),'FontSize',12); 
        if nvar == 2; tightfig; end
    end
    
    D1 = -2*mean(reshape(LsAll, [1,numel(LsAll)])); 
    D2 = -2*sum(getLogLik(Xmeans, E, 1:3, nvar));
    pD = D1 - D2;
    DIC = 2*D1 - D2;
    disp(num2str([D1, D2, pD, DIC]))
    DICs(nvar,:) = [D1, D2, pD, DIC];
       
    mat = nan(npara,4);
    for i = 1:npara
        mat(i,1) = mean(matParas(:,i)); mat(i,2) = std(matParas(:,i));
        %[lb, ub] = FindHPDset(matParas(:,i), 0.95, []);
        tmpq = quantile(matParas(:,i), [0.025, 0.975]);
        lb = tmpq(1); ub = tmpq(2);
        mat(i,3) = lb(1); mat(i,4) = ub(1);
    end
    mat = [mat, flags', translatePara(E.lbs)', translatePara(E.ubs)', R'];
    
    J = sum(flags==6); 
    J = sqrt(J); 
    Pi = reshape(mat(flags==6,1), J, J); 
    
    % if not converge, save as initial values
    % initvec = mat(1:(size(mat,1)-1), 1)'; inittau2 = mat(size(mat,1), 1);
    % save(strcat('inits_',num2str(nvar),'.mat'),'initvec','inittau2')
    
    % disp(num2str(mat))
    save(strcat(dirs, 'paras',num2str(nvar),'.mat'),'mat','nam')
end

save(strcat(dirs,'DICs.mat'),'DICs')
end

function [y, flags, Xmean] = translatePara(x)
%% translate x into y. CAUTION: differ from projects
global E
Xmean = mean(x,1); % posterior mean for original parameterss
y = x;
for k = [2,4,5,9]; y(:, E.flags==k) = exp(y(:, E.flags==k)); end
flags = E.flags;
for k = [6,7,8]; y(:, flags==k) = []; flags(flags==k)=[]; end

N_j1k = exp(x(:, E.flags==7));
N_1ak = exp(x(:, E.flags==8));
Pi0 = [exp(x(:, E.flags==6)), ones(size(x,1), E.K)];
for k = 1:E.K
    inds = (k-1)+(1:E.K:size(Pi0,2));
    Pi0(:,inds) = Pi0(:,inds)./ repmat(sum(Pi0(:,inds), 2), [1, numel(inds)]);
end
y = [y, N_j1k, N_1ak, Pi0];
flags = [flags, 7*ones(1,sum(E.flags==7)), 8*ones(1,sum(E.flags==8)), 6*ones(1,size(Pi0,2))];
end

function [LBout,UBout] = FindHPDset(Samples,p,npoints)
%%Function to find the 100p % HPD set based on Samples
if isempty(npoints)
    npoints=200;
end

[f,x] = ksdensity(Samples,'npoints',npoints); N = size(Samples,1); maxf = max(f);
step = maxf/npoints;
HPDdensity = (step:step:maxf);
NHPD = size(HPDdensity,2);
LB = cell(1,NHPD); UB = cell(1,NHPD); Prob = zeros(1,NHPD);
for i=1:NHPD
    indices0 = find(HPDdensity(NHPD-i+1) < f);
    if ~isempty(indices0)
        indices1 = find(diff(indices0)> 1);
        if isempty(indices1)
            LB{i} = x(indices0(1)); UB{i} = x(indices0(end));
        elseif (size(indices1,1)==1)
            LB{i} = [x(indices0(1)) x(indices0(indices1(1)+1))];
            UB{i} = [x(indices0(indices1(1))) x(indices0(end))];
        else
            LB{i} = x(indices0(1)); UB{i} = [];
            for j=1:(size(indices1,2)-1)
                LB{i} = [LB{i} x(indices0(indices1(j)+1))];
                UB{i} = [UB{i} x(indices0(indices1(j)))];
            end
            UB{i} =[UB{i} x(indices0(end))];
        end
    end
    Ns = size(LB{i},2);
    count = zeros(1,Ns);
    for j=1:Ns
        count(j) = sum((LB{i}(j) <= Samples).*(Samples <= UB{i}(j)));
    end
    Prob(i) = sum(count)/N;
end
[minval indexmin] = min(abs(Prob - p));
LBout = LB{indexmin};
UBout = UB{indexmin};
end

function [R,neff,V,W,B] = psrf(varargin)
%PSRF Potential Scale Reduction Factor
%
%   [R,NEFF,V,W,B] = PSRF(X) or
%   [R,NEFF,V,W,B] = PSRF(x1,x2,...,xs)
%   returns "Potential Scale Reduction Factor" (PSRF) for
%   collection of MCMC-simulations. X is a NxDxM matrix
%   which contains M MCMC simulations of length N, each with
%   dimension D. MCMC-simulations can be given as separate
%   arguments x1,x2,... which should have the same length.
%
%   Returns
%     R     PSRF in a row vector of length D
%     neff  estimated effective number of samples M*N*V/B
%     V     estimated mixture-of-sequences variances
%     W     estimated within sequence variances
%     B     estimated between sequence variances
%
%   The idea of the PSRF is that if R is not near 1 (below 1.1 for
%   example) one may conclude that the tested samples were not from
%   the same distribution (chain might not have been converged
%   yet).
%
%   If only one simulation is given, the factor is calculated
%   between first and last third of the chain. Note that use of
%   only one chain will produce over-optimistic result.
%
%   Method is from:
%      Brooks, S.P. and Gelman, A. (1998) General methods for
%      monitoring convergence of iterative simulations. Journal of
%      Computational and Graphical Statistics. 7, 434-455. Note that
%      this function returns square-root definiton of R (see Gelman
%      et al (2003), Bayesian Data Analsyis p. 297).
%
%   See also
%     CPSRF, MPSRF, IPSRF

% Copyright (C) 1999 Simo Särkk?
% Copyright (C) 2003 Aki Vehtari
%
% This software is distributed under the GNU General Public
% Licence (version 2 or later); please refer to the file
% Licence.txt, included with the software, for details.

% 2004-01-22 Aki.Vehtari@hut.fi Added neff, R^2->R, and cleaning

% In case of one argument split to two halves (first and last thirds)
onechain=0;
if nargin==1
    X = varargin{1};
    if size(X,3)==1
        n = floor(size(X,1)/3);
        x = zeros([n size(X,2) 2]);
        x(:,:,1) = X(1:n,:);
        x(:,:,2) = X((end-n+1):end,:);
        X = x;
        onechain=1;
    end
elseif nargin==0
    error('Cannot calculate PSRF of scalar');
else
    X = zeros([size(varargin{1}) nargin]);
    for i=1:nargin
        X(:,:,i) = varargin{i};
    end
end

[N,D,M]=size(X);

if N<1
    error('Too few samples');
end

% Calculate means W of the variances
W = zeros(1,D);
for n=1:M
    x = X(:,:,n) - repmat(mean(X(:,:,n)),N,1);
    W = W + sum(x.*x);
end
W = W / ((N-1) * M);

% Calculate variances B (in fact B/n) of the means.
Bpn = zeros(1,D);
m = mean(reshape(mean(X),D,M)');
for n=1:M
    x = mean(X(:,:,n)) - m;
    Bpn = Bpn + x.*x;
end
Bpn = Bpn / (M-1);

% Calculate reduction factors
S = (N-1)/N * W + Bpn;
R = (M+1)/M * S ./ W - (N-1)/M/N;
V = R .* W;
R = sqrt(R);
B = Bpn*N;
neff = min(M*N*V./B,M*N);
if onechain & (nargout>1)
    neff=neff*3/2;
end
end

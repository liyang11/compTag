function [loglik] = getLogLik(x, E, comp, model)
% get negative loglikelihood
%% Translate parameter vector into individual parameters
paras = mat2cell(x, 1, E.dims);

lambda_jk = reshape(paras{1}, [E.J2, E.K2]);

q_jk = exp(paras{2});
if E.parametrizations(1)==1; q_jk = q_jk*ones(E.J,E.K1); end
if E.parametrizations(1)==2; q_jk = repmat(q_jk',[1,E.K1]); end
if E.parametrizations(1)==3; q_jk = repmat(q_jk,[1,E.J]); end
if E.parametrizations(1)==4; q_jk = reshape(q_jk, [J,E.K1]) ; end

tau = paras{3};
eta = exp(paras{4});

Mtmp = exp(paras{5});
M_jak = zeros(E.K1,E.J,E.A);
if E.parametrizations(2)==1; M_jak =Mtmp*ones(E.K1,E.J,E.A); end
if E.parametrizations(2)==2; for j=1:E.J; M_jak(:,j,:) = Mtmp(j); end; end
if E.parametrizations(2)==3; for a=1:E.A; M_jak(:,:,a) = Mtmp(a); end; end
if E.parametrizations(2)==4; for k=1:E.K1; M_jak(k,:,:) = Mtmp(k); end; end

% change 2016-10-6
M_jak(:,1,:) = M_jak(:,1,:)/12; 

d3 = 1; if E.parametrizations(3)==2; d3 = E.A; end; 
H0 = reshape(paras{6}, [E.K, E.K1-1, d3]);
R_kj = exp(reshape(paras{7}, [E.K0, E.J]));
N0_ak = exp(reshape(paras{8}, [E.K0, E.A-1]));
nu0 = paras{9};

Pi0 = zeros([E.K, E.K1, d3]); for a=1:d3; for k=1:E.K; Pi0(k,:,a) = exp([H0(k,:,a),0]);
        Pi0(k,:,a) = Pi0(k,:,a)/sum(Pi0(k,:,a)); end; end

%% Introduce intermediate parameteres
Pi = zeros(E.K,E.K1,E.A,E.J);
if E.parametrizations(3)==1; for a = 1:E.A; for j = 1:E.J; Pi(:,:,a,j) = Pi0; end; end; end
if E.parametrizations(3)==2; for j=1:E.J; Pi(:,:,:,j) = Pi0; end; end
% if model == 2; Pi = Pi(1:E.K0,:,:,:); end

fac = 10^eta*exp(-tau*10);
v_a = E.agevec.^eta.*exp(-tau*E.agevec)/fac;
Fs = zeros(E.K1,E.J,E.A,E.G);
for i = 1:E.ncat; Fs(E.cat_KJG(i,1), E.cat_KJG(i,2),:,E.cat_KJG(i,3)) = ...
        v_a * (q_jk(E.cat_KJG(i,2),E.cat_KJG(i,1)) * E.cat_effort(i)); end
Fs_sumOverG = sum(Fs, 4);
S = exp(-(Fs_sumOverG + M_jak));
nu = nu0*ones(E.G,E.J,E.K);  %this is CV
loglik = zeros(1,4);

%% compute likelihood 1
if any(comp==1)
    %============================ compute L1
    P = zeros(E.I,E.A,E.K,E.G,E.J,E.K1);
    for m0 = 1:E.ntag %only need for cases with released fish, otherwise both R and T_noRec are 0
        i = E.tag_dbIAK(m0,1); a = E.tag_dbIAK(m0,2); k = E.tag_dbIAK(m0,3); 
        for m1 = 1:E.ncat % only need for cases with catching efforts, otherwise P=0
            g = E.cat_KJG(m1,3); j = E.cat_KJG(m1,2); k1 = E.cat_KJG(m1,1);
            if j >= i 
                a1 = min( a+ j-i, E.A );
                i1 =1;
                if k1<=3; i2 = 1; elseif k1<=5; i2 = 2; elseif k1<=8; i2 = 3; elseif k1==9; i2 = 4; else i2 = 5; end
                delta_gjk0 = E.cat_delta_gjk(m1); 
                phi_ij = phiFun( (j - i)*12 + 1 );
                P(i,a,k,g,j,k1) = phi_ij*(delta_gjk0+(1-delta_gjk0)*lambda_jk(i1,i2));
                P(i,a,k,g,j,k1) = P(i,a,k,g,j,k1) * Fs(k1,j,a1,g) *(1-S(k1,j,a1))/(Fs_sumOverG(k1,j,a1)+M_jak(k1,j,a1));
                if j==i
                    P(i,a,k,g,j,k1) = P(i,a,k,g,j,k1)*Pi(k,k1,a1,j);
                else
                    for j0 = 0:(j-i-1)
                        a0 = min(a+j0, E.A);
                        if model == 1
                            if j0==0; Delta_iaj = Pi(1:E.K0,:,a0,i+j0) * diag(S(:,i+j0,a0)) ;
                            else Delta_iaj = Delta_iaj * (Pi(:,:,a0,i+j0) * diag(S(:,i+j0,a0)) ); end
                        elseif model == 2
                            if j0==0; Delta_iaj = diag( Pi(1:E.K0,:,a0,i+j0) *S(:,i+j0,a0) ) ;
                            else Delta_iaj = diag( Delta_iaj * (Pi(1:E.K0,:,a0,i+j0) * S(:,i+j0,a0)) ); end
                        end
                    end
                    if model == 1; Delta_iaj = Delta_iaj * Pi(:,:,min(a0+1, E.A), j);
                    elseif model == 2; Delta_iaj = Delta_iaj * Pi(1:E.K0,:,min(a0+1, E.A), j); end
                    P(i,a,k,g,j,k1) = P(i,a,k,g,j,k1) * Delta_iaj(k,k1);
                end
            end
        end
    end
    tmpP = sum(sum(sum(P,6),5),4);
    loglik(1) = sum(E.rec_R.*log(P(E.ind_rec_i6)))   +   sum(E.tag_noRec.*log(1-tmpP(E.ind_tag_i6))); 
    
    
    %*** old code which is wrong, miss some P values for no-recapture
    % ** we shall not use one single for loop for E.nrec. Rather, we need
    % two loops for E.ntag and E.ncat respectively. See new codes above
    % 
    %     P = E.rec_phi_ij.*(E.rec_delta_gjk+(1-E.rec_delta_gjk).*lambda_jk(E.ind_rec_JK));
    %     %-------------------- compute P using model 1
    %     P = P.*( Fs(E.ind_rec_KJAG).*(1-S(E.ind_rec_KJA))./...
    %         (Fs_sumOverG(E.ind_rec_KJA)+M_jak(E.ind_rec_KJA)) );
    %     for i0 = 1:E.nrec
    %         i = E.rec_IJKK1AG(i0,1); j = E.rec_IJKK1AG(i0,2); k = E.rec_IJKK1AG(i0,3);
    %         k1 = E.rec_IJKK1AG(i0,4); a = E.rec_IJKK1AG(i0,5);
    %         a1 = min( a+ E.rec_lags(i0), E.A );
    %
    %         if E.rec_lags(i0)==0;
    %             P(i0) = P(i0)*Pi(k,k1,a1,j);
    %         else
    %             for j0 = 0:(E.rec_lags(i0)-1)
    %                 a0 = min(a+j0, E.A);
    %                 if model == 1
    %                     if j0==0; Delta_iaj = Pi(1:E.K0,:,a0,i+j0) * diag(S(:,i+j0,a0)) ;
    %                     else Delta_iaj = Delta_iaj * (Pi(:,:,a0,i+j0) * diag(S(:,i+j0,a0)) ); end
    %                 elseif model == 2
    %                     if j0==0; Delta_iaj = diag( Pi(1:E.K0,:,a0,i+j0) *S(:,i+j0,a0) ) ;
    %                     else Delta_iaj = diag( Delta_iaj * (Pi(1:E.K0,:,a0,i+j0) * S(:,i+j0,a0)) ); end
    %                 end
    %             end
    %             if model == 1; Delta_iaj = Delta_iaj * Pi(:,:,min(a0+1, E.A), j);
    %             elseif model == 2; Delta_iaj = Delta_iaj * Pi(1:E.K0,:,min(a0+1, E.A), j); end
    %             P(i0) = P(i0) * Delta_iaj(k,k1);
    %         end
    %     end
    %     tmpP = [accumarray(E.rec_iak, P); 0];
    %     loglik(1) = sum(E.rec_R.*log(P)) + sum( E.tag_noRec.*log(1-tmpP(E.ind_rec2tag)) );
    
end

%% compute likelihood 2
if any(comp==2)
    lambda = lambda_jk(E.ind_cat_JK)';
    delta = E.cat_delta_gjk(E.ind_cat_positive, :);
    loglik(2) = sum(    E.cat_R_um(E.ind_cat_positive).*log(lambda) -  ...
        E.cat_R_tot(E.ind_cat_positive).*log( delta + (1-delta).*lambda )   );
end

%% compute likelihood 3
if any(comp==3)
    %============================ compute L_CT
    N = nan(E.K0, E.J, E.A);
    for j = 1:E.J; for k = 1:E.K0; N(k,j,1) = R_kj(k,j); end; end
    for a = 2:E.A; for k = 1:E.K0; N(k,1,a) = N0_ak(k,a-1); end; end
    for a = 2:E.A; for j = 2:E.J; for k = 1:E.K0;  N(k,j,a) = sum(N(k,j-1,a-1)*Pi(k,:,a-1,j-1)*S(:,j-1,a-1)); end; end; end
    EC = zeros(E.K1, E.J, E.G, E.A);
    
    %useinds = find(E.cat(:,3) <= E.K0)'; 
    l_ct = 0.0;
    for m0 = 1:E.ncat %useinds
        j = E.cat_KJG(m0,2); k1 = E.cat_KJG(m0,1); g = E.cat_KJG(m0,3);
        EC(k1,j,g,:) = squeeze( Fs(k1,j,:,g).*(1-S(k1,j,:))./(Fs_sumOverG(k1,j,:)+M_jak(k1,j,:)) .*...
            sum(N(:,j,:).*Pi(1:E.K0,k1,:,j),1) );
        mu = log(sum(EC(k1,j,g,:)));
        
        %VC = (nu(g,j,k1)*mu)^2;
        VC = log(nu(g,j,k1)^2+1);
        
        l_ct = l_ct - 0.5*(log(VC) + (E.log_cat_catch(m0)-mu)^2/VC);
    end
    loglik(3) = l_ct;
    
    sumEC = sum(EC, 4);
    tmp = EC(E.ind_age_KJGA)./sumEC(E.ind_age_KJG); 
    inds = find(tmp>0); 
    tmp = E.age_n_gjka(inds).*log(tmp(inds)); 
    loglik(4) = sum(tmp); %l_cp %(E.induse_age)
end

%%
% loglik = sum(loglik);
end

function [y] = phiFun(x)
y = 1-0.31./(1+exp(8.52-double(x))); 
end
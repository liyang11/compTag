function [] = mainProFreqHPC(ID)

repl = str2double(num2str(ID));

modID = 2;  
load(strcat('Ds',num2str(repl),'q1m1p1.mat'));
A = []; b = []; Aeq = []; beq = []; nonlcon = [];
lbs = E.lbs;
ubs = E.ubs;
options = optimoptions('fmincon','Display','iter','Algorithm','sqp','MaxFunEvals',5e5,'MaxIter',2e4);%'Display','iter',
x0 = x';
num2str(getLogLik(x0, E, 1:4, 1))

% fixInds = find(E.flags~=5&E.flags~=1&E.flags~=2);
fixInds = find(E.flags==9);
freeInds = 1:numel(x0); freeInds(fixInds) = [];

x00 = x0(freeInds);
objfunc = @(x) - sum(getLogLikWrapper(x, freeInds, fixInds, x0(fixInds), E, 1:3, modID));
objfunc(x00)

tic
[x1,fval,~,~,~,~,~] = fmincon(objfunc, x00, A,b,Aeq,beq,lbs(freeInds), ubs(freeInds), nonlcon,options);
runtime = toc/60;
fprintf('\nElapsed time %.2f minutes.\n', runtime)

disp(fval)

xtmp = x0; xtmp(freeInds) = x1;
Mest = xtmp(E.flags==5);
disp(exp(Mest))

num2str(getLogLik(xtmp, E, 1:4, 1))

a = [x,xtmp',E.lbs,E.ubs]; a(freeInds,:)

save(strcat('out',num2str(repl),'.mat'),'Mest','fval','xtmp')

end

function [] = comb()
nsim = 200; 
for i = 1:nsim
    load(strcat('out',num2str(i),'.mat'))
    if i== 1; Ms = zeros(nsim, numel(Mest)); end
    Ms(i,:) = Mest;
end
quantile(Ms,[0.025,0.5,0.975],1)
end
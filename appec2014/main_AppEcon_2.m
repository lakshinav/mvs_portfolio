%% Set options for optimizer
opt  =  optimset('fmincon');
opt=optimset(opt, 'Display', 'final');
opt=optimset(opt, 'TolX', 1e-6);
opt=optimset(opt, 'TolFun', 1e-5);
opt=optimset(opt, 'TolCon', 1e-6);
opt=optimset(opt, 'MaxFunEvals', 8e4);
opt=optimset(opt, 'MaxIter', 1e5);
% opt=optimset(opt, 'FinDiffType', 'central'); % 'central' more accurate but more time-consuming
opt=optimset(opt, 'Algorithm', 'interior-point'); %interior-point sqp active-set trust-region-reflective
% opt=optimset(opt, 'InitTrustRegionRadius',0.01);
opt=optimset(opt, 'AlwaysHonorConstraints', 'none');

%% Create structure to store the results
vc=struct('x0',[],'x',[],'fv',[],'ef',[],'et',[],'op',[],'grad',[],'hess',[],...
    'matr',[],'stderr',[],'signif',[],'mes',[]); 

%% Estimate GO-GARCH
warning off all
%{ 
  Types:
  for GO-GARCH:
  2 - homogeneous with spatial effects (nparams=2*n+3)
  3 - group-homogeneous with spatial effects (nparams=2*n+3*5)
  4 - heterogeneous with spatial effects (nparams=5*n)
  5 - original (nparams=n*(n-1)/2+3*n)
  6 - scalar (nparams=n*(n-1)/2+2+n)
%}
type=2:6;
fileName='estGoGarch2';
vc(1,1).x0=[0.1*ones(1,n) 0.5*ones(1,n) 0.1*ones(1,n) 0.02 0.01 0.03];
vc(2,1).x0=[0.1*ones(1,n) 0.5*ones(1,n) 0.1*ones(1,n) 0.02*ones(1,5) 0.01*ones(1,5) 0.03*ones(1,5)];
vc(3,1).x0=[0.1*ones(1,n) 0.5*ones(1,n) 0.1*ones(1,n) 0.02*ones(1,n) 0.01*ones(1,n) 0.05*ones(1,n)];
vc(4,1).x0=[0.05*ones(1,n*(n-1)/2) 0.5*ones(1,3*n)];
vc(5,1).x0=[0.05*ones(1,n*(n-1)/2) 0.5*ones(1,2+n)];

lb{1}=[zeros(1,3*n) -Inf*ones(1,3)];
lb{2}=[zeros(1,3*n) -Inf*ones(1,3*5)];
lb{3}=[zeros(1,3*n) -Inf*ones(1,3*n)];
lb{4}=-10*ones(1,n*(n-1)/2+3*n);
lb{5}=-10*ones(1,n*(n-1)/2+2+n);
for k=5
    f=@(u) cmptLLspGoGarch(u,y,W,type(k));
    b=@(v) checkSigmaPd(v,y,W,type(k));
tic;
try
[vc(k,1).x,vc(k,1).fv,vc(k,1).ef,vc(k,1).op,~,vc(k,1).grad,vc(k,1).hess]=...
    fmincon(f,vc(k,1).x0,[],[],[],[],lb{k},[],b,opt); % strerr computation
vc(k,1).matr=(vc(k,1).hess^(-1))*(vc(k,1).grad*vc(k,1).grad')*(vc(k,1).hess^(-1)); % covariation matrix
vc(k,1).stderr=sqrt(diag(vc(k,1).matr)'/T); % asymptotic s. e.
vc(k,1).signif=(1-normcdf(abs(vc(k,1).x./vc(k,1).stderr),0,1))*2;
catch err
    e{k}=err;
    disp('error');
    vc(k).fv=NaN;
end
vc(k,1).et=toc;
save(fileName);
disp(k);
end
clear k type f b lb err;
save(fileName);

%% Estimate full GO-GARCH by MFE toolbox
opt=optimset(optimset('fmincon'), 'MaxFunEvals', 5e5);
% [PARAMETERS,LL,HT,VCV,SCORES] = gogarch(DATA,P,Q,GJRTYPE,TYPE,STARTINGVALS,OPTIONS)
i=4;
tic;
[vc(i).x,vc(i).fv,~,vc(i).matr] = gogarch(y',1,1,[],[],[],opt);
vc(i).et=toc;
clear i;

%% Check if Sigma p.d. in GO_GARCH
[~,sigma]=cmptLLspGoGarch(vc(3).x,y,W,2);
p=-1*ones(1,T);
for i=1:T
    [~,p(i)]=chol(sigma(:,:,i));
end
if all(not(p))
    disp('Sigma is p.d.');
else
    disp('NOT');
end
clear i p;
plot(squeeze(sigma(5,5,50:end)));
clear sigma;

%% Estimate CCC
warning off all
%{ 
  Types:
  for CCC:
  2 - homogeneous with spatial effects (nparams=4*n+3)
  3 - group-homogeneous with spatial effects (nparams=4*n+3*5)
  4 - heterogeneous with spatial effects (nparams=7*n)
  5 - original
  6 - diagonal
  7 - scalar
%}
type=[ones(1,6) 4];
fileName='estCcc';
% vc(1,1).x0=[2e-5*ones(1,n) 0.1*ones(1,n) 0.1*ones(1,n) 0.1*ones(1,n) 0.02 0.01 0.03];
% vc(2,1).x0=[2e-5*ones(1,n) 0.1*ones(1,n) 0.1*ones(1,n) 0.1*ones(1,n) 0.02*ones(1,5) 0.01*ones(1,5) 0.03*ones(1,5)];
vc(7,1).x0=[2e-5*ones(1,n) 0.1*ones(1,n) 0.1*ones(1,n) 0.1*ones(1,n) 0.02*ones(1,n) 0.01*ones(1,n) 0.03*ones(1,n)];
% vc(6,1).x0=0.02*ones(1,2*n*n+n+n*(n+1)/2);
% vc(7,1).x0=0.02*ones(1,2*n+n+n*(n+1)/2);
% vc(8,1).x0=0.02*ones(1,2+n+n*(n+1)/2);
for k=7
    f=@(u) cmptLLspCcc(u,y,W,type(k));
    b=@(v) checkSigmaPd(v,W,type(k));
tic;
try
[vc(k,1).x,vc(k,1).fv,vc(k,1).ef,vc(k,1).op,~,vc(k,1).grad,vc(k,1).hess]=...
    fmincon(f,vc(k,1).x0,[],[],[],[],setBndsCcc(type(k),n),[],b,opt); % strerr computation
vc(k,1).matr=(vc(k,1).hess^(-1))*(vc(k,1).grad*vc(k,1).grad')*(vc(k,1).hess^(-1)); % covariation matrix
vc(k,1).stderr=sqrt(diag(vc(k,1).matr)'/T); % asymptotic s. e.
vc(k,1).signif=(1-normcdf(abs(vc(k,1).x./vc(k,1).stderr),0,1))*2;
catch err
    e{k}=err;
    disp('error');
    vc(k).fv=NaN;
end
vc(k,1).et=toc;
save(fileName);
disp(k);
end
clear k type f b err;
save(fileName);

%% Estimate full CCC by MFE toolbox
opt=optimset(optimset('fmincon'), 'MaxFunEvals', 5e5);
% [PARAMETERS,LL,HT,VCV,SCORES] = ccc_mvgarch(DATA,DATAASYM,P,O,Q,GJRTYPE,STARTINGVALS,OPTIONS)
i=10;
tic;
[vc(i).x,vc(i).fv,~,vc(i).matr] = ccc_mvgarch(y',[],1,0,1,[],[],opt);
vc(i).et=toc;
clear i opt;
save(fileName);

%% Calculate BIC
i=1:6;
vc=vcGg;
for j=i
    l(j)=length(vc(j).x);
end
[~,b]=aicbic(-[vc(i).fv],l(l>0),T); 
b=b';
clear j i;

%% Calculate BIC 2
nam={'Bekk', 'Gg', 'Ccc'};
s=[6 3];
fv=zeros(s);
l=zeros(s);
a=zeros(s);
b=zeros(s);
for i=1:3
    tmp=eval(['vc' char(nam{i})]);
    fv(:,i)=-[tmp.fv]';
    for j=1:6
        l(j,i)=length(tmp(j).x);
    end
    [a(:,i),b(:,i)]=aicbic(fv(:,i),l(:,i),T*ones(s(1),1));
end
clear i j tmp nam s;

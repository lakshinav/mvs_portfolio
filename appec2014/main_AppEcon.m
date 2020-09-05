%% Download
% 01/01/2011 - 22/02/2014
% 5 industries: oil and gas; banks; utilities; industrial machinery; IT. From SnP500
names={'cvx','xom','chk','mur',...
    'c','jpm','pbct','sti',...
    'cms','dte','fe','nrg',...
    'dhr','dov','itw','joy',...
    'csco','symc','orcl','ads'};
n=20;
T=789;
%% Make cell array (first element - last observation)
d=cell(n,1);
for i=1:n
    d{i}=eval(['d' num2str(i)]);
end
clear i;

%% Check dates
date=zeros(T,n);
for i=1:n
    date(:,i)=d{i}(:,1);
end
clear i;
disp(sum(sum(date,2)./(date(:,1)*n)~=1));
% It's ok!

%% Make panel
y=zeros(T,n);
for i=1:n
    y(:,i)=d{i}(:,2);
end
clear i;
y=log(y(1:end-1,:)./y(2:end,:)); % log returns
y=y(end:-1:1,:); % first element - first observation
y=y'; % panel n x T
T=T-1;
yWithCondMean=y;

%% Estimate and subtract conditional mean
f=struct('arma',[],'vc',[],'ll',[],'signif',[],'sum',[]);
sp=[[0 0]; [1 0]; [0 1]; [1 1]; [1 2]; [2 1]; [2 2]];
bmin=zeros(1,n);
for i=1:n
    for j=1:7
        [f(i,j).arma,f(i,j).vc,f(i,j).ll,~,~,f(i,j).sum] = garchfit(garchset('R',sp(j,1),'M',sp(j,2),'VarianceModel','Constant'),...
            yWithCondMean(i,:));
    end
    [~,bic]=aicbic([f(i,:).ll],[1 2 2 3 4 4 5],T*ones(1,7));
    bmin(i)=find(bic==min(bic));
    switch bmin(i)
        case 1
            f(i,1).signif=[f(i,bmin(i)).arma.C/f(i,bmin(i)).vc.C,...
                f(i,bmin(i)).arma.K/f(i,bmin(i)).vc.K];
        case 2
            f(i,1).signif=[f(i,bmin(i)).arma.C/f(i,bmin(i)).vc.C,...
                f(i,bmin(i)).arma.K/f(i,bmin(i)).vc.K,...
                f(i,bmin(i)).arma.AR./f(i,bmin(i)).vc.AR];
        case 3
            f(i,1).signif=[f(i,bmin(i)).arma.C/f(i,bmin(i)).vc.C,...
                f(i,bmin(i)).arma.K/f(i,bmin(i)).vc.K,...
                f(i,bmin(i)).arma.MA./f(i,bmin(i)).vc.MA];
        otherwise
            f(i,1).signif=[f(i,bmin(i)).arma.C/f(i,bmin(i)).vc.C,...
                f(i,bmin(i)).arma.K/f(i,bmin(i)).vc.K,...
                f(i,bmin(i)).arma.AR./f(i,bmin(i)).vc.AR,...
                f(i,bmin(i)).arma.MA./f(i,bmin(i)).vc.MA];
    end
    [~,~,cm]=garchsim(f(i,bmin(i)).arma,T,1);
    y(i,:)=yWithCondMean(i,:)-cm';
end
clear i j sp cm bmin f bic;

%% Weight matrix
cap=[ones(1,4) 2*ones(1,4) 3*ones(1,4) 4*ones(1,4) 5*ones(1,4)];
W=zeros(n,n);
for i=1:n
    W(i,:)=cap'==cap(i);
    W(i,i)=0;
    W(i,:)=W(i,:)/sum(W(i,:));
end
clear i cap;

%% Set options for optimizer
opt  =  optimset('fmincon');
opt=optimset(opt, 'Display', 'final');
% opt=optimset(opt, 'TolX', 1e-12);
% opt=optimset(opt, 'TolFun', 1e-5);
% opt=optimset(opt, 'TolCon', 1e-8);
opt=optimset(opt, 'MaxFunEvals', 5e5);
opt=optimset(opt, 'MaxIter', 1e5);
% opt=optimset(opt, 'FinDiffType', 'central'); % 'central' more accurate but more time-consuming
opt=optimset(opt, 'Algorithm', 'interior-point'); %interior-point sqp active-set trust-region-reflective
opt=optimset(opt, 'InitTrustRegionRadius',0.05);
opt=optimset(opt, 'AlwaysHonorConstraints', 'none');

%% Create structure to store the results
vc=struct('x0',[],'x',[],'fv',[],'ef',[],'et',[],'op',[],'grad',[],'hess',[],...
    'matr',[],'stderr',[],'signif',[],'mes',[]);

%% Calculate BIC
i=[1:3 6:8];
for j=i
    l(j)=length(vc(j).x);
end
[~,b]=aicbic(-[vc(i).fv],l(l>0),T); 
b=b';
clear j i;

%% Estimate BEKK
warning off all
%{
  Types:
  for BEKK:
  2 - homogeneous with spatial effects (nparams=3*n+3)
  3 - group-homogeneous (W) with spatial effects (nparams=3*n+3*5)
  4 - heterogeneous with spatial effects (nparams=6*n)
  9 - unrestricted targeted (nparams=n(n+1)/2+2*n^2)
  10 - scalar targeted (nparams=2)
  11 - diagonal targeted (nparams=2*n)
%}
type=[2:4 9 10 11];
fileName='estBekk';
% vc(1,1).x0=[2e-5*ones(1,n) 0.2*ones(1,n) 0.95*ones(1,n) 0.02 0.03 0.02];
% vc(2,1).x0=[2e-5*ones(1,n) 0.2*ones(1,n) 0.95*ones(1,n) 0.02*ones(1,5) 0.03*ones(1,5) 0.02*ones(1,5)];
% vc(3,1).x0=[2e-5*ones(1,n) 0.2*ones(1,n) 0.95*ones(1,n) 0.02*ones(1,n) 0.03*ones(4444444444444444444444444444441,n) 0.02*ones(1,n)];
% vc(4,1).x0=0.02*ones(1,2*n*n+n*(n+1)/2);
% vc(6,1).x0=0.02*ones(1,2*n*n);
vc(7,1).x0=0.2*ones(1,2);
vc(8,1).x0=0.2*ones(1,2*n);
for k=7:8
    f=@(u) cmptLL(u,y,W,type(k-2),1e-3);
%     b=@(v) checkSinvert(v,W,type(k-2));
    tic;
    try
        [vc(k,1).x,vc(k,1).fv,vc(k,1).ef,vc(k,1).op,~,vc(k,1).grad,vc(k,1).hess]=...
            fmincon(f,vc(k,1).x0,[],[],[],[],setBnds(type(k-2),n),[],[],opt); % strerr computation
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

%% Estimate scalar, diagonal and full BEKK by MFE toolbox
opt=optimset(optimset('fmincon'), 'MaxFunEvals', 9e5,'MaxIter',5e3);
% [PARAMETERS,LL,HT,VCV,SCORES] = bekk(DATA,DATAASYM,P,O,Q,TYPE,STARTINGVALS,OPTIONS)
type{1}='Diagonal';
type{2}='Full';
for i=5:6
    tic;
    [vc(i).x,vc(i).fv,~,vc(i).matr] = bekk(y',[],1,0,1,type{i-4},[],opt);
    vc(i).et=toc;
    vc(i).mes=[type{i-4} '_bekk'];
    save(fileName);
end
clear i opt type;
save(fileName);

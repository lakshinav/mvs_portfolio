%% Set options for optimizer
options  =  optimset('fmincon');
options=optimset(options, 'Display', 'final');
options=optimset(options, 'TolX', 1e-12);
% options=optimset(options3, 'TolFun', 1e-5);
options=optimset(options, 'TolCon', 1e-8);
options=optimset(options, 'MaxFunEvals', 5e5);
options=optimset(options, 'MaxIter', 3e4);
% options=optimset(options, 'FinDiffType', 'central'); % 'central' more accurate but more time-consuming
options=optimset(options, 'Algorithm', 'interior-point'); %interior-point sqp active-set trust-region-reflective
% options=optimset(options, 'InitTrustRegionRadius',0.05);
% options=optimset(options, 'AlwaysHonorConstraints', 'none');

%% Create structure for result storage
vc=struct('x0',[],'x',[],'fv',[],'ef',[],'et',[],'op',[],'grad',[],'hess',[],...
    'matr',[],'stderr',[],'signif',[],'mes',[]); 

%% Estimation (spatial GO-GARCH)
warning off all
type=5;
f=@(u) cmptLLspGoGarch(u,y,W,type);
% lb=[zeros(1,2*n) -Inf*ones(1,3*4)];
lb=zeros(1,n*(n-1)/2+2*n);
% b=@(v) checkSigmaPd(v,W,type);
k=3;
% vc(k,1).x0=[0.5*ones(1,n) 0.1*o nes(1,n) 0.02 0.01 0.03];
% vc(k,1).x0=[0.05*ones(1,n) 0.1*ones(1,n) 0.02*ones(1,n) 0.01*ones(1,n) 0.05*ones(1,n)];
% vc(k,1).x0=[0.5*ones(1,n) 0.1*ones(1,n) 0.02*ones(1,4) 0.01*ones(1,4) 0.05*ones(1,4)];
vc(k,1).x0=[0.5*ones(1,n*(n-1)/2) 0.05*ones(1,n) 0.05*ones(1,n)];
tic;
[vc(k,1).x,vc(k,1).fv,vc(k,1).ef,vc(k,1).op]=fmincon(f,vc(k,1).x0,[],[],[],[],lb,[],[],options);
vc(k,1).fv(1,1)=-vc(k,1).fv(1,1);
[vc(k,1).fv(2,1),vc(k,1).fv(3,1)]=aicbic(vc(k,1).fv(1,1),length(vc(k,1).x),T);
disp(vc(k,1).fv); % manual grid-search

% [vc(k,1).x,vc(k,1).fv,vc(k,1).ef,vc(k,1).op,~,vc(k,1).grad,vc(k,1).hess]=...
%     fmincon(f,vc(k-1,1).x,[],[],[],[],lb,[],b,options); % strerr computation
% vc(k,1).fv(1,1)=-vc(k,1).fv(1,1);
% [vc(k,1).fv(2,1),vc(k,1).fv(3,1)]=aicbic(vc(k,1).fv(1,1),length(vc(k,1).x),T);
% vc(k,1).matr=(vc(k,1).hess^(-1))*(vc(k,1).grad*vc(k,1).grad')*(vc(k,1).hess^(-1)); % covariation matrix
% vc(k,1).stderr=sqrt(diag(vc(k,1).matr)'/T); % asymptotic s. e.
% vc(k,1).signif=(1-normcdf(abs(vc(k,1).x./vc(k,1).stderr),0,1))*2;

vc(k,1).et=toc;
beep
clear f b;
disp(vc(k).signif);
[~,sigma]=cmptLLspGoGarch(vc(k).x,y,W,type);
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
plot(squeeze(sigma(5,5,:)));
% clear sigma;

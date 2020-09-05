function [LL,sigma]=cmptLLspGoGarch(par,y,W,type)
%{
Likelihood function for a spatial GO-GARCH model (see Caporin, Paruolo, 2009)
par - parameter vector
y - data
W - weight matrix
%}
n=size(y,1); % number of assets in portfolio
sigma=zeros([n,n,length(y)]);

switch type % all specs are targeted!
    case 1 % without spatial effects
        disp('No sense');
    case 2 % homogeneous with spatial effects (nparams 3*n+3)
        c=par(1:n)';
        a=diag(par(n+1:2*n))+diag(par(3*n+2)*ones(1,n))*W;
        b=diag(par(2*n+1:3*n))+diag(par(3*n+3)*ones(1,n))*W;
        s=eye(n)-diag(par(3*n+1)*ones(1,n))*W;
    case 3 % group-homogeneous with spatial effects (nparams=3*n+3*5)
%         ord=[4	2	3	1	2	3	3	1	2	4	2	1	4	1];
        ord=[ones(1,4) 2*ones(1,4) 3*ones(1,4) 4*ones(1,4) 5*ones(1,4)];
        tmp1=[par(3*n+ord); par(3*n+5+ord); par(3*n+10+ord)];
        c=par(1:n)';
        a=diag(par(n+1:2*n))+diag(tmp1(2,:))*W;
        b=diag(par(2*n+1:3*n))+diag(tmp1(3,:))*W;
        s=eye(n)-diag(tmp1(1,:))*W;
        clear i tmp1 ord;
    case 4 % heterogeneous with spatial effects (nparams=6*n)
        c=par(1:n)';
        a=diag(par(n+1:2*n))+diag(par(4*n+1:5*n))*W;
        b=diag(par(2*n+1:3*n))+diag(par(5*n+1:6*n))*W;
        s=eye(n)-diag(par(3*n+1:4*n))*W;
    case 5 % original go-garch (nparams=n*(n-1)/2+3*n)
        var_covar=cov(y');
        [P, Lam]=eig(var_covar);
        u=umat_(par(1:n*(n-1)/2),n);
        x=P*sqrtm(Lam)*u;
        s=eye(n)/x;
        a=diag(par(n*(n-1)/2+1:n*(n-1)/2+n));
        b=diag(par(n*(n-1)/2+n+1:n*(n-1)/2+2*n));
        c=par(n*(n-1)/2+2*n+1:n*(n-1)/2+3*n)';
    case 6 % scalar (nparams=n*(n-1)/2+n+2)
         var_covar=cov(y');
        [P, Lam]=eig(var_covar);
        u=umat_(par(1:n*(n-1)/2),n);
        x=P*sqrtm(Lam)*u;
        s=eye(n)/x;
        a=par(n*(n-1)/2+1);
        b=par(n*(n-1)/2+2);
        c=par(n*(n-1)/2+3:n*(n-1)/2+2+n)';
    otherwise
        disp('Unexpected type of model.');
end

sigma(:,:,1)=eye(n)/(s'*s);
ll=zeros(size(y,2),1);
ll(1)=-0.5*log(det(sigma(:,:,1)))-0.5*y(:,1)'/sigma(:,:,1)*y(:,1);
v1=ones(n,1);
for j=2:length(ll)
%     v2=(eye(n)-a-b)*ones(n,1)+a*(y(:,j-1).*y(:,j-1))+b*v1;
    v2=c+a*(y(:,j-1).*y(:,j-1))+b*v1; % for spec without targeting
    if v1>0
    sigma(:,:,j)=(eye(n)/s)*diag(v2)*(eye(n)/s)';
    v1=v2;
    ll(j)=-0.5*log(det(sigma(:,:,j)))-0.5*y(:,j)'/sigma(:,:,j)*y(:,j); % likelihood function
    else
            ll(j)=Inf;
    end
end
LL=-sum(ll);
end
function [LL,sigma]=cmptLLspCcc(par,y,W,type)
%{
Likelihood function for a spatial CCC model (see Caporin, Paruolo, 2009)
par - parameter vector
y - data
W - weight matrix
%}
n=size(y,1); % number of assets in portfolio
sigma=zeros([n,n,length(y)]);

switch type % all specs are targeted!
    case 1 % without spatial effects
        disp('No sense');
    case 2 % homogeneous with spatial effects (nparams=4*n+3)
        v=diag(par(1:n));
        a=diag(par(n+1:2*n))+diag(par(3*n+2)*ones(1,n))*W;
        b=diag(par(2*n+1:3*n))+diag(par(3*n+3)*ones(1,n))*W;
        s=eye(n)-diag(par(3*n+1)*ones(1,n))*W;
        c=par(3*n+4:4*n+3)';
        r=(eye(n)/s)*v*(eye(n)/s)';
    case 3 % group-homogeneous with spatial effects (nparams=4*n+3*5)
        v=diag(par(1:n));
%         ord=[4	2	3	1	2	3	3	1	2	4	2	1	4	1];
        ord=[ones(1,4) 2*ones(1,4) 3*ones(1,4) 4*ones(1,4) 5*ones(1,4)];
        tmp1=[par(3*n+ord); par(3*n+5+ord); par(3*n+10+ord)];
        a=diag(par(n+1:2*n))+diag(tmp1(2,:))*W;
        b=diag(par(2*n+1:3*n))+diag(tmp1(3,:))*W;
        s=eye(n)-diag(tmp1(1,:))*W;
        c=par(3*n+16:4*n+15)';
        r=(eye(n)/s)*v*(eye(n)/s)';
        clear i tmp1 ord;
    case 4 % heterogeneous with spatial effects (nparams=7*n)
        v=diag(par(1:n));
        a=diag(par(n+1:2*n))+diag(par(4*n+1:5*n))*W;
        b=diag(par(2*n+1:3*n))+diag(par(5*n+1:6*n))*W;
        s=eye(n)-diag(par(3*n+1:4*n))*W;
        c=par(6*n+1:7*n)';
        r=(eye(n)/s)*v*(eye(n)/s)';
    case 5 % full (nparams=2*n*n+n+n*(n+1)/2)
        a=reshape(par(1:n*n),n,n);
        b=reshape(par(n*n+1:2*n*n),n,n);
        c=par(2*n*n+1:2*n*n+n)';
        r=lowerMatr(par(2*n*n+n+1:end),n);
        r=r*r';
    case 6 % diagonal (nparams=3*n+n*(n+1)/2)
        a=diag(par(1:n));
        b=diag(par(n+1:2*n));
        c=par(2*n+1:3*n)';
        r=lowerMatr(par(3*n+1:end),n);
        r=r*r';
    case 7 % scalar (nparams=n+2+n*(n+1)/2)
        a=par(1);
        b=par(2);
        c=par(3:2+n)';
        r=lowerMatr(par(3+n:end),n);
        r=r*r';
    otherwise
        disp('Unexpected type of model.');
end
sigma(:,:,1)=diag(sqrt(c))*(diag((diag(r)).^(-1/2))*r*diag((diag(r)).^(-1/2)))*diag(sqrt(c));
ll=zeros(size(y,2),1);
ll(1)=-0.5*log(det(sigma(:,:,1)))-0.5*y(:,1)'/sigma(:,:,1)*y(:,1);
v1=c; 
for j=2:length(ll)
    v2=c+a*(y(:,j-1).*y(:,j-1))+b*v1;
    if v1>0
    sigma(:,:,j)=diag(sqrt(v2))*(diag((diag(r)).^(-1/2))*r*diag((diag(r)).^(-1/2)))*diag(sqrt(v2));
    v1=v2;
    ll(j)=-0.5*log(det(sigma(:,:,j)))-0.5*y(:,j)'/sigma(:,:,j)*y(:,j); % likelihood function
    else
            ll(j)=Inf;
    end
end
LL=-sum(ll);
end
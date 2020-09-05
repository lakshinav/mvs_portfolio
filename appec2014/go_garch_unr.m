function [ LL ] = go_garch_unr( b, y, cd )
% likelihood function for a GO-GARCH_unrestricted model
% See Weide (2002)
ll=zeros(size(y,1),1);
var_covar=cov(y);
[P, Lam]=eig(var_covar);
umat=umat_(b(1:cd*(cd-1)/2), cd);
xmat=P*sqrtm(Lam)*umat;
amat=diag(b(cd*(cd-1)/2+1:length(b)-cd));
bmat=diag(b(length(b)-cd+1:length(b)));
cmat=(eye(cd)-amat-bmat)*ones(cd,1);
v_1=cmat;
ll(1)=-0.5*log(det(xmat*diag(v_1)*xmat'))-0.5*y(1,:)/(xmat*diag(v_1)*xmat')*y(1,:)';
for j=2:length(ll)   
   v_t=cmat+amat*(y(j,:).*y(j,:))'+bmat*v_1;
   v_1=v_t;
   ll(j)=-0.5*log(det(xmat*diag(v_t)*xmat'))-0.5*y(j,:)*(xmat*diag(v_t)*xmat')*y(j,:)'; 
end
LL=-sum(ll);
end
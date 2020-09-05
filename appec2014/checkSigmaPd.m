function [c,ceq] = checkSigmaPd(par,y,W,type)
n=size(W,1);
switch type
    case 2 % spatial GO-GARCH / CCC homogeneous
        s_=diag(par(3*n+1)*ones(1,n))*W;
    case 3 % spatial GO-GARCH / CCC group-homogeneous
%         ord=[4	2	3	1	2	3	3	1	2	4	2	1	4	1];
        ord=[ones(1,4) 2*ones(1,4) 3*ones(1,4) 4*ones(1,4) 5*ones(1,4)];
        s_=diag(par(3*n+ord))*W;
        clear ord i tmp1;
    case 4 % spatial GO-GARCH / CCC heterogeneous
        s_=diag(par(3*n+1:4*n))*W;
    case 5 % original
%         s_=0.5*eye(n); % dummy
        var_covar=cov(y');
        [P, Lam]=eig(var_covar);
        u=umat_(par(1:n*(n-1)/2),n);
        s_=P*sqrtm(Lam)*u;
    case 6 % diagonal
%         s_=0.5*eye(n); % dummy
        var_covar=cov(y');
        [P, Lam]=eig(var_covar);
        u=umat_(par(1:n*(n-1)/2),n);
        s_=P*sqrtm(Lam)*u;
    case 7 % scalar
        s_=0.5*eye(n); % dummy
    otherwise
        disp('Unexpected type of model.');
end
c=[];
ceq=double(any(abs(eig(s_))>1)); % equalities

end


function [lb]=setBndsCcc(type,n)
switch type
    case 1 % without spatial effects (A1=B1=S1=0)
        disp('No sense');
    case 2 % homogeneous with spatial effects (nparams=4*n+3)
        lb=[zeros(1,n) -Inf*ones(1,3*n+3) zeros(1,3*n+4:4*n+3)];
    case 3 % group-homogeneous (W_ind) with spatial effects (nparams=4*n+3*4)
        lb=[zeros(1,n) -Inf*ones(1,3*n+12) zeros(1,3*n+13:4*n+12)];
    case 4 % heterogeneous with spatial effects (nparams=7*n)
        lb=[zeros(1,n) -Inf*ones(1,5*n) zeros(1,n)];
    case 5 % full
        lb=[-Inf*ones(1,2*n*n) zeros(1,n) -Inf*ones(1,n*(n+1)/2)];
    case 6 % diagonal
        lb=[-Inf*ones(1,2*n) zeros(1,n) -Inf*ones(1,n*(n+1)/2)];
    case 7 % scalar
        lb=[-Inf*ones(1,2) zeros(1,n) -Inf*ones(1,n*(n+1)/2)];
    otherwise
        disp('Unexpected type of model.');
end
end
n=6;
th=0.14*ones(n*(n-1)/2,1);
rg=zeros(1,n); % rotation generator (http://dxdy.ru/topic56481.html)
rm=zeros(1,n); % rotation matrices (http://en.wikipedia.org/wiki/Rotation_matrix)
k=1;
for i=1:n-1
    for j=2:n
    if j>i
   rg_=zeros(n,n);
   rg_(j,i)=(-1)^(i+j);
   rg_(i,j)=-(-1)^(i+j);
     
 rm_=eye(n);
 rm_(i,i)=cos(th(k)); 
 rm_(j,j)=cos(th(k));
 rm_(i,j)=-sin(th(k));
 rm_(j,i)=sin(th(k));
    
    rg=vertcat(rg,rg_);
    rm=vertcat(rm,rm_);
   k=k+1;
    end
end
end
rg=rg(2:end,:);
rm=rm(2:end,:);

R_=zeros(n,n);
for i=1:n*(n-1)/2
    R_=R_+rg((i-1)*n+1:i*n,:);
end

Urg1=eye(n,n);
Urg2=eye(n,n);
Urm=eye(n,n);

for i=1:n*(n-1)/2
    [V,D] = eig(th(i)*R_); 
    Urg1=Urg1*(V*diag(exp(diag(D)))/V);  % second-best
    Urg2=Urg2*expm(th(i)*R_); % first-best
    Urm=Urm*rm((i-1)*n+1:i*n,:);  % from Wiki
end

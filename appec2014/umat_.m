function [ Urm ] = umat_( th, cd )
% rotation matrices (http://en.wikipedia.org/wiki/Rotation_matrix)
rm=zeros(1,cd); 
k=1;
for i=1:cd-1
for j=2:cd
    if j>i
 rm_=eye(cd);
 rm_(i,i)=cos(th(k)); 
 rm_(j,j)=cos(th(k));
 rm_(i,j)=-sin(th(k));
 rm_(j,i)=sin(th(k));
    
 rm=vertcat(rm,rm_);
 k=k+1;
    end
end
end
rm=rm(2:end,:);

Urm=eye(cd,cd);
for i=1:cd*(cd-1)/2
    Urm=Urm*rm((i-1)*cd+1:i*cd,:);
end
end


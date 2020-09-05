function l=lossFun(true,forecast)
T=length(true);
n=size(true,1);
l=zeros(1,T);
for i=1:T
%     l(1,i)=trace((true(:,:,i)-forecast(:,:,i))'*(true(:,:,i)-forecast(:,:,i)));
%     l(1,i)=trace((eye(n)/forecast(:,:,i))*true(:,:,i))-log(det((eye(n)/forecast(:,:,i))*true(:,:,i)))-n;
    l(1,i)=1/6*trace(true(:,:,i)^3-forecast(:,:,i)^3)-1/2*trace(forecast(:,:,i)^2*(true(:,:,i)-forecast(:,:,i)));
end
end
function sd=specDenEst(d)
T=length(d);
dbar=mean(d);
tau=-(T-1):T-1;
truncLag=fix(12*(T/100)^0.25);
% truncLag=14;
for i=1:length(tau)
    for t=abs(tau(i))+1:T
        c1(t)=(d(t)-dbar)*(d(t-abs(tau(i)))-dbar);
    end
    gamma(i)=1/T*sum(c1);
%         c2(i)=exp(-(tau(i)/truncLag)^2/2)*gamma(i); % Gauss lag window function
    if abs(tau(i)/truncLag)<=1
        c2(i)=gamma(i);
    else
        c2(i)=0;
    end
end
sd=sum(c2);
end
%% Preliminaries
vc(1)=vc_bekk; % hetero
vc(2)=vc_go_garch(1); % homo
vc(3)=vc_go_garch(2); % group homo
vc(4)=vc_go_garch(4); % hetero
vc(5)=vc_ccc(2); % homo
vc(6)=vc_ccc(5); % group homo
vc(7)=vc_ccc(6); % hetero

% vc(1).mes='bekk_hetero';
% vc(2).mes='go_garch_hetero';
% vc(3).mes='ccc_group_homo';

% Covariances
sTrue=zeros(n,n,T); % True covariance
for i=n+1:T
    sTrue(:,:,i)=y(:,i-n:i)*y(:,i-n:i)';
end
clear i;
[~,s{1}]=cmptLL(vc(1).x,y,W,3,1e-3); % spBekk covariance
[~,s{2}]=cmptLLspGoGarch(vc(2).x,y,W,2); % spGoGarch covariance
[~,s{3}]=cmptLLspGoGarch(vc(3).x,y,W,4); % spGoGarch covariance
[~,s{4}]=cmptLLspGoGarch(vc(4).x,y,W,3); % spGoGarch covariance
[~,s{5}]=cmptLLspCcc(vc(5).x,y,W,2); % spCcc covariance
[~,s{6}]=cmptLLspCcc(vc(6).x,y,W,3); % spCcc covariance
[~,s{7}]=cmptLLspCcc(vc(7).x,y,W,4); % spCcc covariance

%% Loss functions differentials
for i=1:7
    for j=i+1:7
        d{i,j}=lossFun(sTrue(:,:,n+1:end),s{i}(:,:,n+1:end))-lossFun(sTrue(:,:,n+1:end),s{j}(:,:,n+1:end));
    end
end
clear i j;

% Equal accuracy test
for i=1:7
    for j=i+1:7
       S(i,j)=mean(d{i,j})/sqrt(specDenEst(d{i,j})/length(d{i,j}));
    end
end
clear i j;
disp(S);

%% Cartoons
ord=[4	2	3	1	2	3	3	1	2	4	2	1	4	1];
% bekk
% a0=vc(1).x(n+1:2*n);
% b0=vc(1).x(2*n+1:3*n);
% a1=vc(1).x(4*n+1:5*n);
% b1=vc(1).x(5*n+1:6*n);

% go-garch
% a0=vc(4).x(1:n);
% b0=vc(4).x(n+1:2*n);
% a1=vc(4).x(3*n+1:4*n);
% b1=vc(4).x(4*n+1:5*n);

% ccc
a0=vc(6).x(n+1:2*n);
b0=vc(6).x(2*n+1:3*n);
% a1=vc(6).x(3*n+1:4*n);
% b1=vc(6).x(4*n+1:5*n);


Agg(1,:)=[mean(a0(ord==1)) mean(a0(ord==2)) mean(a0(ord==3)) mean(a0(ord==4))];
Agg(2,:)=[mean(b0(ord==1)) mean(b0(ord==2)) mean(b0(ord==3)) mean(b0(ord==4))];
% Agg(3,:)=[mean(a1(ord==1)) mean(a1(ord==2)) mean(a1(ord==3)) mean(a1(ord==4))];
% Agg(4,:)=[mean(b1(ord==1)) mean(b1(ord==2)) mean(b1(ord==3)) mean(b1(ord==4))];

clear ord a0 b0 a1 b1;
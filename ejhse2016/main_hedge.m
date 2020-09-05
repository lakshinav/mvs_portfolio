%% 10 Russian companies, blue chips from Micex index
% 01/01/2011 - 31/12/2013
names={'gazp','gmkn','hydr','lkoh',...
    'nvtk','rosn','sber','sngs','tatn','vtbr'};
names2={'gazr','gmkr','hydr','lkoh','notk','rosn','sbrf','sngr','tatn','vtbr'};
n=length(names);
tmp=cell(n,2);
sizes=zeros(n,1);
for i=1:n
    for j=1:2
        switch j
            case 1
                fid=fopen([upper(names{1,i}), '_110101_131231.txt']);
            case 2
                fid=fopen(['SPFB.', upper(names2{1,i}), '_110101_131231.txt']);
        end
        dd=textscan(fid, '%s %*f %*f %*f %*f %f %*f', 'delimiter', ',', 'headerlines', 1);
        fclose(fid);
        tmp{i,j}=[datenum(dd{1,1},'dd/mm/yy') dd{2}];
        clear dd fid;
    end
end
clear i j names2;

%% Make panel
for i=[3 5 9]
    l=[length(tmp{i,1}) length(tmp{i,2})];
    y=zeros(min(l),3); % {date, underlying, futures}
    if l(1)==l(2)
        y(:,1)=tmp{i,1}(:,1); % date
        y(:,2)=tmp{i,1}(:,2); % prices of the underlying
        y(:,3)=tmp{i,2}(:,2); % prices of futures
    else
        for j=1:min(l)
            ind=find(tmp{i,1}(:,1)==tmp{i,2}(j,1),1,'first');
            if not(isempty(ind))
            y(j,1)=tmp{i,1}(ind,1); % date
            y(j,2)=tmp{i,1}(ind,2); % prices of the underlying
            y(j,3)=tmp{i,2}(j,2); % prices of futures
            end
        end
    end
    eval(['yLev' upper(names{1,i}(1)) names{1,i}(2:end) '=y;']);
    y=[];
end
clear l i j ind y tmp; % 

%% Calculate log returns
for i=1:n
   y{i,1}=[y{i,1}(2:end,1) log(y{i,1}(2:end,2:3)./y{i,1}(1:end-1,2:3))];
end
clear i;

%% Plotting
set(0,'DefaultFigureWindowStyle','docked')
for i=1:n
    figure;
    subplot(2,1,1);
    plot(1:length(y{i,1}(:,2)),y{i,1}(:,2))
    subplot(2,1,2);
    plot(1:length(y{i,1}(:,3)),y{i,1}(:,3))
end
clear i;

%% Correlation between underlying and futures
C=zeros(n,2);
for i=1:n
   [R,P]=corrcoef(y{i,1}(:,2:3));
   C(i,1)=R(2,1);
   C(i,2)=P(2,1);
end
clear i R P;
% Huge significant correlation!

%% Subtract conditional mean, model is selected by BIC
yWithCondMean=y;
u=struct('arma',[],'vc',[],'ll',[],'innov',[],'sum',[]); % for the underlying
f=struct('arma',[],'vc',[],'ll',[],'innov',[],'sum',[]); % for the futures
sp=lags(0:3,2);
bminU=zeros(1,n); bminF=bminU;
for i=1:n
    for j=1:length(sp)
        [u(i,j).arma,u(i,j).vc,u(i,j).ll,u(i,j).innov,~,u(i,j).sum] = garchfit(garchset('R',sp(j,1),'M',sp(j,2),'VarianceModel','Constant'),...
            yWithCondMean{i,1}(:,2));
        [f(i,j).arma,f(i,j).vc,f(i,j).ll,f(i,j).innov,~,f(i,j).sum] = garchfit(garchset('R',sp(j,1),'M',sp(j,2),'VarianceModel','Constant'),...
            yWithCondMean{i,1}(:,3));
    end
    T=length(yWithCondMean{i,1}(:,3));
    [~,bicU]=aicbic([u(i,:).ll],(sum(sp,2)+1)',T*ones(1,length(sp)));
    [~,bicF]=aicbic([f(i,:).ll],(sum(sp,2)+1)',T*ones(1,length(sp)));
    bminU(i)=find(bicU==min(bicU));
    bminF(i)=find(bicF==min(bicF));
    
    y{i,1}(:,2)=u(i,bminU(i)).innov;
    y{i,1}(:,3)=f(i,bminF(i)).innov;
    underly(i,1)=u(i,bminU(i));
    futur(i,1)=f(i,bminF(i));
end
clear i j T bicU bicF u f;
% save('est_arma');

%%%TODO: calc signif

%% Descriptive statistics
% export DS to *.xlsx: xlswrite('filename', M)
ds=zeros(5,n);
k=3;
for i=1:n
    ds(1,i)=mean(y{i,1}(:,k));
    ds(2,i)=std(y{i,1}(:,k));
    ds(3,i)=median(y{i,1}(:,k));
    ds(4,i)=skewness(y{i,1}(:,k));
    ds(5,i)=kurtosis(y{i,1}(:,k));
end
clear i k;

%% Box plots
zz=zeros(600,n);
for i=1:n
    zz(:,i)=y{i,3}(1:600,3);
end
boxplot(zz,'plotstyle','traditional','labels',names)
fs=17;
set(gca,'FontSize',fs)
title('Распределение лог доходностей цен фьючерсов');
xlabel('Тикер акции');
ylabel('Лог доходность');
txt = findobj(gca,'Type','text');
set(txt,'FontSize',fs-4)
set(txt,'VerticalAlignment', 'bottom');
saveas(gcf,'boxplot_futures.eps') 
clear i zz txt fs;

%% Create structure to store the results
vc=struct('x',[],'fv',[],'et',[],'parsVcov',[],'Sigma_t',[],'op',[],...
    'stderr',[],'signif',[],'mes',[]);

%% Estimate CCC, DCC, GO-GARCH and RARCH by MFE toolbox
tt=0; % elapsed time
fileName='est';
opt=optimset(optimset('fmincon'), 'MaxFunEvals', 9e5,'MaxIter',5e3);
% matrix_garch(?), rcc(na figa?)

%!!!!!!!!!!!!! NADO: DECO, GO-GARCH, RARCH, CCC!!!!!!!!!!!!!!

% type={'Scalar', 'Diagonal', 'Full'}; % bekk
%     [vc(i).x,vc(i).fv,vc(i).Sigma_t,vc(i).parsVcov] = bekk(y',[],1,0,1,type{i},[],opt);

type = {'tarch','gjr'}; % ccc, dcc, gogarch, rcc
% type={'Scalar','CP','Diagonal'}; % rarch, rcc
oos=50; % for out-of-sampling
for i=1:n
    for j=1:4
        switch j
            case 1
    tic;
    [vc(i,j).x,vc(i,j).fv,vc(i,j).Sigma_t,vc(i,j).parsVcov] = ccc_mvgarch(y{i,1}(1:end-oos,2:3),[],...
        1,0,1,...
        2,[],opt);
    vc(i,j).mes=[type{2} '_ccc'];
            case 2
    [vc(i,j).x,vc(i,j).fv,vc(i,j).Sigma_t,vc(i,j).parsVcov,~,vc(i,j).op] = dcc(y{i,1}(1:end-oos,2:3),[],...
        1,0,1,...
        1,0,1,...
        2,'2-stage',[],[],opt); vc(i,j).x=vc(i,j).x';
    vc(i,j).mes=[type{2} '_dcc'];
            case 3
    [vc(i,j).x,vc(i,j).fv,vc(i,j).Sigma_t,vc(i,j).parsVcov]=gogarch(y{i,1}(1:end-oos,2:3),...
        1,1,2,[],[],opt);
    vc(i,j).mes=[type{2} '_gogarch'];
            case 4
    [vc(i,j).x,vc(i,j).fv,vc(i,j).Sigma_t,vc(i,j).parsVcov]=rarch(y{i,1}(1:end-oos,2:3),...
        1,1,'Scalar','2-stage',[],opt);
    vc(i,j).mes='scalar_rarch';
        end
    vc(i,j).et=toc;
    vc(i,j).stderr=sqrt(diag(vc(i,j).parsVcov)/length(y{i,1}));
    vc(i,j).signif=(1-normcdf(abs(vc(i,j).x./vc(i,j).stderr)))*2;
    save(fileName);
    disp([i,j]);
    tt=tt+vc(i,j).et;
    disp([num2str(tt/60) ' min.']);
    end
end
clear i j opt type;
save(fileName);
% beep

%% Log likelihoods or BICs for slides
ll=zeros(n,4);
for i=1:n
    for j=1:4
        [~,ll(i,j)]=aicbic(vc(i,j).fv,length(vc(i,j).x),length(y{i,1}));
    end
end
clear i j;

%% LR test
% only ccc and dcc are nested (perhaps rarch and gogarch too)
lr=zeros(n,2);
[~,lr(:,2),lr(:,1)] = lratiotest([vc(:,2).fv],[vc(:,1).fv],...
    length(vc(1,2).x)-length(vc(1,1).x));

%% Calculate the hedge ratios (IN-SAMPLE!!!)
hr=cell(n,size(vc,2)); % hedge ratios
for i=1:n
    for j=1:size(vc,2)
        hr{i,j}=squeeze(vc(i,j).Sigma_t(1,2,:)./vc(i,j).Sigma_t(2,2,:));
    end
end
clear i j;

%% Kernel densities of returns for hedged and unhedged investments (in-sample)
% set(0,'DefaultFigureWindowStyle','docked')
for i=8:8
    figure;
    hold on;
    cc=hsv(5);
    [kd(:,1),kd(:,2)] = ksdensity(y{i,1}(1:end-oos,2));
    plot(kd(:,2),kd(:,1),'--','LineWidth',2.5)
    t=normpdf(kd(:,2),mean(y{i,1}(1:end-oos,2)),std(y{i,1}(1:end-oos,2)));
    plot(kd(:,2),t,'k','LineWidth',2.5)
    kd=[];
    for j=1:5
        [kdHedged(:,1),kdHedged(:,2)] = ksdensity(rHedged{i,j});
        plot(kdHedged(:,2),kdHedged(:,1),'color',cc(j,:),'LineWidth',2)
        kdHedged=[];
    end
    set(gca,'FontSize',16)
    xlabel('SGNS')
    ylabel('Доходность')
    legend('Без хеджирования','Нормальное распределение',...
        'CCC','DECO','GO-GARCH','RARCH','OLS')
    hold off;
end
clear i j cc kdHedged kd t;

%% Constant conditional correlation test
% calculated in R; results are in EngShepTest.mat

%% Constant hedge ratio estimated by OLS (IN-SAMPLE!!!)
hrOls=zeros(n,3); % hr, its t-stat and mean of time-variant CCC hr
for i=1:n
    [p,R]=polyfit(y{i,1}(1:end-oos,3),y{i,1}(1:end-oos,2),1);
    hrOls(i,1)=p(1);
    R=(inv(R.R)*inv(R.R)')*R.normr^2/R.df; % var-covar matrix
    hrOls(i,2)=p(1)/sqrt(R(1,1));
    hrOls(i,3)=mean(hr{i,1});
end
clear i p R;
%TODO: plots

%% Investor's log return if hedged (IN-SAMPLE!!!)
k=5; % # of models: 4 MGARCH and OLS
rHedged=cell(n,k);
for i=1:n
    for j=1:k
        if j<k
            l=length(hr{i,j});
            rHedged{i,j}=y{i,1}(1:l,2)-hr{i,j}.*y{i,1}(1:l,3);
        else
            rHedged{i,j}=y{i,1}(1:l,2)-hrOls(i,1).*y{i,1}(1:l,3);
        end
    end
end
clear i j l k;

%% Hedge ratio effectiveness and Profit-Loss (IN-SAMPLE!!!)
k=5; % # of models: 4 MGARCH and OLS
hrEff=zeros(n,k);
plHedged=zeros(n,k);
plUnHedged=zeros(n,1);
for i=1:n
    for j=1:k
        l=length(rHedged{i,j});
        hrEff(i,j)=var(y{i,1}(1:l,2))./var(rHedged{i,j});        
        plHedged(i,j)=sum(rHedged{i,j});
    end
    plUnHedged(i,1)=sum(y{i,1}(1:l,2));
end
clear i j l k;

%% Plot hrEff
hold on;
cc=hsv(5);
for i=1:5
    plot(1:n,hrEff(:,i),'-o','color',cc(i,:),'LineWidth',1.5,...
        'MarkerFaceColor',cc(i,:),...
        'MarkerEdgeColor','k',...
        'MarkerSize',8)
end
plot(1:n,zeros(1,10),'--r')
set(gca,'XTickLabel',names,'FontSize',13)
legend('CCC','DECO','GO-GARCH','RARCH','OLS')
% title('Эффективность хеджирования');
xlabel('Тикер');
ylabel('Эффективность');
% findobj(gca,'Type','line')
hold off;
clear i cc;

%% Plot profit - loss (out-of-sample)
t=[plHedgedOos plUnHedged];
plot(1:4,sum(t,1),'-o','LineWidth',3,...
    'MarkerFaceColor','r',...
    'MarkerEdgeColor','k',...
    'MarkerSize',10)
set(gca,'XTick',1:4,'XTickLabel',{'DECO','GO-GARCH','OLS','No hedge'},...
    'FontSize',16)
xlim([0.6 4.4])
legend('Суммарная прибыль от всех активов')
% title('Эффективность хеджирования');
xlabel('Тикер');
ylabel('Прибыль');
% findobj(gca,'Type','line')
clear i t;

%%
coeff=garchfit(y{1,1}(:,2));
[SigmaForecast,MeanForecast] = garchpred(coeff,y{1,1}(:,2),20);

%% Out-of-sample
%% 1-step-ahead dynamic forecast for oos periods; got from rmgarch (GO-GARCH and DCC only)
rHedgedOos=cell(n,3);
for i=1:n
    rHedgedOos{i,1}=y{i,1}(end-oos+1:end,2)-squeeze(frcDcc(2,1,:,i)./frcDcc(2,2,:,i)).*y{i,1}(end-oos+1:end,3);
    rHedgedOos{i,2}=y{i,1}(end-oos+1:end,2)-squeeze(frcGg(2,1,:,i)./frcGg(2,2,:,i)).*y{i,1}(end-oos+1:end,3);
    rHedgedOos{i,3}=y{i,1}(end-oos+1:end,2)-hrOls(i,1).*y{i,1}(end-oos+1:end,3);
end
clear i;

%% Kernel densities of returns for hedged and unhedged investments (out-of-sample)
% set(0,'DefaultFigureWindowStyle','docked')
for i=8:8
    figure;
    hold on;
    cc=hsv(5);
    [kd(:,1),kd(:,2)] = ksdensity(y{i,1}(1:end-oos,2));
    plot(kd(:,2),kd(:,1),'--','LineWidth',2.9)
    t=normpdf(kd(:,2),mean(y{i,1}(1:end-oos,2)),std(y{i,1}(1:end-oos,2)));
    plot(kd(:,2),t,'k','LineWidth',2.9)
    kd=[];
    for j=1:2
        [kdHedged(:,1),kdHedged(:,2)] = ksdensity(rHedgedOos{i,j});
        plot(kdHedged(:,2),kdHedged(:,1),'color',cc(j,:),'LineWidth',2.5)
        kdHedged=[];
    end
    set(gca,'FontSize',16)
    xlabel('SGNS')
    ylabel('Доходность')
    legend('Без хеджирования','Нормальное распределение',...
        'DECO','GO-GARCH')
    hold off;
end
clear i j cc kdHedged kd t;

%% Hedge ratio effectiveness and Profit-Loss (out-of-sample)
k=3;
hrEffOos=zeros(n,k);
plHedgedOos=zeros(n,k);
for i=1:n
    for j=1:k
        l=length(rHedgedOos{i,j});
        hrEffOos(i,j)=var(y{i,1}(1:l,2))./var(rHedgedOos{i,j});        
        plHedgedOos(i,j)=sum(rHedgedOos{i,j});
    end
end
clear i j l k;

%% Plot hrEff (out-of-sample)
hold on;
cc=hsv(3);
for i=1:3
    plot(1:n,hrEffOos(:,i),'-o','color',cc(i,:),'LineWidth',1.5,...
        'MarkerFaceColor',cc(i,:),...
        'MarkerEdgeColor','k',...
        'MarkerSize',8)
end
set(gca,'XTickLabel',names,'FontSize',13)
legend('DECO','GO-GARCH','OLS')
% title('Эффективность хеджирования');
xlabel('Тикер');
ylabel('Эффективность');
% findobj(gca,'Type','line')
hold off;
clear i cc;

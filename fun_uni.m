 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Ниже представлена функция, рассчитанная на построение векторной        %
 % авторегрессии с двумя переменными. Функция берёт на вход 3 аргумента:  % 
 % 1) число - номер столбца из датасета, который отвечает объёмам,        %
 % 2) число - номер стобца из датасета, который отвечает ценам,           %
 % 3) Название excel файла в кавычках - датасет в формате excel.          %
 % Возвращает функция таблицы с историчексими декомпозициями объёмов      %  
 % и цен и записывает их на лист в excel. Первый столбец в excel - вклад  %
 % шока спроса, второй - вклад шока предложения, третий - вклад константы,% 
 % четвёртый - вклад первоначальных значений.                             % 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function[hd_q, hd_p]=fun(w,e, xlsx) 
data1=xlsread(xlsx); 
d1=data1(: ,w);
d2=data1(: ,e);
data=[d1,d2];
% exo=xlsread('const_prod.xlsx'); 
N=size(data,2); %Число эндогенных переменных
m=1; % Число экзогенных переменных
L=6;   %number of lags in the VAR
Y=data;
k=N*L+m;   % Число оцениваемых коэффциентов в одном уравнении
y_init=data';
y_init=y_init(:, 1:L);
X=[];
for j=1:L
X=[X lag0(data,j) ];
end
X=[X ones(size(Y,1),1)];
%X=[X ones(size(Y,1),1) exo];
Y=Y(L+1:end,:);
X=X(L+1:end,:);
T=rows(X);
irf_length=T+L;

%% Compute standard deviation of each series residual via an ols regression to be used in setting the prior
%first variable
y=Y(:,1);
x=X(:,[N*L+1 1]); % Из X берём первую переменную (то есть лаг первой переменной) и константу, которая по номеру идёт N*L+1
b0=inv(x'*x)*(x'*y);
s1=sqrt(((y-x*b0)'*(y-x*b0))/(rows(y)-2));  %std of residual standard error
%second variable
y=Y(:,2);
x=X(:,[N*L+1 2]); 
b0=inv(x'*x)*(x'*y);
s2=sqrt(((y-x*b0)'*(y-x*b0))/(rows(y)-2));  

%% Parameters to control the prior
lamda1=0.1;  %tightness prior on the AR coefficients
lamda2=0.5;
lamda3=2;   %tightness of prior on higher lags 
lamda4=100;  %tightness of prior on the constant term

%specify the prior mean of the coefficients of the Two equations of the VAR

B0=zeros((N*L+m),N);
for i=1:N
    B0(i,i)=1;
end
% Пока для экзогенных переменных и для константы предпологаем нулевое
% среднее
B0(N*L+1,1)=0; 
B0(N*L+1,2)=0;
B0=vec(B0);

%% Определим априорную ковариоционную матрицу для параметров BVAR
 % Код для ковариоционной матрицы распредления Миннесоты взять из кодов
 % BEAR Toolbox, разработанным ЕЦБ.Номера уравнений из Technical guide
 % пакетв
H=zeros(k*N,k*N);
arvar=[s1, s2]';
% set the variance on coefficients trelated to own lags, using (1.3.5)
for ii=1:N
   for jj=1:L
   H((ii-1)*k+(jj-1)*N+ii,(ii-1)*k+(jj-1)*N+ii)=(lamda1/jj^lamda3)^2;
   end
end


%  set variance for coefficients on cross lags, using (1.3.6)
for ii=1:N
   for jj=1:L
      for kk=1:N
      if kk==ii
      else
      H((ii-1)*k+(jj-1)*N+kk,(ii-1)*k+(jj-1)*N+kk)=(arvar(ii,1)/arvar(kk,1))^2*(((lamda1*lamda2)/(jj^lamda3))^2);
      end
      end
   end
end

% finally set the variance for exogenous variables, using (1.3.7)
for ii=1:N 
   for jj=1:m
   H(ii*k-m+jj,ii*k-m+jj)=(arvar(ii,1))^2*((lamda1*lamda4)^2);
   end
end


%prior scale matrix for sigma the VAR covariance
S=eye(N);
%prior degrees of freedom
alpha=N+1;

%starting values for the Gibbs sampling algorithm
Sigma=eye(N);
betaols=vec(inv(X'*X)*(X'*Y));
REPS=2000;
BURN=1000;

out1=zeros(REPS-BURN,irf_length,N); % Для irf на шок спроса
out2=zeros(REPS-BURN,irf_length,N); % Для irf на шок предложения
out3=zeros(REPS-BURN,irf_length,N); % Для irf на шок экзогенной переменной 1
out4=zeros(REPS-BURN,irf_length,N); % Для irf на константу
out6=zeros(REPS-BURN,1,N*k); % Для хранения сэмплов коэффициентов beta;
out7=zeros(REPS-BURN,N,N); % Для хранения сэмплов коэффициентов дисперсий;
fevd_record1=zeros(REPS-BURN,12,N); % Для FEVD вклад шока спроса
fevd_record2=zeros(REPS-BURN,12,N);% Для FEVD вклад шока предложения
%fevd_record3=zeros(REPS-BURN,12,N);  %Для FEVD вклад щока экзогенной переменной 1
%fevd_record4=zeros(REPS-BURN,12,N);  %Для FEVD вклад щока экзогенной переменной 2
hd_record1=zeros(REPS-BURN,irf_length-L,N); % Для HD вклад шока спроса
hd_record2=zeros(REPS-BURN,irf_length-L,N); % Для HD вклад шока предложения
%hd_record3=zeros(REPS-BURN,irf_length-L,N); % Для HD вклад шока экзогенной переменной 1
hd_record4=zeros(REPS-BURN,irf_length-L,N); % Для HD вклад константы
%hd_record5=zeros(REPS-BURN,irf_length-L,N); % Для HD вклад константы
hd_record6=zeros(REPS-BURN,irf_length-L,N); % Для HD вклад первоначальных значений
likelihood=[];
%% Гиббс сэмплинг

jj=1;
for i=1:REPS
M=inv(inv(H)+kron(inv(Sigma),X'*X))*(inv(H)*B0+kron(inv(Sigma),X'*X)*betaols);
V=inv(inv(H)+kron(inv(Sigma),X'*X));
%check for stability of the VAR
check=-1;
while check<0 
beta=M+(randn(1,N*(N*L+m))*chol(V))';
CH=stability(beta,N,L,m);
if CH==0
    check=10;
end
end
%draw sigma from the IW distribution
e=Y-X*reshape(beta,N*L+m,N);
%scale matrix
scale=e'*e+S;
Sigma=IWPQ(T+alpha,inv(scale));
    
    if i>=BURN
        %impose sign restrictions
        chck=-1;
        while chck<0
            K=randn(N,N);
            Q=getqr(K);
            A0hat=chol(cov(e));
            A0hat1=(A0hat*Q);  %candidate draw
            %check signs
            e1=A0hat1(1,1)>0;  %отклик q на шок спроса
            e2=A0hat1(1,2)>0;  %отклик q на шок предложения
            e3=A0hat1(2,1)>0;  %отклик p на шок спроса
            e4=A0hat1(2,2)<0;  %отклик p на шок предложения
            if e1+e2+e3+e4==4
                chck=10;
            else
                e1=-A0hat1(1,1)>0;
                e2=-A0hat1(1,2)>0;  
                e3=-A0hat1(2,1)>=0;
                e4=-A0hat1(2,2)<0; 
                if e1+e2+e3+e4==4;
                    A0hat1(1,1:N)=-A0hat1(1,1:N);
                    A0hat1(2,1:N)=-A0hat1(2,1:N);
                    chck=10;
                end
           end
        end
        
 % Далее считаем имульсные отклики       
yhat1=zeros(irf_length,N); %+1 для экзо
vhat1=zeros(irf_length,N); %+1 для экзо
vhat1(L+1,1)=1; %Шок спроса

for j=L+1:irf_length
 yhat1(j,:)=[reshape(wrev(yhat1(j-L:j-1,:))', 1, L*N) 0 ]*reshape(beta,N*L+m,N)+vhat1(j,:)*A0hat1;
end



yhat2=zeros(irf_length,N);
vhat2=zeros(irf_length,N);
vhat2(L+1,2)=1; % Шок предложения
for j=L+1:irf_length
 yhat2(j,:)=[reshape(wrev(yhat2(j-L:j-1,:))', 1, L*N) 0]*reshape(beta,N*L+m,N)+vhat2(j,:)*A0hat1;
end



yhat4=zeros(irf_length,N);
vhat4=zeros(irf_length,N);
vhat4=zeros(irf_length,1);
vhat4(L+1, 1)=1; % "Шок" константы

for j=L+1:irf_length
 yhat4(j,:)=[reshape(wrev(yhat4(j-L:j-1,:))', 1, L*N)  vhat4(j, :) ]*reshape(beta,N*L+m,N);
end




 
 out1(jj,:,:)=yhat1; % Отклики на шок спроса, 
 %1-й индекс в out1 отвечает за итерацияю сэмплирования Гибса, второй за
 %горизонт прогноза отклика, 3 за номер переменной
 out2(jj,:,:)=yhat2; % Отклики на шок предложения
 out4(jj,:,:)=yhat4; % Отлклик на "шок" константы
 out6(jj,:,:)=beta;
 out7(jj,:,:)=Sigma;
 l=loglik(reshape(beta,N*L+m,N),Sigma,Y,X);
 likelihood=[likelihood;l];
 %fevd
 ETA=(inv(A0hat1)*e')'; %Матрица структурных шоков T на N
 sigma_1=var(ETA(:, 1));
 sigma_2=var(ETA(:, 2));
 fevd1=zeros(12,N);
 fevd2=zeros(12,N);

for h=L+1:L+12
    fevd1(h-L, :)=[sigma_1*sum(out1(jj, L+1:h, 1).^2), sigma_1*sum(out1(jj, L+1:h, 2).^2)]; % вклады  шока спроса в FEVD для первой и второй переменной
    fevd2(h-L, :)=[sigma_2*sum(out2(jj, L+1:h, 1).^2), sigma_2*sum(out2(jj, L+1:h, 2).^2)]; % вклады шока предложения в FEVD для первой и второй переменной
end

hd1=zeros(irf_length-L, N); %Для каждой переменной считаем накопленный вклад шока спроса
hd2=zeros(irf_length-L, N); %Для каждой переменной считаем накопленный вклад  шока предложения
hd4=zeros(irf_length-L, N); %Для каждой переменной считаем накопленный вклад  константы
hd6=zeros(irf_length-L, N); %Для каждой переменной считаем накопленный вклад  первоначальных значений

%beta1=beta[1:N*L+2, 1];
%beta2=beta[N*L+2+1:lengh(beta), 1];
%for i=1:L
  %  A
for t=L+1:irf_length
    hd1(t-L, :)=[ETA(1:t-L, 1)'*wrev(out1(jj, L+1:t, 1))', ETA(1:t-L, 1)'*wrev(out1(jj, L+1:t, 2))']; %Вклады шока спроса в HD первой и второй переменной
    hd2(t-L, :)=[ETA(1:t-L, 2)'*wrev(out2(jj, L+1:t, 1))', ETA(1:t-L, 2)'*wrev(out2(jj, L+1:t, 2))']; % Вклады шока предложения в HD первой и второй переменной
    hd4(t-L, :)=[X(1:t-L,N*L+1)'*wrev(out4(jj, L+1:t, 1))', X(1:t-L, N*L+1)'*wrev(out4(jj, L+1:t, 2))']; %Вклады константы в HD первой и второй переменной
    beta1=beta(1:N*L, 1);
q1=N*L+m+1;
q2=(N*L+m)*N-m;
beta2=beta(q1:q2, 1);
beta3=[beta1'; beta2'];
%t=5;
A=zeros(2, 1);
for i=1:L
    q1=2*i-1;
    q2=2*i;
    A=A+beta3(:, q1:q2)^(t-L)*y_init(:, -i+1+L);
end
     hd6(t-L, :)=A'; 
end

 hd_record1(jj, :, :)=hd1;
 hd_record2(jj, :, :)=hd2;
 hd_record4(jj, :, :)=hd4;
 hd_record6(jj, :, :)=hd6;
 fevd_record1(jj, :, :)=fevd1;
 fevd_record2(jj, :, :)=fevd2;
   jj=jj+1; 

 
   end 
end
betam=squeeze(mean(out6,1));
sigmam=squeeze(mean(out7,1));
lik=loglik(reshape(betam,N*L+m,N),sigmam,Y,X);
BIC=log(N)*(N*L+m)-2*lik;


  
%% Возможности для графического вывода
%figure(1);
 % subplot(4,3,1)
%  temp=out1(:,:,1);
 % temp1=squeeze(prctile(temp,[50 16 5 95 84],1))';
%  plot(temp1(L+1:L+13,:), '--', 'Color',"#7E2F8E", 'LineWidth',1);
% title('Q, demand shock');
 % axis tight
  
  
  
 % subplot(4,3,2)
 % temp=out1(:,:,2);
 % temp1=squeeze(prctile(temp,[50 16 5 95 84],1))';
 % plot(temp1(L+1:L+20,:), '--', 'Color',"#7E2F8E", 'LineWidth',1 );
 % title('P, demand shock');
 %axis tight
  
 %  subplot(4,3,3)
% temp=out2(:,:,1);
% temp1=squeeze(prctile(temp,[50 16 5 95 84],1))';
% plot(temp1(L+1:L+13,:), '--', 'Color',"#7E2F8E", 'LineWidth',1);
% title('Q, supply shock');
 % axis tight
  
 %   subplot(4,3,4)
 % temp=out2(:,:,2);
% temp1=squeeze(prctile(temp,[50 16 5 95 84],1))';
 % plot(temp1(L+1:L+13,:), '--', 'Color',"#7E2F8E", 'LineWidth',1);
 % title('P, supply shock');
 %axis tight
  
  
 % subplot(4,3,5)
 % temp=out3(:,:,1);
 % temp1=squeeze(prctile(temp,[50 16 84],1))';
 % plot(temp1(L+1:L+13,:));
 % title('Q, exchange rate shock');
 % axis tight
  
  
   % subplot(4,3,6)
%  temp=out3(:,:,2);
 % temp1=squeeze(prctile(temp,[50 16 84],1))';
 % plot(temp1(L+1:L+13,:));
 % title('P, exchange rate shock');
 % axis tight

  
  
     %subplot(4,3,7)
  %h=[1;2;3;4;5;6;7;8;9;10;11;12];
 % temp1=fevd_record1(:,:,1);
  %temp2=fevd_record2(:, :, 1);
  %temp3=fevd_record3(:, :, 1);
  %c=squeeze(prctile(temp1,[50],1))'+squeeze(prctile(temp2,[50],1))';
 %temp11=squeeze(prctile(temp1,[50],1))'./c;
 %temp21=squeeze(prctile(temp2,[50],1))'./c ;
 % temp31=squeeze(prctile(temp3,[50],1))'./c;
 % bar(h, [temp11, temp21], "stacked");
 % title('FEVD Q');
 % axis tight
  
   % subplot(4,3,8)
 % h=[1;2;3;4;5;6;7;8;9;10;11;12];
 % temp1=fevd_record1(:,:,2);
 % temp2=fevd_record2(:, :, 2);
 % temp3=fevd_record3(:, :, 2);
 % c=squeeze(prctile(temp1,[50],1))'+squeeze(prctile(temp2,[50],1))';
 % temp11=squeeze(prctile(temp1,[50],1))'./c;
 % temp21=squeeze(prctile(temp2,[50],1))'./c ;
 % temp31=squeeze(prctile(temp3,[50],1))'./c;
 %bar(h, [temp11, temp21], "stacked");
 % title('FEVD P');
 % axis tight
  
%% Исторические декомпозиции: графики и сами декмопозиции в виде таблицы для выгрузки 
%figure(2);
A = datetime(2010,2+L,01);
B = datetime(2022,10,01);
C=A:calmonths(1):B;
  temp1=hd_record1(:,:,1);
  temp2=hd_record2(:, :, 1);
  temp4=hd_record4(:, :, 1);
  temp6=hd_record6(:, :, 1);
  temp11=squeeze(prctile(temp1,[50],1))';
  temp21=squeeze(prctile(temp2,[50],1))' ;
  temp41=squeeze(prctile(temp4,[50],1))';
  temp61=squeeze(prctile(temp6,[50],1))';
 % bar(C, [temp11, temp21, temp31, temp41, temp51], "stacked");
  %hold on
  %plot(C, c, 'k','LineWidth',1.6);
 % hold on
 % plot(C, data(L+1:153, 1), 'k','LineWidth',1.6);
 % title('Q');
 % axis tight;
 % legend('demand','supply', "exch", "const", "init");
 % hold off;
  
hd_q=[temp11 temp21 temp41 temp61];

  
  %figure(3);
  A = datetime(2010,2+L,01);
  B = datetime(2022,10,01);
  C=A:calmonths(1):B;
  temp1=hd_record1(:,:,2);
  temp2=hd_record2(:, :, 2);
  temp4=hd_record4(:, :, 2);
  temp6=hd_record6(:, :, 2);
  temp11=squeeze(prctile(temp1,[50],1))';
  temp21=squeeze(prctile(temp2,[50],1))' ;
  temp41=squeeze(prctile(temp4,[50],1))';
  temp61=squeeze(prctile(temp6,[50],1))';
  %bar(C, [temp11, temp21, temp31, temp41, temp51], "stacked");
  %hold on
  %plot(C, c, 'k','LineWidth',1.6);
 % hold on
 % plot(C, data(L+1:153, 2), 'k','LineWidth',1.6);
 % title('P');
 % axis tight;
 % legend('demand','supply', "exch", "const", "init");
 % hold off;

  
  hd_p=[temp11 temp21 temp41 temp61];
    
 %% Расчёты DIC и BIC 
  betam=squeeze(mean(out6,1));
sigmam=squeeze(mean(out7,1));
D_mean=-2*loglik(reshape(betam,N*L+m,N),sigmam,Y,X);
D=squeeze(mean(-2*likelihood));% the mean of the likelihood evaluated at each saved draw
params=D-D_mean;

%Calculate the DIC 
dic=D+2*params;
BIC=log(N)*(N*L+m)+(squeeze(prctile(-2*likelihood,[50],1)));
dic
 
end
 


function [y,x] = create_dummies(lamda,tau,delta,epsilon,p,mu,sigma,n)
% Creates matrices of dummy observations [...];
%lamda tightness parameter
%tau  prior on sum of coefficients
%delta prior mean for VAR coefficients
% epsilon tigtness of the prior around constant
% mu sample mean of the data
% sigma AR residual variances for the data




% Initialise output (necessary for final concatenation to work when tau=0):
x = [];
y = [];
yd1 = [];
yd2 = [];
xd1 = [];
xd2 = [];

%% Get dummy matrices in equation (5) of Banbura et al. 2007:
if lamda>0
    if epsilon >0
	yd1=[diag(sigma.*delta)./lamda;
         zeros(n*(p-1),n);
         diag(sigma);
         zeros(1+1,n)]; % короче вот эта одна строчка для константы походу, поэтому добавил 1+1, теперь две строчки одна из которых для экзогенки
     
	jp=diag(1:p);
    
	xd1=[kron(jp,diag(sigma)./lamda) zeros((n*p),1+1); % +1 для экзо
         zeros(n,(n*p)+2); %здесь было n*p+1, сделали на n*p+2 для экзо
         zeros(1,n*p) epsilon 0 ; 
         zeros(1,n*p) 0 epsilon]; %добавил последнюю строчку для экзогенной переменной
   

    else % а далее идёт код, если бы не было константы
    
   yd1=[diag(sigma.*delta)./lamda;
         zeros(n*(p-1),n);
         diag(sigma);
         zeros(1,n)]; %zeros(1,n) добавлена строчка для экзо
     
	jp=diag(1:p);
    
	xd1=[kron(jp,diag(sigma)./lamda);
         zeros(n,(n*p+1))]; %n*p+1 для экзо
    end
end
%% Get additional dummy matrices - see equation (9) of Banbura et al. 2007:
if tau>0
    if epsilon>0
	yd2=diag(delta.*mu)./tau;
	xd2=[kron(ones(1,p),yd2) zeros(n,1+1)]; %+1 для экзо
    else
      yd2=diag(delta.*mu)./tau;
	xd2=[kron(ones(1,p),yd2) zeros(n,1)];  %zeros(n,1) для экзо
    end
end
     
%% 
y=[yd1;yd2];
x=[xd1;xd2];
 
         
 
 
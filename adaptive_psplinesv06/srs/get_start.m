

x=unifrnd(0,1,1000,1);
e = normrnd(0,sqrt(0.02),1000,1);
t = exp(-400*(x-.6).^2) + (5/3)* exp(-500*(x-0.75).^2) + 2* exp(-500*(x-.9).^2);

y=t+e;



%load data.txt;
%x = data(:,1);
%y = data(:,2);
degree  = 2;
nknots = 20;

fit3 = srslocal(x,y,2,20);

disp('Estimate of alpha');

disp(exp(fit3.alpha0));

function X = gen_pspline_basis(x,knots,degree);

% Generates P-spline basis

nknots=length(knots);
n=length(x);

for i=1:nknots,
 temp=(x-knots(i));
 temp(temp<0)=0;
 X(:,i)=temp.^degree;
 clear temp;
end;

power_series = ones(n,1);

for i = 1:degree,
    power_series =  [power_series x.^i];
end;

X=[power_series X];  


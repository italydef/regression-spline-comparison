function fmin=spline(x,y,alpha)

% Pre:
% Takes the vector of y values (the observed responses), and 
% the vextor of x values (the explanatory variable)
% and returns the vector fmin of fitted values using fixed 
% knots natural cubic smoothing splines and GCV to choose 
% smoothing parameter.  Alpha is the grid of values for the
% grid search.

% Note: the vector A in this program allows you to track the 
% GCV scores corresponding to the grid of alpha values.  The first 
% has the alpha value used in smoothing, and the associated GCV 
% score is in the second column.  

n=length(x);
m=length(alpha);
gcv=1e10;
fmin=[];

A=zeros(m,2);

for i=1:m

[Q,B,f]=smooth(x,y,alpha(i));

diff=(f-y);
ss=diff'*diff;

A(i,1:2)=[alpha(i),gcvtest(Q,B,ss,n,alpha(i))];

 if A(i,2)<gcv
    fmin=f;
    gcv=A(i,2);
    aopt=A(i,1);
 end;

end;







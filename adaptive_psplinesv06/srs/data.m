function [X,Y] = data()
% reads points from screen and stores them in matrices X, Y, sorted by rows;


colormap([1 1 1]);
image(1);
axis([0 10 0 10]);
axis xy;
black=[0 0 0];
i=1;
[X(i),Y(i),button]=ginput(1);
text(X(i),Y(i),'o','color',black);
while button < 3,
    i=i+1;
   [a,b,button]=ginput(1);
   if button < 3,
      X(i) = a; Y(i) = b;
      text(X(i),Y(i),'o','color',black);
   end;
end;
X = X';  Y = Y';
A = sortrows([X Y]);
X = A(:,1);
Y = A(:,2);

function r=ranks(y) ;
[n,m]=size(y) ;
for j = 1:m ;
temp = sortrows([y(:,j) (1:n)'],1) ;
temp = [temp (1:n)'] ;
temp = sortrows(temp,2) ;
r(:,j) = temp(:,3) ;

end ;

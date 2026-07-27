

 function [c,d]=binning(X1,X2,Y,n,A1,A2,M1,M2,L1,L2,de1,de2)
  
  % Performs linear binning on the set of observations.
  % We pad the matrices c and d with zeros, so that when
  % using the macro convreg.m we do not have to verify if
  % the element in the matrix c or d is defined at a 
  % particular index. 
  c=zeros(M1+2*L1+1,M2+2*L2+1);
  d=zeros(M1+2*L1+1,M2+2*L2+1);

  for k=1:n
   
   g1=floor((X1(k)+A1/2)/de1)+1+L1;
   g2=floor((X2(k)+A2/2)/de2)+1+L2;

   length1=X1(k)-de1*(g1-M1/2-L1-1);
   length2=X2(k)-de2*(g2-M2/2-L2-1);

   % The weights of the surrounding gridpoints are computed.
   area1=(de1-length1)*(de2-length2)/(de1*de2);
   area2=length1*(de2-length2)/(de1*de2);
   area3=length1*length2/(de1*de2);
   area4=(de1-length1)*length2/(de1*de2);

   % We update the matrices c and d.
   c(g1,g2)=c(g1,g2)+area1;
   d(g1,g2)=d(g1,g2)+Y(k)*area1;

   c(g1+1,g2)=c(g1+1,g2)+area2;
   d(g1+1,g2)=d(g1+1,g2)+Y(k)*area2;

   c(g1+1,g2+1)=c(g1+1,g2+1)+area3;
   d(g1+1,g2+1)=d(g1+1,g2+1)+Y(k)*area3;

   c(g1,g2+1)=c(g1,g2+1)+area4;
   d(g1,g2+1)=d(g1,g2+1)+Y(k)*area4;
 
  end 



clc;
clear all;


AA=[1 0 0 1 0 0 1 0 0 0 0 0 -1 0 0 ;
   0 1 0 0 1 0 0 1 0 0 0 0 0 -1 0 ;
   0 0 1 0 0 1 0 0 1 0 0 0 0 0 -1 ;
   -1 0 0 0 0 0 0 0 0 1 0 0 0 0 0 ;
   0 -1 0 0 0 0 0 0 0 1 0 0 0 0 0 ;
   0 0 -1 0 0 0 0 0 0 1 0 0 0 0 0 ;
   0 0 0 -1 -1 0 0 0 0 0 1 0 0 0 0 ;
   0 0 0 -1 0 -1 0 0 0 0 1 0 0 0 0 ;
   0 0 0 0 -1 -1 0 0 0 0 1 0 0 0 0 ;
   0 0 0 0 0 0 -1 -1 -1 0 0 1 0 0 0 ;[-eye(9),zeros(9,6)]]



%Four      r2      r3      r4      H1      R1
A=[1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 -1 0 0 0;
   0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 -1 0 0;
   0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 -1 0;
   0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 1 0 0 0 0 0 0 0 -1;
   -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0;
   0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0;
   0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0;
   0 0 0 -1 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0 0;
   0 0 0 0 -1 -1 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0;
   0 0 0 0 -1 0 -1 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0;
   0 0 0 0 -1 0 0 -1 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0;
   0 0 0 0 0 -1 0 -1 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0;
   0 0 0 0 0 0 -1 -1 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0;
   0 0 0 0 0 -1 -1 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0;
   0 0 0 0 0 0 0 0 -1 -1 -1 0 0 0 0 0 0 0 1 0 0 0 0 0;
   0 0 0 0 0 0 0 0 -1 -1 0 -1 0 0 0 0 0 0 1 0 0 0 0 0;
   0 0 0 0 0 0 0 0 -1 0 -1 -1 0 0 0 0 0 0 1 0 0 0 0 0;
   0 0 0 0 0 0 0 0 0 -1 -1 -1 0 0 0 0 0 0 1 0 0 0 0 0;
   0 0 0 0 0 0 0 0 0 0 0 0 -1 -1 -1 -1 0 0 0 1 0 0 0 0;[-eye(24)]]%,zeros(20,4)]] 
 %A=[AA(:,13:16),AA(:,1:12),AA(:,17:24)]
%A=[AA(:,7:9),AA(:,1:6),AA(:,10:15)];
A=[ 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0;
    0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0;
    0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 0 0 -1 0;
    0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 -1;
     -1 0 0 0  0 0 0 0  0 0 0 0  0 0 0 0  1 0 0 0  0 0 0 0;
     0 0 0 0  -1 0 0 0  0 0 0 0  0 0 0 0  1 0 0 0  0 0 0 0;
     0 0 0 0  0 0 0 0  -1 0 0 0  0 0 0 0  1 0 0 0  0 0 0 0;
     0 0 0 0  0 0 0 0  0 0 0 0  -1 0 0 0  1 0 0 0  0 0 0 0;
     0 -1 0 0  0 -1 0 0  0 0 0 0  0 0 0 0  0 1 0 0  0 0 0 0;
     0 0 0 0  0 -1 0 0  0 -1 0 0  0 0 0 0  0 1 0 0  0 0 0 0;
     0 -1 0 0  0 0 0 0  0 -1 0 0  0 0 0 0  0 1 0 0  0 0 0 0;
     0 -1 0 0  0 0 0 0  0 0 0 0  0 -1 0 0  0 1 0 0  0 0 0 0;
     0 0 0 0  0 -1 0 0  0 0 0 0  0 -1 0 0  0 1 0 0  0 0 0 0;
     0 0 0 0  0 0 0 0  0 -1 0 0  0 -1 0 0  0 1 0 0  0 0 0 0;
     0 0 0 0  0 0 -1 0  0 0 -1 0  0 0 -1 0  0 0 1 0  0 0 0 0;
     0 0 -1 0  0 0 0 0  0 0 -1 0  0 0 -1 0  0 0 1 0  0 0 0 0;
     0 0 -1 0  0 0 -1 0  0 0 0 0  0 0 -1 0  0 0 1 0  0 0 0 0;
     0 0 -1 0  0 0 -1 0  0 0 -1 0  0 0 0 0  0 0 1 0  0 0 0 0;
     0 0 0 -1  0 0 0 -1  0 0 0 -1  0 0 0 -1  0 0 0 1  0 0 0 0;[-eye(16),zeros(16,8)]]

%  r11      r21     r31     r41     R1      Rs1
%  A=[1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0 0;
%     0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0;
%     0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 0 0 0 -1 0;
%     0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 0 0 0 0 0 0 0 -1;
%     -1 0 0 0 -1.5 0 0 0 -2 0 0 0 -.25 0 0 0 1 0 0 0 0 0 0 0;
%     0 -1/3 0 0 0 -1 0 0  0 -2/3 0 0 0 -1/6 0 0 0 1 0 0 0 0 0 0;
%     0 0 -1/3 0 0 0 -1 0 0 0 -1 0 0 0 -1/4 0 0 0 1 0 0 0 0 0;
%     0 0 0 -1 0 0 0 -3 0 0 0 -3 0 0 0 -1 0 0 0 1 0 0 0 0;[-eye(16),zeros(16,8)]]
% 

%  A=[1 1 1 0 0 0 0 0 0 0 0 0 -1 0 0;
%     0 0 0 1 1 1 0 0 0 0 0 0 0 -1 0;
%     0 0 0 0 0 0 1 1 1 0 0 0 0 0 -1;
%     -1 0 0 -1 0 0 -1/3 0 0 1 0 0 0 0 0;
%     0 -.5 0 0 -1 0 0 -1/3 0 0 1 0 0 0 0;
%     0 0 -1 0 0 -2 0 0 -1 0 0 1 0 0 0;[-eye(9),zeros(9,6)]]

% 
sinks=4;



% 


% 
% b=zeros(6,1);
% d=[A,b]
% a1 = reduce_rows(d)
% [X,Y]= fourmotz(A,b,6)

%The program Does forier Motzkin elimination
%This program puts Variable which should be eliminated in first rows
%Number of eliminated variables e
%A is in the form [A|b]

%A=[1 0 -1  0;1  -1 0 0;-1 2 1 0;0 -1 0 -1];
e=16;   %for 4
%e=9    %for 3
[m n]=size(A);


for i=1:e
rp=find(A(:,1)>0);
rn=find(A(:,1)<0);

[s s1]=size(rn);
[t t1]=size(rp);
%temp1=zeros(1,n);
for ii=1:s
    for j=1:t
    %newsize=[m-s+s*t n-1]
%     if A(i,1)>0
%         temp1=[temp1;A(i,:)];

        A=[A;A(rn(ii),:)+abs(A(rn(ii),1))/A(rp(j),1)*A(rp(j),:)];
    end
end

%A(rp,:)=0;
%A(rn,:)=0;
A([rp;rn],:)=[];
%A(rn,:)=[];
A(:,1)=[];
%A(rp(i,1),:)=[];
 threshold=size(A)*[1 0]';
 size(A)
 
 for y=1:threshold
     
     if A(y,1)~=0
         A(y,:)=A(y,:)/abs(A(y,1));
     elseif (A(y,1)==0)&(max(A(y,:))~=0)
         A(y,:)=A(y,:)/max(abs(A(y,:)));
     end
 end
 % while threshold==s*t
if A(threshold-1,:)<=A(threshold,:)
    A(threshold-1,:)=[];
end 
threshold=threshold-1;
end

A=sortrows(A);
   % A=reduce_rows(A);
   
  size(A)  
for q=0:size(A)*[1 0]'-1
   
    for k=1:(size(A)*[1 0]'-q-1)
       
        
       if (A(k,:)<=A(size(A)*[1 0]'-q,:))==ones(1,size(A)*[0 1]')
            A(k,:)=[];
    
       end
   % A=reduce_rows(A);
       if k>=size(A)*[1 0]'-q-1
        break;
       end
    end
end
[r c]=size(A);


if r>36
    
    flg=0;
    i=1;
  while i< (size(A)*[1;0])
    
    Temp=A;
       
       if flg==1
           i=i-1;
           flg=0;
       end
       Temp(i,:)=[];
       [x,val,exitflag]=linprog(zeros(1,c)',Temp,zeros(r-1,1),A(i,:),0.0001);
       %[y,valy,exitflagy]=linprog(zeros(1,c),Temp,zeros(r-1,1),A(i,:),5);

      if (exitflag==-2 |exitflag==-5)
        %if (x>=0.000001 | x<=-0.000001)& (exitflag~=1)
           A(i,:)=[];
           flg=1;
       end
      i=i+1; 
      [r c]=size(A);

  end
 end

% size(A)
B=[];
A=sortrows(A);
% b=zeros(size(A)*[1 0]',1);
% [A,b] = fourmotz(A,b,8);
for l=1:size(A)*[1 0]'
    if (A(l,1:sinks)==-A(l,sinks+1:2*sinks))==ones(1,sinks)
        B=[B;A(l,:)];
            
    end
    
end
            
%printing rational form
d=sinks+1;
[r c]=size(A);
tempi=zeros(1,sinks);
defecit=0
while d<size(A)*[1 0]'
   for i=sinks+1:r 
    if (i~=d&(A(i,sinks+1:2*sinks)<0)&(A(d,sinks+1:2*sinks)<0)==ones(1,sinks) )   % then not all the A(i,5:8)>=A(d,5:8), because otherwise we have removed inthe algorithm
        ratio=min(A(d,sinks+1:2*sinks)./A(i,sinks+1:2*sinks));
        tempi= ratio*A(i,sinks+1:2*sinks);
        tempd=A(d,sinks+1:2*sinks)
        %cd=find(A(d,sinks+1:2*sinks)>A(i,sinks+1:2*sinks));
        ci=find(A(i,sinks+1:2*sinks)>A(d,sinks+1:2*sinks));
        %defecit1=sum(A(d,cd)-A(i,cd))
        defecit=sum(tempi(ci)-tempd(ci))
            if (defecit+ratio>=1) & (ratio*A(i,2:sinks)>=A(d,2:sinks)==ones(1,sinks-1))
                 A(d,:)=[];
                 d=d-1;
                 break;
            end
    end
    
   end
  [r c]=size(A);

    d=d+1
end
     
k=zeros(1,8)
for i=1:size(A)*[1;0]
    minA=min(abs(A(i,:)))
    if minA~=0
      A(i,:)=A(i,:)/minA
     
    end
    
end
A=sortrows(A)

 % A = reduce_rows(A);
 for q=1:size(A)*[1 0]'
    
    flag=0;
    for k=1:(size(A)*[1 0]'-q-1)
       if flag==1
           k=k-1;
           flag=0;
       elseif k>=q
           break;
       end 
       if (A(k,:)<=A(size(A)*[1 0]'-q,:))==ones(1,size(A)*[0 1]')
      %  if (k~=i & A(k,:)<=A(q,:))==ones(1,size(A)*[0 1]')
            A(k,:)=[];
            flag=1;
       else
           flag=0;
    
       end
   % A=reduce_rows(A);
       if k>=size(A)*[1 0]'-q-1
        break;
       end
    end
end


            
            


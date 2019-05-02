E=[0,1,0,0,1,1,1,1,1,0,0,1,1];
R=[0,0,1];
disp('E quotient')
disp(fliplr(E))
disp('R remainder')
disp(fliplr(R))
% E=[1  1 1]
 
% R=[1 1 1 1];
Ak=[1,0,1,1,0,1,1];
B=[0,1,0,0,1];
 
%Step 1:E*B using booth algorithm for multiplication
Q=B;
M=E;
n=length(Q);
Q=fliplr(Q) ;
M=fliplr(M);
A=zeros(1,length(M));
C=0;
I=n;
 
while(I~=0)  
    if(Q(length(Q))==1)
        op=adds(A,M);
        C=op(1);
        A=op(2:end);
    end
    J=[C,A,Q];
    J=[0,J(1:end-1)];
    C=J(1);
    A=J(2:1+length(M));
    Q=J(length(M)+2:length(M)+length(Q)+1);
    I=I-1;
    
end
mul=[A Q];
mul=fliplr(mul);
if(length(mul)>length(R))
    R(numel(mul))=0;
else
    if(length(mul)<length(R))
    mul(numel(R))=0;
    end 
end
R=fliplr(R);
mul=fliplr(mul);
fop=adds(mul,R);
F= fop(length(fop)-length(Q)-length(A)+1:end);
disp('Multiplying and adding remainder')
disp(F)
 
%2's
F=fliplr(F);
for i=1:length(F)
    if(F(i)==1)
        i=i+1;
        while(i<=length(F))
            F(i)=~F(i);
            i=i+1;
        end
        break;
    end
end
disp('2 complement')
F=fliplr(F);
disp(F)
 
%subtract A
F=fliplr(F);
if(length(F)>length(Ak))
    Ak(numel(F))=0;
else
    if(length(F)<length(Ak))
    F(numel(Ak))=0;
    end 
end
Ak=fliplr(Ak);
AC=zeros(1,length(Ak));
    for i=length(Ak):-1:1
    if(Ak(i)==1)
        AC(i)=Ak(i);
        i=i-1;
        while(i>0)
            AC(i)=~Ak(i);
            i=i-1;
        end
        break;    
    end
    AC(i)=Ak(i);
    end
    F=fliplr(F);
    AC;
    TEMP=adds(F,AC);
    TEMP=TEMP(2:end);
    disp('Subtract A')
disp(TEMP)
    %TEMP= TEMP(find(TEMP==1,1):end
    TEMP=fliplr(TEMP);
    disp('Reverse')
disp(TEMP)
    ones(1,length(TEMP));
    TEMP=adds(TEMP,ones(1,length(TEMP)));
    TEMP=TEMP(2:end);
    disp('Subtract 1')
disp(TEMP)
    fliplr(TEMP);
    disp('Decrypted output')
    disp(TEMP)
 
 
 
    
function [op]=adds(x,y)
z(length(x)+1)=0;
for i=length(x):-1:1
    s(i)=xor(xor(x(i),y(i)),z(i+1));
    z(i)=or(or(and(z(i+1),x(i)),and(z(i+1),y(i))),and(x(i),y(i)));
end
c=z(1);
op=[c s];
end
    
    


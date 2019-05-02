close all;
clear all;
 
%Cryptographic encryption 
a=fliplr(huffadapt());
n=length(a);
%a=[1,0,1,0,0,1,0,0,0,1]; %encoded input
A=[1,0,1,1,0,1,1];
B=[0,1,0,0,1];
disp('Encoded input:')
disp(fliplr(a))
%STEP 1: Add 1 to a
for i=1:n
    if(a(i)==0)
    a(i)=1;
    break;
    else
    a(i)=0;
    if(i==5)
        n=n+1;
        a(n)=1;
    end
    end
end 
a;
disp('Adding 1')
disp(fliplr(a))
a=fliplr(a);
 
%Step 2:Reversing a
disp('Reversing')
disp(fliplr(a))
 
%Step 3:Add A
 
A;
if(length(a)>length(A))
    A(numel(a))=0;
else
    if(length(a)<length(A))
    a(numel(A))=0;
    end 
end
 
OP=adds(fliplr(A),fliplr(a));
disp('Adding private key A')
disp(OP)
 
% for i=1:length(a)
%     s1(i)=xor(c(i),xor(a(i),A(i)));
%     c(i+1)=or(or(and(c(i),a(i)),and(c(i),A(i))),and(a(i),A(i)));
% end
% s1
% cout=c(length(a)+1)
% OP=[s1,cout]
% c1=and(a,A)
% s=xor(s1,c1)
% c2=or(c1,and(s1,c1))
 
% for i=1:length(a)
%     if(a(i)==0 && A(i)==0 && c==0)
%         continue;
%     elseif(a(i)==0 && A
 
%Step 4-2s complement
OP=fliplr(OP);
for i=1:length(OP)
    if(OP(i)==1)
        i=i+1;
        while(i<=length(OP))
            OP(i)=~OP(i);
            i=i+1;
        end
        break;
    end
end
Q=fliplr(OP);
disp('Taking 2''s complement')
disp(Q)
 
 
 
%Step 6:divide by B
 
Q=fliplr(Q);
M=  B;
n=length(Q);
if(length(M)>length(Q))
    Q(numel(M))=0;
else
    if(length(M)<length(Q))
    M(numel(Q))=0;
    end 
end
Q=fliplr(Q);
M=fliplr(M);
A=zeros(1,length(Q));
I=n;
%J=[A,Q]
while(I~=0)
    J=[A,Q];
    J=[J(2:end),zeros(1,1)];
    A=J(1:length(Q));
    Q=J(length(Q)+1:end);
    MC=zeros(1,length(M));
    for i=length(M):-1:1
    if(M(i)==1)
        MC(i)=M(i);
        i=i-1;
        while(i>0)
            MC(i)=~M(i);
            i=i-1;
        end
        break;    
    end
    MC(i)=M(i);
    end
    MC;
    T=A;
    A=adds(MC,A);
    if(A(1)==1)
       Q(end)=0;
        A=T;
    else
        Q(end)=1;
    end
    I;
    A;
    Q;
    I=I-1;
end
disp('Quotient')
E= Q(find(Q==1,1):end);
disp(E)
disp('Remainder')
R= A(find(A==1,1):end);
disp(R)
disp('Encrypted output')
disp([E R])
 
 
 
function [op]=adds(x,y)
z(length(x)+1)=0;
for i=length(x):-1:1
    s(i)=xor(xor(x(i),y(i)),z(i+1));
    z(i)=or(or(and(z(i+1),x(i)),and(z(i+1),y(i))),and(x(i),y(i)));
end
c=z(1);
op=[s];
end
 
 
function [encod]=huffadapt()
clc
fid=fopen('seq.txt','r'); %Reading text input from .txt file
seq=fread(fid,'*char');
symbols=unique(seq); %Finding the possible source symbols
p=zeros(1,length(unique(seq)));
for i=1:length(symbols)
    for j=1:length(seq)
        if isequal(symbols(i),seq(j))
            p(i)=p(i)+1;
        end
    end
    p(i)=p(i)/length(seq); %Calculation of probability of ocurence of ith source symbol 
end
 
%Arranging probability in ascending order and arranging source symbols in
%corresponding order
for i=1:length(p)
    for j=i+1:length(p)
        if p(i)>p(j)
            temp=p(i);
            p(i)=p(j);
            p(j)=temp;
            
            temp1=symbols(i);
            symbols(i)=symbols(j);
            symbols(j)=temp1;
        end
    end
end
 
sym_num=[1:length(symbols)]; %Integer placeholders for each source symbol
input_str=zeros(1,length(seq)); 
for i=1:length(input_str)
    for j=1:length(symbols)
        if isequal(seq(i), symbols(j))
            input_str(i)=sym_num(j); %Storing each source symbol with its integer equivalent
        end
    end
end
disp('Input text sequence:');
fprintf(seq)
fprintf('\n')
disp('Unique characters:');
fprintf(symbols)
fprintf('\n')
disp('Corresponding numeric sequence:');
disp(sym_num)
disp('Respective probabilities:');
disp(p)
fprintf('\n')
 
%fprintf('\n')
%disp(input_str)
%fprintf('\n')
 
[dict,avglen]=huffmandict(sym_num,p);
temp=dict; 
t=dict(:,2);
for i=1:length(temp)
    temp{i,2}=num2str(temp{i,2});
end
disp('The huffman code dict:');
disp(temp) %Final Huffman Dictionary
 
%HUFFMAN ENCODING PROCESS WITH INBUILT FUNCTION
% encod=huffmanenco(input_str,dict);
% disp('The encoded output:');
% disp(encod);
 
%HUFFMAN ENCODING WITHOUT INBUILT FUNCTION
ptr=1;encod=[];
for i=1:length(input_str)
    bin_val=str2num(temp{input_str(i),2});
    for j=ptr:length(bin_val)+ptr-1
        encod(j)=bin_val(j-ptr+1);
    end
    ptr=ptr+length(bin_val);
end
disp('Encoded output:');
disp(encod)
 
%HUFFMAN DECODING PROCESS WITHOUT INBUILT FUNCTION
ptr=1;i=1;decod=[];index=1;flag=1;
while i<=length(temp)
    temp_check=str2num(temp{i,2});
    len=length(temp_check);
    flag=1;
    if ptr+len-1>length(encod)
        i=i+1;
        continue;
    end
    for j=ptr:ptr+len-1
        if encod(j)==temp_check(j-ptr+1)
            flag=1;
        else
            flag=0;
            break
        end
    end
    if flag==1
        decod(index)=symbols(i);
        index=index+1;
        ptr=ptr+len;
        i=1;
    else
        i=i+1;
    end
 
end
disp('Decoded sequence:');
disp(native2unicode(decod,'ASCII'))
end
 
%HUFFMAN DECODING WITH INBUILT FUNCTION
% % bits=input('Enter the bit stream in[];');
% decod=huffmandeco(encod,dict);
% disp('The symbols are:');
% disp(decod); 



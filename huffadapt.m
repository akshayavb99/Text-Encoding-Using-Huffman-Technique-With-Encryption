function huffadapt()
clc
fid=fopen('seq.txt','r');
seq=fread(fid,'*char');
symbols=unique(seq);
p=zeros(1,length(unique(seq)));
for i=1:length(symbols)
    for j=1:length(seq)
        if isequal(symbols(i),seq(j))
            p(i)=p(i)+1;
        end
    end
    p(i)=p(i)/length(seq);
end

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
sym_num=[1:length(symbols)];
input_str=zeros(1,length(seq));
for i=1:length(input_str)
    for j=1:length(symbols)
        if isequal(seq(i), symbols(j))
            input_str(i)=sym_num(j);
        end
    end
end
fprintf(seq)
fprintf('\n')
fprintf(symbols)
fprintf('\n')
disp(p)
fprintf('\n')
disp(sym_num)
fprintf('\n')
disp(input_str)
fprintf('\n')
[dict,avglen]=huffmandict(sym_num,p);
temp=dict;
t=dict(:,2);
for i=1:length(temp)
    temp{i,2}=num2str(temp{i,2});
end
disp('The huffman code dict:');
disp(temp)

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
disp(native2unicode(decod,'ASCII'))
end

%HUFFMAN DECODING WITH INBUILT FUNCTION
% % bits=input('Enter the bit stream in[];');
% decod=huffmandeco(encod,dict);
% disp('The symbols are:');
% disp(decod); 

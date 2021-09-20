clear all; close all; clc;
N = 203;
Nr= 100;
G=[1,1,1,1;1,0,0,1];
SNRdB = 1:1:9;        %Signal to Noise Ratio in dB
SNR = 10.^(SNRdB/10);  %Signal to Noise Ratio in Linear Scale
for k = 1:length(SNRdB) % Monte Carlo estimation of BER
    error_decoding=0;
    error_encoding=0;
    error_nocoding=0;
    for kk=1:Nr
        data = round( rand(1,200) );% input bits generation
        data = [data,zeros(1,3)];
        data_nocoding=data*2-1;
        data_encoding=coding(data);
        data_bpsk=data_encoding*2-1;
        noi = randn(1,406);
        noise = noi*sqrt(1/(2*SNR(k)));
        noise_encoding=noi*sqrt(1/(2*(SNR(k)-10*log10(2))));
        out_awgn = data_bpsk + noise;
        out_awgn_encoding= data_bpsk+noise_encoding;
        out_awgn_nocoding = data_nocoding + noise(1:N);
        out_awgn(out_awgn<0)=0;
        out_awgn(out_awgn>0)=1;
        out_awgn_nocoding(out_awgn_nocoding<0)=0;
        out_awgn_nocoding(out_awgn_nocoding>0)=1;
        out_awgn_encoding(out_awgn_encoding<0)=0;
        out_awgn_encoding(out_awgn_encoding>0)=1;
        out_decoding=viterbi2(out_awgn,G,2,1,3);
        error_decoding=error_decoding+(N-length(find(out_decoding == data)));
        error_encoding=error_encoding+(2*N-length(find(out_awgn_encoding == data_encoding)));
        error_nocoding=error_nocoding+(N-length(find(out_awgn_nocoding == data)));
    end
    BER_decoding(k) = error_decoding/(N*Nr);
    BER_encoding(k) = error_encoding/(2*N*Nr);
    BER_nocoding(k) = error_nocoding/(N*Nr);
end

semilogy(SNRdB,BER_decoding,'r', 'linewidth' ,1.0);
hold on
semilogy(SNRdB,BER_encoding,'b','linewidth',1.0);
hold on
semilogy(SNRdB,BER_nocoding,'g','linewidth',1.0); %Theoritical Bit Error Rate
hold on
line([1 9],[0.01 0.01],'linestyle','--', 'Color','r', 'LineWidth', 1);
hold on 
line([1 9],[0.001 0.001],'linestyle','--', 'Color','r', 'LineWidth', 1);
hold on 
line([1 9],[0.0001 0.0001],'linestyle','--', 'Color','r', 'LineWidth', 1);
title('BPSK over AWGN Simulation');xlabel('Eb/N0 in dB');ylabel('BER');
legend('BER(decoding)','BER(encoding)','BER(nocoding)')
axis tight
grid

function y=coding(x)
N = 203;
States = {   1   5     [0 0] [1 1]; ... state 000
    1   5     [1 1] [0 0]; ... state 001
    2   6     [1 0] [0 1]; ... state 010
    2   6     [0 1] [1 0]; ... state 011
    3   7     [1 0] [0 1]; ... state 100
    3   7     [0 1] [1 0]; ... state 101
    4   8     [0 0] [1 1]; ... state 110
    4   8     [1 1] [0 0]; ... state 111
    
    };
%               col 1:  next state if input = 0
%               col 2:  next state if input = 1
%               col 3:  output value if input = 0
%               col 4:  output value if input = 1

current_state = States(1,:);                % begins with state 1
for i = 1:N
    output(:,i) = current_state{x(i)+3}; % output under current state
    next_state = current_state{x(i)+1};  % index of next state
    current_state = States(next_state,:);   % update current state
end
out_encoding=[];
for i=1:N
    out_encoding=[out_encoding,output(1,i),output(2,i)];
end
y=out_encoding;
end
function y=viterbi2(x,G,n,k,m)
%x is received signal           
%G is generate matrix
%（n,k,m）convolutional code
a=size(x);
s=a(2)*k/n;%after decoding
r=zeros(1,s);  %final result
ra=zeros(2^m,s+1);%2^m routine
tempra=zeros(2^m,s+1);%minimum routine each time
Fa=zeros(2^m,1);
%transition matrix
g=size(G);
q=g(2)-1;
%transition matrix of b1,b2
T1=inf(2^m,2^m);
T2=inf(2^m,2^m);
for i=0:2^m-1
    z=ten2d(i,q);
    T1(i+1,floor((i+8)/2)+1)=mod(1+sum(z.*G(1,2:end)),2);
    T1(i+1,floor((i+0)/2)+1)=mod(0+sum(z.*G(1,2:end)),2);
    T2(i+1,floor((i+8)/2)+1)=mod(1+sum(z.*G(2,2:end)),2);
    T2(i+1,floor((i+0)/2)+1)=mod(0+sum(z.*G(2,2:end)),2);
end
temp=0;
%8 states：
%s0-s7:000,001,010,011,100,101,110,111
%ra is the states of each time（decimal）
%first state is s0
next_state=zeros(2^m,n);%the next probable state of each route
ra_temp=zeros(2,s+1);
for j=1:s
    if temp >=m
        temp = temp-1;
        %after the fourth comparison
        %delete the biggest four 
        for jj=1:4
            [~,db]=max(tempra(:,j));
            tempra(db,:)=[];
            ra(db,:)=[];
        end
        %after the fourth time
        for i=1:2^temp
            next_state(i,:)=find(T1(ra(i,j)+1,:)~=inf)-1;%next possible state（each one is to two probabilty）
            %calculate the first distance
            Fa(i,1)=dis(x(2*j-1),x(2*j),T1(ra(i,j)+1,next_state(i,1)+1),T2(ra(i,j)+1,next_state(i,1)+1));
            %the second one
            Fa(i,2)=dis(x(2*j-1),x(2*j),T1(ra(i,j)+1,next_state(i,2)+1),T2(ra(i,j)+1,next_state(i,2)+1));
            ra(i+2^temp,:)=ra(i,:);
            tempra(i+2^temp,:)=tempra(i,:);
            ra(i,j+1)=next_state(i,1);
            ra(i+2^temp,j+1)=next_state(i,2);
            %calculate the total distance
            tempra(i,j+1)=tempra(i,j)+Fa(i,1);
            tempra(i+2^temp,j+1)=tempra(i+2^temp,j)+Fa(i,2);
        end
        %after the fourth comparison
        temp=temp+1;
    else
        for i=1:2^temp
            
            %adding up the first four
            %current state：ra(i,j)，this probability's minimum routine :tempra(i,j)
            next_state(i,:)=find(T1(ra(i,j)+1,:)~=inf)-1;%next possible state（each one is to two probabilty）
            %calculate the first distance
            Fa(i,1)=dis(x(2*j-1),x(2*j),T1(ra(i,j)+1,next_state(i,1)+1),T2(ra(i,j)+1,next_state(i,1)+1));
            %calculate the second distance
            Fa(i,2)=dis(x(2*j-1),x(2*j),T1(ra(i,j)+1,next_state(i,2)+1),T2(ra(i,j)+1,next_state(i,2)+1));
            ra(i+2^temp,:)=ra(i,:);
            tempra(i+2^temp,:)=tempra(i,:);
            ra(i,j+1)=next_state(i,1);
            ra(i+2^temp,j+1)=next_state(i,2);
            %calculate the total distance
            tempra(i,j+1)=tempra(i,j)+Fa(i,1);
            tempra(i+2^temp,j+1)=tempra(i+2^temp,j)+Fa(i,2);
        end
        temp=temp+1;
    end
end
%select the best routine from the final routines
for jj=1:2^temp-1
    [~,db]=max(tempra(:,j+1));
    tempra(db,:)=[];
    ra(db,:)=[];
end
r=ra(2:end);%states after the first time，the first state is 0
for i=1:size(r,2)
    y(:,i)=ten2d(r(i),m);
end
y=y(1,:);
end
function y=ten2d(x,n)
%change decimal to binary
y=zeros(1,n);
for i=1:n
    y(i)=rem(x,2);
    x=floor(x/2);
end
y=fliplr(y);
end
function d=dis(m,n,g,h)
d=xor(m,g)+xor(n,h);
end
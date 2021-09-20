clear all; close all; clc;
N = 203;
Nr= 10000;
G=[1,1,1,1;1,0,0,1];
SNRdB = 1:1:9;        %Signal to Noise Ratio in dB
SNR = 10.^(SNRdB/10);  %Signal to Noise Ratio in Linear Scale
for k = 1:length(SNRdB) % Monte Carlo estimation of BER
    error_il=0;
    error_noil=0;
    for kk=1:Nr
        data = round( rand(1,200) );% input bits generation
        data = [data,zeros(1,3)];
        data_encode=coding(data);
        for i=1:14
            for ii=1:29
                data_encoding(ii,i)=data_encode(ii+29*(i-1));
            end
        end
        for h=1:29
            for hh=1:14
                data_encoding_il(hh+14*(h-1))=data_encoding(h,hh);
            end
        end
        data_bpsk=data_encoding_il*2-1;
        data_noil=data_encode*2-1;
        noi1 = randn(1,100);
        noi2 = randn(1,10);
        noi3 = randn(1,296);
        noise = [noi1*sqrt(1/(2*SNR(k))),noi2*sqrt(1/(2*SNR(k)*0.1)),noi3*sqrt(1/(2*SNR(k)))];
        out_awgn = data_bpsk + noise;
        out_noil = data_noil + noise;
        out_noil(out_noil<0)=0;
        out_noil(out_noil>0)=1;
        out_awgn(out_awgn<0)=0;
        out_awgn(out_awgn>0)=1;
        for m=1:29
            for mm=1:14
                out_awgn_il(m,mm)=out_awgn(mm+14*(m-1));
            end
        end
        for g=1:14
            for gg=1:29
                out_il(gg+29*(g-1))=out_awgn_il(gg,g);
            end
        end
        out_il=viterbi2(out_il,G,2,1,3);
        out_noil=viterbi2(out_noil,G,2,1,3);
        error_il=error_il+(N-length(find(out_il == data)));
        error_noil=error_noil+(N-length(find(out_noil == data)));
    end
    BER_il(k) = error_il/(N*Nr);
    BER_noil(k) = error_noil/(N*Nr);
end

semilogy(SNRdB,BER_il,'--r', 'linewidth' ,1.0);
hold on
semilogy(SNRdB,BER_noil,'b--','linewidth',1.0);
hold on
title('BPSK over AWGN Simulation');xlabel('Eb/N0 in dB');ylabel('BER');
legend('BER(il)','BER(noil)')
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
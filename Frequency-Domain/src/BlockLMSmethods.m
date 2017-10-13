clear all;

%% Nested

n=2^10*100; % u length
M=2^10; %filter coefficients
L=M; %block size=filter size
kmax=n/M; %blocks number
mu=6.0698e-04;%4.56e-04;
%white noise with s=0.57 and mean =0
v = sqrt(0.57)*randn(n,1); 
v = v - mean(v);  

u = zeros(n,1);
u(1) = v(1);
for i=2:n
  u(i) = -0.34 * u(i-1) + v(i);
end
d=plant(u')'; %given signal

e=zeros(n,1); %error
w=zeros(L,1); %coefficients
y=zeros(n,1); %output
J=zeros(kmax,1); %learning curve


tic; %timer
for k=2:kmax 
     %following BlockLMS.pdf page 5
     fi=zeros(L,1);
     for i=1:L 
         for j=1:M
             y((k-1)*L+i)=y((k-1)*L+i)+w(j)*u((k-1)*L+i-j+1);%output
         end         
         e((k-1)*L+i)=d((k-1)*L+i)-y((k-1)*L+i); %error
         fi=fi+mu*e((k-1)*L+i)*u((k-1)*L+i:-1:(k-2)*L+i+1); %phi
         J(k)=J(k)+e((k-1)*L+i)^2; %learning curve
     end
     w=w+fi; %filter
end
time=toc;
fprintf('Time for nested method : %0.4f  \n',time);
%ploting learning curves
z=ones(100,1)*(10^(-5));
figure
semilogy(J);
hold on;
plot(z);
xlabel('Iterations')
ylabel('Error')
title('Learning curves with method "nested loops"');
hold off;


%% arrays
mu=6.0698e-04;%3.56e-04;
%white noise with s=0.57 and mean =0

d=plant(u')'; %given signal

w=zeros(M,1); %filter coefficient
y=zeros(n,1); %output
e=zeros(n,1); %error
J=zeros(kmax,1); %learning curve

tic; %timer
for k=1:kmax-1 
    U=toeplitz(u(k*M:1:(k+1)*M-1),u(k*M:-1:(k-1)*M+1)); %toeplitz u
    D=d(k*M:1:(k+1)*M-1); %d vector
    Y=U*w; %output vector
    Err=D-Y; %erro vector
    e(k*M:1:(k+1)*M-1)=Err; %error
    fi=U.'*Err; 
    w=w+mu*fi; %filter
    J(k)=J(k)+sum(Err.^2);  %learning curve
end
time=toc;
fprintf('Time for arrays method : %0.4f  \n',time);

% learning curves 
z=ones(100,1)*(10^(-5));
figure
semilogy(J);
hold on;
plot(z);
xlabel('Iterations')
ylabel('Error')
title('Learning curves with array opperations');
hold off;

%% FFT

kmax=n/M; %blocks number
mu=6.0698e-04;%3.56e-04;
%white noise with s=0.57 and mean =0
v = sqrt(0.57)*randn(n,1); 
v = v - mean(v);  
u = zeros(n,1);
u(1) = v(1);
for i=2:n
  u(i) = -0.34 * u(i-1) + v(i);
end

d=plant(u')'; %given signal

a=0.5;  %step size
g=0.7;  %forgetting factor
d=d(:); %vectorizing d
u=u(:); %vectorizing u
e=d; %error
P=1*ones(2*M,1); %energy
w=zeros(2*M,1);  %filter
J=zeros(kmax,1); %learning curves
tic;
for k=1:kmax-1
    %following BlockLMS.pdf page 17
     U=fft([u((k-1)*M+1:(k+1)*M)],2*M); %2 last blocks fft
     Y=ifft(U.*w); %inverse fft in order to find output
     Y=Y(M+1:2*M,1); %last block
     D=d(k*M+1:(k+1)*M); %desired signal
     e(k*M+1:(k+1)*M,1)=D-Y; %error
     Err=fft([zeros(M,1);e(k*M+1:(k+1)*M)],2*M); %fft while padding
     P=g*P+(1-g)*abs(U).^2; 
     Dvector=1./P;
     fi=ifft(Dvector.*conj(U).*Err,2*M); 
     fi=fi(1:M);
     J(k)=J(k)+sum(real(D-Y).^2); %learning cuves
     w=w+a*fft([fi;zeros(M,1)],2*M); %filter
end
 e=real(e(:)); %hold only real values
 w=ifft(w); %inverse back in time region
 w=real(w(1:length(w)/2));
 time=toc;
 fprintf('Time for method using fft : %0.4f  \n',time);

%learning curves 
z=ones(100,1)*(10^(-5));
figure
semilogy(J);
hold on;
plot(z);
xlabel('Iterations')
ylabel('Error')
title('Learning curves with fft');
hold off;

%% Unconstrained


mu=3.7e-04;%3.56e-04;
k_max=n/M; %number of blocks
%white noise with s=0.57 and mean =0
v = sqrt(0.57)*randn(n,1); 
v = v - mean(v);  
u = zeros(n,1);
u(1) = v(1);
for i=2:n
  u(i) = -0.34 * u(i-1) + v(i);
end

 d=plant(u')'; %desired signal
 d=d(:); %vectorizing d
 u=u(:); %vectorizing u
 e=d; %error
 w=zeros(2*M,1); %filter
 J=zeros(k_max,1); %learning curves
 tic;
 for k=1:k_max-1 

    %following BlockLMS.pdf page 17 without striped box
     U=fft([u((k-1)*M+1:(k+1)*M)],2*M);%2 last blocks fft
     Y=ifft(U.*w);%inverse fft in order to find output
     Y=Y(M+1:2*M,1);%last block
     D=d(k*M+1:(k+1)*M); %desired signal
     e(k*M+1:(k+1)*M,1)=D-Y;%error
     Err=fft([zeros(M,1);e(k*M+1:(k+1)*M)],2*M);%fft while padding
     J(k)=J(k)+sum(real(D-Y).^2); %learning curves
     w=w+mu*conj(U).*Err; %filter
 end
 e=real(e(:));%real values only
 w=ifft(w);%back in time origin
 w=real(w(1:length(w)/2));
 time=toc;
 fprintf('Time for unconstrained method : %0.4f \n',time);
 
% learning curves 

z=ones(100,1)*(10^(-5));
figure
semilogy(J);
hold on;
plot(z);
xlabel('Iterations')
ylabel('Error')
title('Learning curves unconstrained method');
hold off;


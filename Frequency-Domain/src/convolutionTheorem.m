clear all
clc

n=2^3;

x = randn(n,1) + 1i*randn(n,1);
y = randn(n,1) + 1i*randn(n,1);
%% correct solution with matlabs built in function conv
ys=conv(x,y);

%% Y*x with Y being Toeplitz array of y
Y= toeplitz([y;zeros(length(x)-1,1)],[y(1);zeros(length(x)-1,1)]);
ytoep=Y*x;

%% C*x with C being Circulant array of Y
C=zeros(2*n-1,2*n-1);
C(1:end,1)=Y(:,1) ;

for j=2:2*n-1
    C(1,j)=C(end,j-1);
    C(2:end,j)=C(1:end-1,j-1);
end


yCirc=C*[x;zeros(2*n-1-length(x),1)];


%% FFT-IFFT
L=length(x)+length(y)-1;
yfft=ifft(fft(x,L).*(fft(y,L)));

e(1)=norm(ytoep -ys);
e(2)=norm(yCirc-ys);
e(3)=norm(yfft-ys);

fprintf('Error between Toeplitz method and conv(): %e\n',e(1));
fprintf('Error between Circulant method and conv(): %e\n',e(2));
fprintf('Error between FFT method and conv(): %e\n',e(3));
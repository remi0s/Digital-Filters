function Y = fft_recursive(x)
n = length(x); %input length
global countMuls; %number of muls
global countSums; %number of adds

if (n == 1)
    
  Y = conj(x);
  
else 
    
  Ye = fftrecursive(x(1:2:n)); %even values
  k = n/2;
  W = exp(2*pi*1i*[0:1:k-1]/n);
  Yo = fftrecursive(x(2:2:n)).*W; %%odd values
  
  % count the number of muls and sums
  if n==2
      countSums=countSums+2; %%T(2)=4
  else
      countMuls=countMuls+n/2;
      countSums=countSums+2*length(Yo);
  end
  
  Y = [Ye+Yo Ye-Yo];
  
end
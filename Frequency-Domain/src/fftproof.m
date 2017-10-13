% "Proof" by MATLAB
% A simple technique to develop and verify the steps of a proof 
% using random data input
%
% N P P
% Cornell U 
% Sept 1992
%

clear

n = 2^3; % any even

% input 
%generated once in order to avoid every time new data
x=[-0.082494425370955 + 0.303520794649354i;-1.933022917850987 - 0.600326562133734i;-0.438966153934773 + 0.489965321173948i;-1.794678841455123 + 0.739363123604474i;0.840375529753905 + 1.711887782981555i;-0.888032082329010 - 0.194123535758265i;0.100092833139322 - 2.138355269439939i;-0.544528929990548 - 0.839588747336614i];
%x = randn(n,1) + 1i*randn(n,1);
% correct answer
ys = fft(x);

% root of unity
w = @(n,e) exp(-2*pi*1i.*e/n);

k = (0:n-1)';

% DFT proof steps
y = zeros(n,1);

for j = 0:n-1
  y(j +1) = sum(w(n,j*k) .* x(k +1));
end

fprintf('DFT : %e\n', norm(y - ys))

% split output top bottom
y = zeros(n,1);

for j = 0:n/2-1
  y(j +1) = sum(w(n,j*k) .* x(k +1));
end
for j = n/2:n-1
  y(j +1) = sum(w(n,j*k) .* x(k +1));
end

fprintf('split output top bottom : %e\n', norm(y - ys))

% split input even odd
y = zeros(n,1);

k = (0:n/2-1)';
for j = 0:n/2-1
  y(j +1) = sum(w(n,j*2*k) .* x(2*k +1)) + sum(w(n,j*(2*k+1)) .* x(2*k+1 +1));
end
for j = n/2:n-1
  y(j +1) = sum(w(n,j*2*k) .* x(2*k +1)) + sum(w(n,j*(2*k+1)) .* x(2*k+1 +1));
end

fprintf('split input even odd : %e\n', norm(y - ys))

%%completing the proof with the last line
k = (0:n/2-1)';

for j=0:(n/2)-1
 fe(j+1) = sum(w(n/2,j*k) .* x(2*k+1));
 fo(j+1) = sum(w(n/2,j*k) .* x(2*k+1+1));
 
end

wfo = w(n,(0:n/2-1)') .* fo.'; 
y = [fe.' + wfo; fe.' - wfo];

fprintf('done : %e\n', norm(y - ys))


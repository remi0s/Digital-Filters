clear all
global countSums;
global countMuls;
countSums=0;
countMuls=0;
x=[-0.082494425370955 + 0.303520794649354i;-1.933022917850987 - 0.600326562133734i;-0.438966153934773 + 0.489965321173948i;-1.794678841455123 + 0.739363123604474i;0.840375529753905 + 1.711887782981555i;-0.888032082329010 - 0.194123535758265i;0.100092833139322 - 2.138355269439939i;-0.544528929990548 - 0.839588747336614i];

y=fftrecursive(x);
compl=2*countSums+6*countMuls; %complexity of my recursive fft
compl_2=T(length(x)); %complexity of fft as given
ys=fft(x); %correct result

if compl==compl_2
    fprintf('Both methods have same complexity\nRecursive Error =: %e\n', norm(y-ys'))
end
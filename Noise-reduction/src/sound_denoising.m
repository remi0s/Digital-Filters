clear all
close all
clc

load('sound.mat','d','fs');
load('noise.mat','u');
fprintf('\t\t<<Playing a sample of the sound given as input>>\n\n');
player = audioplayer(d(1600000:1850000), fs);
%play sound
playblocking(player)

n=length(d);

% Create u(n-1) and u(n-2) vectors
uDelay = [0 ; u(1 : n - 1)];
uDelay2 = [0 ; 0 ; u(1 : n - 2)];

%% Calculate auto-correlation R array
r0 = mean(u .* u);
r1 = mean(u .* uDelay); 
r2 = mean(u .* uDelay2);
r=[r0;r1;r2];
R = [r0 r1 r2; r1 r0 r1; r2 r1 r0]; 

fprintf('Auto-correlation R array\nR=');
fprintf([repmat('\t%0.4f\t', 1, size(R, 2)) '\n'], R');

p0 = mean(u .* d);
p1 = 0;%mean(uDelay .* d);
p2 = 0;%mean(uDelay2 .* d);
P=[p0 ; p1;p2];
% a = xcorr(u,u,3-1,'unbiased');
% a = a(3:(2*3-1));
% R = toeplitz(a);
% P=[0.72; 0 ; 0];
w0 = R \ P;

%% Calculate the range of the coefficient ì 
minlR = 0;
maxlR = 2 / max(eig(R));
fprintf('Range of coefficient ì : 0 < ì < %f \n', maxlR);


%% steepest descent method
w=[1;1;1];
w_old = [0;0;0];
m = 0.01;
while max(abs(w - w_old)) > 1.0e-7
    w_old=w;
    w = w + m*(P-R*w);
    
end
T = [u [0; u(1:n-1)] [0; 0; u(1:n-2)]]; 
y = T*w;
e = d - y;

fprintf('Stepest descent method Weiner coefficients\n w= ');
fprintf([repmat('%0.4f\t', 1, size(R, 2)) '\n'], w');

fprintf('\n\t\t<<Now playing the same sample after filtering>>\n');
player = audioplayer(e(1600000:1850000), fs);
%play filtered sample
playblocking(player)
fprintf('\t\t\t\t\tAll the noise is gone!\n');


%% Steepest Descent for many steps
mus = 0.01:0.001:(2/max(eig(R)));
iterations = zeros(1,length(mus));
for index = 1:length(mus)
    mu = mus(index);
    wst=[1;1;1];
    w_old = [0;0;0];
    i = 0;
    while max((abs(wst(:,end) - w_old)) > 1.0e-7) && (i < 5000)
        w_old = wst;
        wst = wst + mu*(P-R*wst);
        i = i + 1;
    end 
    iterations(index) = i;
end

% Plot
figure()
plot(mus,iterations,'LineWidth',1);
grid
title('Steepest Descent Convergence')
xlabel('Step $$\mu$$','interpreter','latex')
ylabel('Iterations')

sound(e,fs);
fprintf('\n\n\t\t\t\t<<Now playing the whole song>>\n');
fprintf('\t\t\tTo stop sound write command "clear sound"\n');



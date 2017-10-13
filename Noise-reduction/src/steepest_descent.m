clear all
close all
clc
%% declerations
n = 2000;
varA=0.15; %variance of random variable A
varV=0.32; %variance of white noise

%define the vectors
u = zeros(n,1); 
d = zeros(n,1);


%% Create variable A with variance 0.15 and mean 0
A=sqrt(0.15).*randn(n,1);
A=A-mean(A); %remove the mean value of A in order to mean(A)=0;
%fprintf('mean(A)=%0.6f var(A)= %0.4f \n',mean(A),var(A));

%% Create white noise with variance 0.32 and mean 0
v = sqrt(0.32) .* randn(n, 1);
v=v-mean(v);

%% Functions Creation

% Create the function x(n) of the signal as an vector of n steps
for j=1:n
    x(j) = A(j)* sin( (pi/8)*j + pi/6 );
end
x=x';

% Create u(n) and d(n) vectors
u(1) = v(1); %init the values that can not be measured.
u(2) = v(2);
for j = 3 : n %start from 3 because of u(i-2)
    u(j) = 0.25 * u(j-1)-0.12*u(j-2) + v(j); %u(n)
    d(j) = x(j) +v(j);
end


% Create u(n-1) and u(n-2) vectors
uDelay = [0 ; u(1 : n - 1)];
uDelay2 = [0 ; 0 ; u(1 : n - 2)];

%% Calculate auto-correlation R array
r0 = mean(u .* u);
r1 = mean(u .* uDelay); 
r2 = mean(u .* uDelay2);
r=[r0 ; r1; r2];
R = [r0 r1 r2; r1 r0 r1; r2 r1 r0]; %R for 3x3 


%% Calculate cross-correletaion vector P
p0 = mean(u .* d);
p1 = 0;%mean(uDelay .* d);
p2 = 0;%mean(uDelay2 .* d);
P=[p0 ; p1 ; p2];
 


%% Calculate the optimal Wiener coefficients
w0 = R \ P;


%% Calculate the range of the coefficient ì 
minlR = 0;
maxlR = 2 / max(eig(R));
fprintf('Range of coefficient ì : 0 < ì < %f \n', maxlR);
%% Calculate Theoretical minimum mean-square error
sigma_d=mean(d.*d);
J_min = sigma_d - P'*w0;
fprintf('Theoretical minimum mean-square error: %f\n', J_min);


%%Steepest descent method
wsm=[0.5;0.5;0.5];
w_old = [0;0;0];
m = 1;
while max(abs(wsm - w_old)) > 1.0e-9
    w_old=wsm;
    wsm = wsm + m*(P-R*wsm);

end
T = [u [0; u(1:n-1)] [0; 0; u(1:n-2)]]; 
y = T*wsm;
e = d - y;




fprintf('r=');
fprintf([repmat('\t%0.4f\t', 1, size(r, 2)) '\n'], r');
fprintf('\nR=');
fprintf([repmat('\t%0.4f\t', 1, size(R, 2)) '\n'], R');
fprintf('\nP=');
fprintf([repmat('\t%0.4f\t', 1, size(P, 2)) '\n'], P');
fprintf('\nw0=');
fprintf([repmat('\t%0.4f\t', 1, size(w0, 2)) '\n'], w0');

mus = [0.001 1 3 4 maxlR+1 maxlR+2];
for j = 1:length(mus)
    % Steepest Descent for sample steps
    mu = mus(j);
    w_old = [0 ; 0 ; 0];
    w = [0.5 ; 0.5 ; 0.5];
    iterations = 0;
  
    for z=1:n
        w_old = w(:,end);
        w = [w (w(:,end) + mu*(P-R*w(:,end)))];
        if (w(:,end) - w_old)~=0
        iterations = iterations + 1;
        end
    end 
    
     % Plots
    figure();
    subplot(2,1,1);
    hold on
    plot(w(1,1),w(2,1),'o','Markersize',10,'MarkerEdgeColor','k','MarkerFaceColor','b');
    plot(w(1,:),w(2,:),'LineWidth',1);
    plot(w0(1),w0(2),'o','Markersize',10,'MarkerEdgeColor','k','MarkerFaceColor','r');
    hold off
    grid
    axis([0 2 -1 1]);
    title(['Steepest Descent for $$\mu = ' num2str(mu) '$$'],'interpreter','latex')
    legend('Start Point',['Iterations = ' num2str(iterations)],'Correct solution')
    xlabel('$$w_0$$','interpreter','latex')
    ylabel('$$w_1$$','interpreter','latex')
    
    subplot(2,1,2);
    wnew = [0.5; 0.5 ;0.5];
    wt = zeros(3, n);
    for i = 1 : n
        wnew = wnew + mu * (P- (R * wnew));
        wt(:, i) = wnew;
    end
    we = (wt - w0 * ones(1, n)) .^ 2;
    e = sqrt(sum(we));
    semilogy(e);
    xlabel('n time steps','interpreter','latex')
    ylabel('Parameter error','interpreter','latex')
    title(['Parameter error for $$\mu = ' num2str(mu) '$$'],'interpreter','latex')
    %saveas(gcf,[num2str(j) '.png']);
end



% Steepest Descent for many steps
mus = 0.01:0.001:(2/max(eig(R)));
iterations = zeros(1,length(mus));
for index = 1:length(mus)
    mu = mus(index);
    w_old = [0 ; 0 ; 0];
    w = [0.5 ; 0.5 ; 0.5];
    j = 0;
    while max((abs(w(:,end) - w_old)) > 1.0e-9) && (j < n)
        w_old = w;
        w = w + mu*(P-R*w);
        j = j + 1;
    end 
    iterations(index) = j;
end

% Plot
figure()
plot(mus,iterations,'LineWidth',1);
grid
title('Steepest Descent Convergence')
xlabel('Step $$\mu$$','interpreter','latex')
ylabel('Iterations')
%saveas(gcf,'Convergence.png');




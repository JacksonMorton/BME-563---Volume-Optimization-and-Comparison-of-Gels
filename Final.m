% Final Project
 % Jackson Bruce Morton II
 % 04-30-14
 % BME 563
 
 % I have adhered to the Duke Community Standard in completing  
 % this assignment. - Jackson Bruce Morton II

% Define problem parameters 
m = [63.0 25.4 48.4 81.6 66.2 92.8 5.70]; % P*s^(n-1)
n = [0.455 0.569 0.518 0.309 0.512 0.450 0.618]; % unitless
T_0 = [2 0 0 20 2 38 0]; % Yield Stress
t = 0:1:7200; % seconds
V = 2:0.1:5; % mL
A_max = [90 100 110]; % cm^2

% Preallocations
A_t1 = zeros(length(V),length(m),length(A_max)); G_t1 = A_t1; A_t2 = A_t1; G_t2 = A_t1;
area = zeros(length(t),length(V),length(m),length(A_max)); h = area;
F_t1 = A_t1; F_t2 = A_t1; 
M_t1 = zeros(length(m),length(A_max)); M_t2 = M_t1; SF_t1_V = M_t1; SF_t2_V = M_t1;

% Complie matricies for A and G
for i=1:length(m)
    for j=1:length(V)
        for k=1:length(A_max)
            [A,~,~,G] = gel_area_final(t,m(i),n(i),T_0(i),V(j),A_max(k));
            A_t1(j,i,k) = A(3601); G_t1(j,i,k) = G(3601);
            A_t2(j,i,k) = A(end);  G_t2(j,i,k) = G(end);
            area(:,j,i,k) = A; h(:,j,i,k) = V(j)./A; G_total(:,j,i,k) = G;
        end
    end
end
 
%Compile matricies for F
for i=1:length(m)
    for j=1:length(V)
        for k=1:length(A_max)
            [F] = gel_spreading(t,area(:,j,i,k),h(:,j,i,k),V(j),A_max(k));
            F_t1(j,i,k) = F(3600); 
            F_t2(j,i,k) = F(end);
        end
    end
end

% Calculate the Scoring Function
SF_t1 = F_t1.*G_t1; SF_t2 = F_t2.*G_t2;

% Calculate "best volume" over the 6 combinations of A_max and time
average = zeros(length(V), length(m));
c_t1 = 2; c_t2 = 1;
for i=1:length(m)
a = SF_t1(:,i,1) + SF_t1(:,i,2) + SF_t1(:,i,3); t1_avg = a./length(A_max);
b = SF_t2(:,i,1) + SF_t2(:,i,2) + SF_t2(:,i,3); t2_avg = b./4;
average(:,i) = (c_t1.*t1_avg + c_t2.*t2_avg)./(c_t1+c_t2);
end
 SF_averaged = max(average);
 
for i=1:length(m)
[~,I]=max(average(:,i)); V_ideal(i) = V(I); 
    for k=1:length(A_max)    
        M_t1(i,k) = F_t1(I,i,k); M_t2(i,k) = F_t2(I,i,k);
        SF_t1_V(i,k) = SF_t1(I,i,k); SF_t2_V(i,k) = SF_t2(I,i,k);
    end
end

figure(2); clf
subplot(1,2,1)
plot(t(1:7200)/60,G_total(:,9,2,2), '-k')
xlabel('time (minutes)'); ylabel('{\itG_O(t)}');
title('New Method for Computing {\itG(t)} (x_1 = 0.95)')
subplot(1,2,2)
plot(area(1:7200,9,2,2),G_total(:,9,2,2), '-k')
xlabel('area (cm^2)'); ylabel('{\itG_O(A)}');
title('New Method for Computing G vs. Area (x_1 = 0.95)')



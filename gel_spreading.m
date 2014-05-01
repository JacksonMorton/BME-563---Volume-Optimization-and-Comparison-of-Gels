% Homework #5
 % Jackson Bruce Morton II
 % 4-14-14
 % BME 563
 
 % I have adhered to the Duke Community Standard in completing  
 % this assignment. - Jackson Bruce Morton II

function [F]=gel_spreading(t,A,h,V,A_max)

global D % A_max
D = 6.*10.^(-6); % cm^2/s
%A_max = 100; % cm^3

inner_int = cumtrapz(t, 1./h.^2);
for i = 1:length(t)
    sum = 0;
    for n = 0:100
        sum = sum + exp(-D.*(n*pi+0.5*pi).^2.*inner_int(i));
    end
    sums(i) = sum;
end
sums = sums(2:end); sums = sums';

% sums = [0; sums];
% y = (A./h).*(A<A_max) + (A_max./h).*(A>=A_max);
% M = cumtrapz(t, 2.*D.*y.*sums);

y = (A(2:end)./h(2:end)).*(A(2:end)<A_max) + (A_max./h(2:end)).*(A(2:end)>=A_max);
M = cumtrapz(t(2:end), 2.*D.*y.*sums);
F = M./V;

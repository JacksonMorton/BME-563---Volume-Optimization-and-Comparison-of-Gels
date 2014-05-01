% Homework #5
 % Jackson Bruce Morton II
 % 4-14-14
 % BME 563
 
 % I have adhered to the Duke Community Standard in completing  
 % this assignment. - Jackson Bruce Morton II

 function [A,h,t_L,G]= gel_area_final(t,m,n,T_0,V,A_max)
 
 global w h_0 %A_max
 w = 2; % cm
 h_0 = 0.1; % cm
 %A_max = 100; % cm^3
 
 %% Determine A(t) and h(t)
 F = 1.26.*V.*10.^4; % N = (kg*m)/s^2
 C = ((2.*w.*m)./(n+2)).*((2.*n+1)./n).^n.*(V./(4.*w)).^(n+2);
 Q = T_0.*(V.^2./(16.*w));
 
 hp=@(t,h)(-1).*(C.^(-1./n).*h.^((3.*n+2)./n).*(F.*h-Q.*h.^(-2)).^(1./n));
 [t,h]=ode45(hp,t,h_0);
 A = V./h;

 
 %% Determine t_L, the time at which the gel starts to leak
 t_L = 0;
 if max(A) <= A_max
     t_L = 7200;
 else
     t_L = t(find(A >= A_max, 1));  
 end
 
 %% Determine G(t)
 V_L = 1; % mL
 A_prime = A_max./(1-(V_L./V));
 
 q = A_max; r = 0.75.*A_max;
 s = 1;     u = 0.5;
 
 z = log(s./u)./(q-r);
 y = s./exp(z.*q);
 
 G = zeros(1,length(t));
 for i=1:length(t)
     if A(i) > A_prime || A(i) < A_max./2
        G(i) = 0;
     elseif A(i) <= A_max
            % G(i) = A(i)./A_max;
            G(i) = y.*exp(z.*A(i));
     else
         G(i) = (A_prime - A(i))./(A_prime - A_max);
     end   
 end
 
%  for i=1:length(t)
%      if A(i) <= A_max
%          G(i) = A(i)./A_max;
%      elseif A(i) > A_prime
%          G(i) = 0;
%      else
%          G(i) = (A_prime - A(i))./(A_prime - A_max);
%      end   
%  end
 
 
 G = G'; 
 G = G(2:end);

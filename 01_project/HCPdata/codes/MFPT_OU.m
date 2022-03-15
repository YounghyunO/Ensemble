%Daniel Charlebois - January 12 2011 - MATLAB v7.11 (R2010b)
%Calculate MFPT of an OU process.

%parameters
k = 100; %number of iterations
a = 2.0; %absorbing barrier (a>=0)
x0 = 0.0; %initial position (x0>=0)

sum_holder = 0; %initialize holder for the sum calculation
for i=1:k
   sum_holder = sum_holder + ((sqrt(2)*a)^i)/factorial(i)*gamma(i/2); %for x0=0, a~=0 
   %sum_holder = sum_holder + (-1)^(i+1)*((sqrt(2)*x0)^i)/factorial(i)*gamma(i/2); %for x0~=0, a=0
end
MFPT = 1/2*sum_holder

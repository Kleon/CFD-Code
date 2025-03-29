%Purpose: calculates integral from x1 to xn of y(x) dx
% using the composite trapezoidal integration rule
%for non equidistant intervals

%Input: 
% x: column vector of n non-equidistant increasing x-values
% y: column vector of n corresponding y-values

%Output: 
% i: integral evaludated using the composite trapezoidal integration rule
% for non equidistant intervals

%Global: 
%h: mesh spacing

function [I] = myTrapezoidal(x,y)
global h

%Preallocate variables: 
I = zeros(size(y));

%Calculate integral from x1 to xn utilizing trapezoidal integration: 
for i = 1:length(x)-1       %Loop over x
    
    h = abs(x(i+1) - x(i)) ;                  %Determine h
    I(i) = (h/2) *(y(i)+y(i+1)) ;       %Area of  Trapezoid
end

%Calculate sum:
I = sum((I));

end
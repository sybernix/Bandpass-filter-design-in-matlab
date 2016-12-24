function [ result ] = myBessel( x )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    k=1;
    result=0;
    term = 10;
    while (term>10^(-6))
        term = (((x/2)^k)/(factorial(k)))^2;
        result = result + term;
        k = k+1;
    end
    
    result = result+1;

end







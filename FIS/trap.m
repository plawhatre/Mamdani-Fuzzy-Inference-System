%% TREPEZOIDAL MEMBERSHIP FUNCTION 
% clc,clear,close all
function [ y ]=trap( x,a,b,c,d )
if x<=a
    y=0;
elseif a<=x & x<=b
    y=(x-a)/(b-a);
elseif b<=x & x<=c
    y=1;
elseif c<=x & x<=d
    y=(d-x)/(d-c);
else
    y=0;
end
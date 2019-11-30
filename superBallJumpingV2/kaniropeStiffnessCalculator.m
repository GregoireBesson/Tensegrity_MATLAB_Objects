% Kanirope Cable Stifness Calculator
% Greg Besson
% 29/08/2019

clc
clear all 
close all

% PE MONOSPAN:
% the tables contain the strength value in N of the cables
% at 25, 50, 75 and 100% of elongation
elongation = [25; 50; 75; 100];
strengthD6 = [25; 34; 42; 54];
strengthD8 = [41; 54; 66; 79];

% passive and active cables have different length
passiveLength = 0.18;
activeLength = 0.33;

stiffnessPassiveD6 = strengthD6./(passiveLength*elongation/100);
stiffnessActiveD6 = strengthD6./(activeLength*elongation/100);
stiffnessPassiveD8 = strengthD8./(passiveLength*elongation/100);
stiffnessActiveD8 = strengthD8./(activeLength*elongation/100);

clc; clear all; close all;    
% This file aims at testing two SNL examples:
%    Square network with 4 fixed anchors, i.e., pro=1 
%    EDM network with m random anchors,   i.e., pro=2  
% You can find more SNL examples to test in 'SNLExamples' file
% If you want to test more difficult problems such as nf=0.5, then try to
% increase range such as range = 0.4 and set pars as 
%    pars.update = 1; 
%    pars.rho    = log(n)/5;
% to make solver render more accurate solutions


pro = 1;
n   = 500;
nf  = 0.1; 

switch pro;   
    case 1
    problem     = 'BLTWY06_Inner';           
    range       = 0.2; 
    m           = 4;
    [D,dim,pars]= SNLExamples(problem,'multiplicative',n,m,nf,range);
    case 2                                  
    range       = 10;  
    m           = 20;
    [D,dim,pars]= EDMExamples(n,m,nf,range,1);   
end

%Call solver SQREDM
if nf>0.2;
pars.update= 1; 
pars.rho   = log(n)/5;
end
pars.draw  = 1;
Outsqredm  = SQREDM(D,dim,pars) ;
 

clc; clear all; close all;    
% This file aims at testing 2 MC problems with PDB data.
% The file 'setupMC' was adopted from Prof. Toh's package for solving MC
% problem in his paper:
%     "K.F. Jiang, D.F Sun and K.C. Toh, Solving nuclear norm regularized
%      and semidefinite matrix least squares problems with linear equality
%      constraints, Discrete Geometry and Optimization, 2013."

pro = 1;
if pro==1; load 1LFB.mat  % n=641
else       load 1RGS.mat  % n=2015      
end

%Generate the problem
Inputs.noiseType     = 'normal'; % uniform 
Inputs.sparseLevel   = 0.5;      % Note: use 0.3-1.0 to make problem easier
Inputs.noiseLevel    = 0.1; 
Inputs.radius        = 6;
Inputs.nDimensions   = 3;
Inputs.minLowerBound = 1.0;      % min distance for bounded atoms           
Inputs.randstate     = 0;             
[PP,D,Lmat,Umat]     = setupMC(porig,Inputs);            
%Set up parameters
[dim,n]              = size(PP);
pars.LOWBD           = full(Lmat);
pars.UPPBD           = full(Umat);
pars.PP              = PP;
D                    = sparse(D);

%Call solver SQREDM 
pars.draw  = 1; 
Out        = SQREDM(D,dim,pars); 

 



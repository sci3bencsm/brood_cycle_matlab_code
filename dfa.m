function [U,V,eigenvals] = DFA(X,group,maxfac)
%[U,V,eigenvals] = DFA(X,group,maxfac)
% Performs DISCRIMINANT FUNCTION ANALYSIS
%
% INPUT VARIABLES
%
% X              = data matrix that contains m groups
%                  Dim(X) = [N x M]. All columns must be independent.
% group          = a vector containing a number corresponding
%                  to a group for every row in X. If you have 
%                  m groups there will be numbers in the range 
%                  1:m in this vector.
% maxfac         = the maximum number of DFA factors extracted
%
% OUTPUT VARIABLES
%
% U              = DFA scores matrix (Dim(U) = [N x maxfac])
%                  the eigenvalues are multiplied with each column
%                  in this matrix.
% V              = DFA loadings matrix, Dim(V) = [M x maxfac]
% eigenvals      = a vector of DFA eigenvalues
%
%
% Copyright, B.K. Alsberg, 1996

% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation version 3 of the License.
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.



[T,W]=TW_gen(X,group);

B = T-W;

invW = inv(W);
P = invW*B;

[vec1,val]=eig(P);
d=(diag(val))';
eigenvals = d(1:maxfac);

% Here we sort the eigenvectors w.r.t. the eigenvalues:
[dummy,idx]=sort(-eigenvals);
vec = vec1(:,idx);
eigenvals = eigenvals(idx);

%% V is the matrix of canonical variates directions %%%
V = vec(:,1:maxfac);

%% U is the matrix of scores %%%
U = X*V;

% new line to multiply eigenvalues to DFA directions:

U = U*diag(eigenvals);
% we need a reference to see how this is done properly


U = real(U);


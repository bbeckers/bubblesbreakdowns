function [ytrend,ycycle]=sparse_hp_filter(y,lambda)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%sparse_hp_filter.m
%
%This program executes a standard HP filter using a sparse matrix
%
%       Input:  y - a Txn data matrix, where T is the number of observations
%                       on n variables  (i.e., data is assumed to be in column
%                       format). Furthermore, it is assumed that T>4
%                   lambda - a scalar. This is a smoothing parameter. Optional: if not entered,
%                        a default value of 1600 will be used.
%
%       Output: ytrend - a (T-discard)xn matrix of extracted trends for
%                       each of the n variables.
%                   ycycle a  (T-discard)xn matrix of deviations from the extracted trends for
%                       each of the n variables. Optional.
%
% The HP filter finds a series {ytrend_t}_{t=1}^T for
% each n that solves the following minimization probelm
% 
%min sum_{t=1}^T(y_t-ytrend_t)
%      +lambda*sum_{t=2}^{T-1}[(ytrend_{t+1}-ytrend_{t})-(ytrend_{t}-ytrend_{t-1})]
%
%Rearanging the first order conditions yields
%  A*ytrend=y
%  where
% A=[   1+lambda    -2*lambda       lambda        0               0           0 ...
%     [    -2*lambda   1+5*lambda   -4*lambda       lambda       0           0 ...
%     [    lambda       -4*lambda        1+6*lambda   -4*lambda  lambda  0 ...
%     [                                         ...
%     [    0      ...   0   lambda       -4*lambda        1+6*lambda   -4*lambda  lambda 0  ... 0]
%     [                                         ...
%     [  0   ...    lambda       -4*lambda        1+6*lambda  -4*lambda  lambda]
%     [  0   ...           0     lambda  -4*lambda       1+5*lambda      -2*lambda ]
%     [  0   ...           0       0        lambda              -2*lambda         1+lambda]
%
%     See also hp_filter.m in Harald Uhlig's toolkit.
%
%   Hodrick, R.J. and E.C.Prescott (1997), "Postwar U.S. Business Cycles:
%   An Empirical Investigation." Jounal of Money, Credit and Banking.
%   29(1), Feb. pp. 1--16.
%
%Copyright: Alexander Meyer-Gohde
%
%You are free to use/modify/redistribute this program so long as original
%authorship credit is given and you in no way impinge on its free
%distribution
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2,  lambda = 1600; end %If the user didn't provide a value for lambda, set it to the default value 1600
[T,n] = size (y);% Calculate the number of periods and the number of variables in the series



%Preliminary calculations
x1=[1+lambda, -2*lambda, lambda];  %The non-zero elements of the first row of A
x2=[-2*lambda, 1+5*lambda, -4*lambda, lambda]; %The non-zero elements of the second row of A
x3=[lambda, -4*lambda, 1+6*lambda, -4*lambda, lambda]; %The non-zero elements of thej'th row of A, 2<j<T-2
x2rev=x2(end:-1:1); %The non-zero elements of the second-to-last row of A (just x2 in reverse)
x1rev=x1(end:-1:1); %The non-zero elements of the last row of A (just x1 in reverse)

%make a list (a column vector) containing at position i the row of the i'th non-zero element in A 
I=3:T-2;
I = I(ones(1,5),:);
I=I(:);
I=[1; 1; 1; 2; 2; 2; 2;I;T-1;T-1;T-1;T-1;T;T;T];

%make a list (a column vector) containing at position i the column of the i'th non-zero element in A 
J=1:T-4;
J = J(ones(1,5),:);
J=J(:);
Temp=0:4;
Temp=Temp(ones(1,T-4),:).';
Temp=Temp(:);
J=J+Temp;
J=[1;2;3;1;2;3;4;J;[T-3:1:T]';[T-2:1:T]'];

%make a list (a column vector) containing at position i the i'th non-zero element in A 
X=x3(ones(1,T-4),:).';
X=X(:);
X=[x1';x2';X;x2rev';x1rev'];
 
%Build the matrix A and solve the system A*ytrend=y
ytrend=sparse(I,J,X)\y;


if nargout==2 %Should the user have requested a second output
    ycycle=y-ytrend; %The second output will be the deviations from the HP trend
end

end
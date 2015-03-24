%% Unconditional VAR forecasts 
%  from the reduced-form model with recursively for h periods

%% INPUT:
% y [T,K]: data set
% p [integer]: VAR order (# of lags)
% A [K,K*p+1]: matrix of estimated VAR coeff.
% h [integer]: forecast horizon

%% OUTPUT:
% y_fc [K,h]: forecasted VAR variables

function y_fc = VARfc(y,p,A,h,c)

if nargin<5
    c = 1;
else c = 0;
end

K = size(y,2);

y_fc = zeros(K,h);
infoset = y(end-p+1:end,:);         % Obtain latest p observation
infoset = infoset(p:-1:1,:);        % Reorder so that newest are on top!
infoset = reshape(infoset',[],1);	% Vectorize
for j=1:h
    if c==1;
        y_fc(:,j) = A(:,1)+A(:,2:end)*infoset;  % Else: no seasonsal dummy
    else y_fc(:,j) = A*infoset;
    end
    infoset = [y_fc(:,j);infoset(1:end-K)];  
end

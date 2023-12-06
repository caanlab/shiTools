function B = shiStatLeastAbsDevReg(Y, X)
%
% performs least absolute deviations regression using iterative reweighted least-squares method
%
% Y - n x k outcomes
% X - n x p predictors
% B - p x k coefficients
%
% zhenhao shi
%
tolerance = 1e-6;
max_iter = 1e5;
delt = 1e-6;

[n, p] = size(X);
[~, k] = size(Y);

XX = kron( speye(k) ,         X ); % kn x kp
WW = kron( speye(k) ,  speye(n) ); % kn x kn
YY = kron( speye(k) , ones(n,1) ); % kn x k
for i = 1:k
    YY(YY(:,i)~=0,i) = Y(:,i);
end
BB = ( XX' * WW * XX ) \ ( XX' * WW * YY ); % kp * k

XXYYrow = kron( (1:k)' , ones(n,1) );  % [ ones(n,1); 2*ones(n,1); ...; k*ones(n,1) ]
BBrow = kron( (1:k)' , ones(p,1) );  % [ ones(p,1); 2*ones(p,1); ...; k*ones(p,1) ]

XXcol = kron( 1:k , ones(1,p) );  % [ ones(1,p), 2*ones(1,p), ..., k*ones(1,p) ]
YYBBcol  = 1:k;                      % [ 1, 2, ..., k ]

B = nan(p,k);

tol = tolerance + 1;

W0ind = kron( speye(k) , ones(n,1) ) > 0; % same size as YY
BBind = kron( speye(k) , ones(p,1) ) > 0; % save size as BB
k0 = k;

% fprintf('%10d',k0);
for i = 1:max_iter
    if max(tol) < tolerance 
        break;
    end
    BB_prev = BB;
    W0 = 1 ./ max(delt, abs(YY - XX * BB));
    WW = spdiags(W0(W0ind),0,k0*n,k0*n);
    BB = ( XX' * WW * XX ) \ ( XX' * WW * YY );
    tol = sum(abs(BB - BB_prev));
    if any(tol<tolerance)
        yind = YYBBcol(tol<tolerance);
        xbb = BB(ismember(BBrow,yind),ismember(YYBBcol,yind));
        xbbind = BBind(ismember(BBrow,yind),ismember(YYBBcol,yind));
        B(:,yind) = reshape(xbb(xbbind),p,numel(yind));
        XX = XX(~ismember(XXYYrow,yind),~ismember(XXcol,yind));
        YY    =    YY(~ismember(XXYYrow,yind),~ismember(YYBBcol,yind));
        W0ind = W0ind(~ismember(XXYYrow,yind),~ismember(YYBBcol,yind));
        BB    =    BB(~ismember(BBrow,yind),~ismember(YYBBcol,yind));
        BBind = BBind(~ismember(BBrow,yind),~ismember(YYBBcol,yind));
        XXYYrow( ismember( XXYYrow, yind ) ) = [];
        BBrow(   ismember( BBrow,   yind ) ) = [];
        XXcol(   ismember( XXcol,   yind ) ) = [];
        YYBBcol( ismember( YYBBcol, yind ) ) = [];
        k0 = k0 - numel(yind);
%         fprintf('\b\b\b\b\b\b\b\b\b\b%10d',k0);
    end
end

% if any(isnan(B(:)))
%     yind = find(any(isnan(B)));
%     xbb = BB(ismember(BBrow,yind),ismember(YYBBcol,yind));
%     xbbind = BBind(ismember(BBrow,yind),ismember(YYBBcol,yind));
%     B(:,yind) = reshape(xbb(xbbind),p,numel(yind));
% end

if max(tol) > tolerance
    warning('IRLS did not converge for %d Y column(s)!',full(sum(tol(:)>tolerance)));
end

fprintf('\n');
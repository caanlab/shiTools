function [r, F, df1, df2, p] = shiStatIcc( M, type, r0 )

% calculates intraclass correlation
%
% Intraclass correlation
% [r, F, df1, df2, p] = shiStatIcc( M, type, r0 )
%
% M is matrix of observations. Each row is an object of measurement and
% each column is a judge or measurement.
%
% 'type' is a string that can be one of the six possible codes for the desired
% type of ICC:
% '1 - 1': The degree of absolute agreement among measurements made on
% randomly seleted objects. It estimates the correlation of any two
% measurements.
% '1 - k': The degree of absolute agreement of measurements that are
% averages of k independent measurements on randomly selected
% objects.
% 'C - 1': case 2: The degree of consistency among measurements. Also known
% as norm - referenced reliability and as Winer's adjustment for
% anchor points. case 3: The degree of consistency among measurements maded under
% the fixed levels of the column factor. This ICC estimates the
% corrlation of any two measurements, but when interaction is
% present, it underestimates reliability.
% 'C - k': case 2: The degree of consistency for measurements that are
% averages of k independent measurements on randomly selected
% onbjectgs. Known as Cronbach's alpha in psychometrics. case 3:
% The degree of consistency for averages of k independent
% measures made under the fixed levels of column factor.
% 'A - 1': case 2: The degree of absolute agreement among measurements. Also
% known as criterion - referenced reliability. case 3: The absolute
% agreement of measurements made under the fixed levels of the column factor.
% 'A - k': case 2: The degree of absolute agreement for measurements that are
% averages of k independent measurements on randomly selected objects.
% case 3: he degree of absolute agreement for measurements that are
% based on k independent measurements maded under the fixed levels
% of the column factor.
%
% ICC is the estimated intraclass correlation. LB and UB are upper
% and lower bounds of the ICC with alpha level of signif icance.
%
% In addition to estimation of ICC, a hypothesis test is performed
% with the null hypothesis that ICC = r0. The F value, degrees of
% freedom and the corresponding p - value of the this test are
% reported.
%
% ( c ) Arash Salarian, 2008
%
% Reference: McGraw, K. O., Wong, S. P., "Forming Inferences About
% Some Intraclass Correlation Coefficients", Psychological Methods,
% Vol. 1, No. 1, pp. 30 - 46, 1996
%


if ~exist( 'r0', 'var' )||isempty( r0 )
    r0 = 0;
end


switch type
    case {'1-1','1','1,1'}
        [r, F, df1, df2, p] = ICC_case_1_1( M, r0 );
    case {'1-k','k','1,k'}
        [r, F, df1, df2, p] = ICC_case_1_k( M, r0 );
    case {'C-1','C,1','3,1'}
        [r, F, df1, df2, p] = ICC_case_C_1( M, r0 );
    case {'C-k','C,k','3,k'}
        [r, F, df1, df2, p] = ICC_case_C_k( M, r0 );
    case {'A-1','A,1','2,1'}
        [r, F, df1, df2, p] = ICC_case_A_1( M, r0 );
    case {'A-k','A,k','2,k'}
        [r, F, df1, df2, p] = ICC_case_A_k( M, r0 );
end


%%
function [r, F, df1, df2, p] = ICC_case_1_1( M, r0 )
[n, k] = size( M );
MSR = var( mean( M, 2 ) ) * k;
MSW = sum( var( M, 0, 2 ) ) / n;

r = ( MSR - MSW ) / ( MSR + ( k - 1 ) * MSW );
F = ( MSR / MSW ) * ( 1 - r0 ) / ( 1 + ( k - 1 ) * r0 );
df1 = n - 1;
df2 = n * ( k - 1 );
p = 1 - fcdf( F, df1, df2 );


%%
function [r, F, df1, df2, p] = ICC_case_1_k( M, r0 )
[n, k] = size( M );
MSR = var( mean( M, 2 ) ) * k;
MSW = sum( var( M, 0, 2 ) ) / n;

r = ( MSR - MSW ) / MSR;
F = ( MSR / MSW ) * ( 1 - r0 );
df1 = n - 1;
df2 = n * ( k - 1 );
p = 1 - fcdf( F, df1, df2 );


%%
function [r, F, df1, df2, p] = ICC_case_C_1( M, r0 )
[n, k] = size( M );
SStotal = var( M( : ) ) * ( n * k - 1 );
MSR = var( mean( M, 2 ) ) * k;
MSC = var( mean( M, 1 ) ) * n;
MSE = ( SStotal - MSR * ( n - 1 ) - MSC * ( k - 1 ) ) / ( ( n - 1 ) * ( k - 1 ) );

r = ( MSR - MSE ) / ( MSR + ( k - 1 ) * MSE );
F = ( MSR / MSE ) * ( 1 - r0 ) / ( 1 + ( k - 1 ) * r0 );
df1 = n - 1;
df2 = ( n - 1 ) * ( k - 1 );
p = 1 - fcdf( F, df1, df2 );


%%
function [r, F, df1, df2, p] = ICC_case_C_k( M, r0 )
[n, k] = size( M );
SStotal = var( M( : ) ) * ( n * k - 1 );
MSR = var( mean( M, 2 ) ) * k;
MSC = var( mean( M, 1 ) ) * n;
MSE = ( SStotal - MSR * ( n - 1 ) - MSC * ( k - 1 ) ) / ( ( n - 1 ) * ( k - 1 ) );

r = ( MSR - MSE ) / MSR;
F = ( MSR / MSE ) * ( 1 - r0 );
df1 = n - 1;
df2 = ( n - 1 ) * ( k - 1 );
p = 1 - fcdf( F, df1, df2 );


%%
function [r, F, df1, df2, p] = ICC_case_A_1( M, r0 )
[n, k] = size( M );
SStotal = var( M( : ) ) * ( n * k - 1 );
MSR = var( mean( M, 2 ) ) * k;
MSC = var( mean( M, 1 ) ) * n;
MSE = ( SStotal - MSR * ( n - 1 ) - MSC * ( k - 1 ) ) / ( ( n - 1 ) * ( k - 1 ) );

r = ( MSR - MSE ) / ( MSR + ( k - 1 ) * MSE + k * ( MSC - MSE ) / n );
a = ( k * r0 ) / ( n * ( 1 - r0 ) );
b = 1 + ( k * r0 * ( n - 1 ) ) / ( n * ( 1 - r0 ) );
F = MSR / ( a * MSC + b * MSE );
a = k * r / ( n * ( 1 - r ) );
b = 1 + k * r * ( n - 1 ) / ( n * ( 1 - r ) );
v = ( a * MSC + b * MSE )^2 / ( ( a * MSC )^2 / ( k - 1 ) + ( b * MSE )^2 / ( ( n - 1 ) * ( k - 1 ) ) );
df1 = n - 1;
df2 = v;
p = 1 - fcdf( F, df1, df2 );


%%
function [r, F, df1, df2, p] = ICC_case_A_k( M, r0 )
[n, k] = size( M );
SStotal = var( M( : ) ) * ( n * k - 1 );
MSR = var( mean( M, 2 ) ) * k;
MSC = var( mean( M, 1 ) ) * n;
MSE = ( SStotal - MSR * ( n - 1 ) - MSC * ( k - 1 ) ) / ( ( n - 1 ) * ( k - 1 ) );

r = ( MSR - MSE ) / ( MSR + ( MSC - MSE ) / n );
c = r0 / ( n * ( 1 - r0 ) );
d = 1 + ( r0 * ( n - 1 ) ) / ( n * ( 1 - r0 ) );
F = MSR / ( c * MSC + d * MSE );
a = k * r / ( n * ( 1 - r ) );
b = 1 + k * r * ( n - 1 ) / ( n * ( 1 - r ) );
v = ( a * MSC + b * MSE )^2 / ( ( a * MSC )^2 / ( k - 1 ) + ( b * MSE )^2 / ( ( n - 1 ) * ( k - 1 ) ) );
df1 = n - 1;
df2 = v;
p = 1 - fcdf( F, df1, df2 );



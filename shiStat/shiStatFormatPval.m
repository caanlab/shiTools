function str = shiStatFormatPval(p)

% formats p value and returns a string (e.g., p=0.01234 -> '=0.012')

if p>1
    error('p must be <1');
elseif p>=0.995
    str = '>0.99';
elseif p>0.0995
    str = sprintf('=%.2f',p);
elseif p>=0.001
    str = sprintf('=%.3f',p);
elseif p>=0
    str = '<0.001';
else
    error('p must be non-negative');
end
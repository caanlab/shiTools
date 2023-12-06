function x_split = shiStatSplit(x,cutoff,dim)

% splits a variable into two halves based on some cutoff criterion
%
% x_split = shiStatSplit(x,cutoff,dim)
%    split x into -1 vs 1 based on a cutoff value (<cutoff vs. >=cutoff)
%    
%    x      - array or matrix
%    cutoff - cutoff value: 'mean' (default),'median',or value (e.g. 3.14)
%    dim    - dimension along which to split
%
% zhenhao shi 2018 4 26

if isequal(x,[]), x_split = []; return; end

if nargin < 2
    cutoff = 'mean';
end
if nargin < 3
    % Figure out which dimension to work along.
    dim = find(size(x) ~= 1, 1);
    if isempty(dim), dim = 1; end
end

if ischar(cutoff)

    if strcmpi(cutoff,'mean')
        M = nanmean(x,dim);
    elseif strcmpi(cutoff(1:3),'med')
        M = nanmedian(x,dim);
    else
        error('Method must be ''mean'', ''median'', or a real value');
    end

    x_split = bsxfun(@minus,x, M);
    x_split = sign(x_split);

else

    x_split = sign(x-cutoff);

end

x_split(x_split==0) = 1;


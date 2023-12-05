function varargout = shiDeal(M)

% distributes each column of a matrix to a new variable
%

if ~ismatrix(M)
    error('input must be at most 2-D');
end

if nargout < size(M,2)
    warning('only %d of the %d columns will be distributed',nargout,size(M,2));
elseif nargout > size(M,2)
    error('too many output variables');
end

for i = 1:nargout
    varargout{i} = M(:,i);
end
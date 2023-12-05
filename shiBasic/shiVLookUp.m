function Return_value = shiVLookUp(Lookup_value,Table_array,Col_index_num,defNA)

% looks for values in the leftmost column of a table, and returns values in the same rows from a specified column, similar to MS Excel function vlookup. look-up stops when first match found
%
% zhenhao shi 2018-5-11

if ~exist('defNA','var')
    if isnumeric(Table_array)
        defNA = NaN;
        Return_value = [];
    elseif iscell(Table_array)
        defNA = {[]};
        Return_value = {};
    else
        error('should be either numeric matrix or cell array');
    end
end

for i = 1:numel(Lookup_value)
    xInd = NaN;
    for j = 1:size(Table_array,1)
        if isequaln(Lookup_value(i),Table_array(j,1))
            xInd = j;
            break;
        end
    end
    if isnan(xInd)
        Return_value(i,:) = defNA;
    else
        Return_value(i,:) = Table_array(xInd,Col_index_num);
    end
end

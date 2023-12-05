function [UDS_sorted,TimeCollect] = shiUds(UDS,TimeCollect_raw,TimeTarget_LowerBound,TimeTarget_UpperBound)

% organizes UDS results according to time stamps
%
% zhenhao shi 2018-6-20


dTime_1 = bsxfun(@minus,TimeCollect_raw(:)',TimeTarget_LowerBound(:));
dTime_2 = bsxfun(@minus,TimeCollect_raw(:)',TimeTarget_UpperBound(:));
Ind = dTime_1>=0 & dTime_2<=0;

UDS_sorted = nan(size(TimeTarget_LowerBound));
TimeCollect = cell(size(TimeTarget_LowerBound));
for i = 1:length(TimeTarget_LowerBound)
    if ~any(Ind(i,:))
        continue;
    end
    UDS_sorted(i) = nanmean(UDS(Ind(i,:)));
    TimeCollect{i} = TimeCollect_raw(Ind(i,:));
end;

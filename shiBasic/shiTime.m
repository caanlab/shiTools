function time = shiTime(Index)

% returns current date/time as a string
%
% time = shiTime(Index)
%
% Example1 :
%   time1 = shiTime(1);
%   time2 = shiTime(2);
%   time3 = shiTime(3);
%   time4 = shiTime(4);
%   time5 = shiTime(5);
%   time6 = shiTime(6);
%   {time1;time2;time3;time4;time5;time6}
%       = {'2014';'12';'24';'20';'53';'57'};
%       
% Example2 :
%   time1 = shiTime(1:3);
%   time2 = shiTime(4:6);
%   {time1,time2}
%       = {'20141224','205357'};
% 
% Example3 :
%   time = shiTime;
%   time
%       = '20141224205357'
% 
%    ###########
% by Zhenhao Shi @ 2014-12-24
%    ###########

if nargin == 0
    Index = 1:6;
end
    
a = clock;
D{1,1} = num2str(a(1),'%.4d');
D{1,2} = num2str(a(2),'%.2d');
D{1,3} = num2str(a(3),'%.2d');
D{1,4} = num2str(a(4),'%.2d');
D{1,5} = num2str(a(5),'%.2d');
D{1,6} = num2str(round(a(6)),'%.2d');

if ischar(Index)
    time = sprintf('%s%s%s%s%s%s%s%s%s%s%s',D{1},Index,D{2},Index,D{3},Index,D{4},Index,D{5},Index,D{6});
elseif min(Index) >= 1 && max(Index) <= 6
    time = cell2mat(D(Index));
else
    time = sprintf('%s-%s-%s %s:%s:%s',D{1},D{2},D{3},D{4},D{5},D{6});
end





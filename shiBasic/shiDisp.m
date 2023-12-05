function S=shiDisp(s,red)

% displays cellstring one line after another
%
% shiDisp(s,red)
%
%    ###########
% by Zhenhao Shi @ 2020-06-28
%    ###########

if nargin < 2
    red = false;
end

if red
    fid = 2;
else
    fid = 1;
end

s = char(s);
[I,J] = size(s);
s = [repmat('#',1,J);repmat(' ',1,J);s;repmat(' ',1,J);repmat('#',1,J)];
s = [repmat('#',I+4,1),['#';repmat(' ',I+2,1);'#'],s,['#';repmat(' ',I+2,1);'#'],repmat('#',I+4,1)];
if nargout == 0
    for i = 1:size(s,1)
        fprintf(fid,'%s\n',s(i,:));
    %     disp(s(i,:));
    end
    return;
else
    S = s;
end
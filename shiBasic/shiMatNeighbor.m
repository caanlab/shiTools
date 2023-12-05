function LATTICE = shiMatNeighbor(Size,Conn)

% returns lattice graph with edges linking neighboring matrix elements (voxels)
%
% [LATTICE] = shiMatNeighbor(Size,Conn)
%   Size - size of matrix
%   Conn - connectivity (6,18(default),26) to define neighbor
%   LATTICE - a sparce adjancency matrix, with neighbors connected
%
% Zhenhao Shi, 2020-05-12
%

switch numel(Size)
    case 1
        Size = [Size;1;1]';
        if Conn~=2 && Conn~=6 && Conn~=18 && Conn~=26
            warning('Conn = 2')
        end
        Conn = 6;
    case 2
        Size = [Size(:);1]';
        if Conn==4
            Conn = 6;
        elseif Conn==8
            Conn = 18;
        elseif Conn~=6 && Conn~=18 && Conn~=26
            warning('Conn = 8')
            Conn = 18;
        end
    case 3
        if Conn~=6 && Conn~=18 && Conn~=26
            warning('Conn = 18')
            Conn = 18;
        end
    otherwise
        error('only supports <=3 dimensions');
end

sizI = Size(1);
sizJ = Size(2);
sizK = Size(3);
cntM = sizI * sizJ * sizK;

[x,y,z]=meshgrid(-1:1,-1:1,-1:1);
GRID = [x(:),y(:),z(:)];
GRID_AbsSum = sum(abs(GRID),2);

switch Conn
    case 6
        GRID = GRID(GRID_AbsSum==1,:);
    case 18
        GRID = GRID(GRID_AbsSum==1|GRID_AbsSum==2,:);
    case 26
        GRID = GRID(GRID_AbsSum==1|GRID_AbsSum==2|GRID_AbsSum==3,:);
end

% xI = GRID(:,1);
% xJ = GRID(:,2);
% xK = GRID(:,3);
% 
% [I,J,K] = ind2sub([sizI,sizJ,sizK],1:cntM);
% 
% II = bsxfun(@plus,I,xI);
% JJ = bsxfun(@plus,J,xJ);
% KK = bsxfun(@plus,K,xK);
% 
% II(II>sizI | II<=0) = NaN;
% JJ(JJ>sizJ | JJ<=0) = NaN;
% KK(KK>sizK | KK<=0) = NaN;
% 
% Ind = sub2ind([sizI,sizJ,sizK],II,JJ,KK);
% 
% NEI = cell(Size);
% 
% for i = 1:cntM
%     xind = Ind(~isnan(Ind(:,i)),i);
%     NEI{i} = sparse(xind,repmat(i,size(xind)),true,cntM,cntM);
% end

M = reshape(1:cntM,[sizI,sizJ,sizK]);
xM = nan([Size,Conn]);

for c = 1:Conn
    xM(:,:,:,c) = shiMatShift(M,NaN,GRID(c,:));
end

xM1 = reshape(xM,cntM,Conn)';
xM0 = repmat(1:cntM,Conn,1);

LATTICE = [xM1(:),xM0(:)];
LATTICE = LATTICE(~any(isnan(LATTICE),2),:);
LATTICE = sparse(LATTICE(:,1),LATTICE(:,2),true,cntM,cntM);

% Mx = repmat(M,1,1,1,Conn);
% 
% % NEI = cell(Size);
% LATTICE = cell(Size);
% for i = 1:cntM
%     xind = xM(Mx==i);
%     xind = xind(~isnan(xind));
%     LATTICE{i} = [xind,repmat(i,size(xind))];
% %     NEI{i} = sparse([xind;repmat(i,size(xind))],[repmat(i,size(xind));xind],true,cntM,cntM);
% end
% 
% LATTICE = cat(1,LATTICE{:});
% LATTICE = sparse(LATTICE(:,1),LATTICE(:,2),true,cntM,cntM);



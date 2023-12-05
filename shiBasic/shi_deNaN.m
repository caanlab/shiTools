function [Img,Vector,AnyMiss] = shi_deNaN(Img,Vector)

if nargin < 2
    Vector = Img;
end

AnyMiss = any(~isfinite(Vector(:)));

RowIndex = ~any(~isfinite(Vector),2);
Img = Img(RowIndex,:);
Vector = Vector(RowIndex,:);

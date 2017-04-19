function [ x,y ] = TriuIdx2Coord( idx)
%TRIUIDX2COORD Summary of this function goes here
%   covernts indexed upper triangular vector to the corresponding X-Y 
%   coordinates in the matrix
%   The index vector represents thr row indices of upper triangular
%   elements for a given matrix. It will be converted in to the
%   corresponding coordiantes.
%   idx: upper triangular indices in column order
%   n: size of the matrix

y=ceil(0.5+sqrt(1+8.*idx)./2-1)+1;
x=idx-(y-1).*(y-2)./2;



end


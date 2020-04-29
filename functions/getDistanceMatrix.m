function distMat = getDistanceMatrix(posVert,posHoriz)
%GETDISTANCEMATRIX Vectorized method to obtain a 2D distance matrix between
%points in N dimensional cartesian coordinated. posVert should be VxN
%matrix, posHoriz a HxN matrix for any N (typically 2 or 3). Returns a VxH
%matrix distMat.

u = reshape(posVert,[size(posVert,1) 1 size(posVert,2)]);
b = reshape(posHoriz,[1 size(posHoriz,1) size(posHoriz,2)]);
% (n_pos_vert)x(n_pos_horiz) matrix containing distances
distMat = sqrt(sum( bsxfun(@minus,u,b).^2, 3));

end
% splits FOV into parts if scanning with high-resolution (>1024 pixels)

function [xFOVs, yFOVs] = jkget_xyFOVs(ops)

% use images only the ones specified as in ops.useX and ops.useY

Ly = ops.Ly;
Lx = ops.Lx;
% iR(1) = size of Y pixel subsets, iR(2) = size of X pixel subsets
iR = ops.splitFOV;

ny = floor(Ly/iR(1));
nx = floor(Lx/iR(2));

% indices for each subset
xFOVs = zeros(nx, iR(2), iR(1));
yFOVs = zeros(ny, iR(2), iR(1));
for i = 1:iR(1)
    for j = 1:iR(2)
        xFOVs(:,j,i) = min(ops.useX) + [1:nx] + (j-1)*nx;
        yFOVs(:,j,i) = min(ops.useY) + [1:ny] + (i-1)*ny;
    end
end

xFOVs = xFOVs(:,:);
yFOVs = yFOVs(:,:);

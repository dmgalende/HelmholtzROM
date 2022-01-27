function X = mesh1Dto2D(x,y)

nx = size(x,1); 
ny = size(y,1);
ini = 1;
fin = nx;
X = zeros(nx*ny,2);
for i = 1:ny
    X(ini:fin,:) = [x y(i)*ones(size(x))];
    ini = fin + 1;
    fin = fin + nx;
end


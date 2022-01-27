function T = createConnectivity_linTri(n,m)

nelems = 2*(n-1)*(m-1);
T = zeros(nelems, 3);

for j = 1:m-1
    for i = 1:n-1
        e1 = i + 2*(j-1)*(n-1);
        e2 = e1 + n - 1;
        v1 = i + (j-1)*n;
        v2 = v1 + 1;
        v3 = i + 1 + j*n;
        v4 = v3 - 1;
        T(e1,1) = v3;
        T(e1,2) = v1;
        T(e1,3) = v2;
        T(e2,1) = v1;
        T(e2,2) = v3;
        T(e2,3) = v4;
    end
end



function f = pointSource1D_TRI(X,T,s,referenceElement)

%Ruled connectivity around the element
f = 2;
h = max( abs( X(T(1,1),:) - X(T(1,2),:) ) );
limits = [ s(1)-f*h, s(1)+f*h, s(2)-f*h, s(2)+f*h ];
pos = createRuledConnectivity(X,T,limits);

%Search for element
elemsvec = 1:size(T,1);
elems = elemsvec(pos);
truelems = zeros(size(elems));
maxelems = length(elems);
in = false;
on = false;
i = 1;
while (~in || on) && (i <= maxelems) %%arreglar!!
    ielem = elems(i);
    [in,on] = inpolygon( s(1), s(2), X(T(ielem,:),1), X(T(ielem,:),2) );
    truelems( i(in) ) = ielem;
    i = i + 1;
end
truelems = truelems( logical(truelems) );

%Vandermonde matrix in the reference coordinates
V = Vandermonde_LP(referenceElement.degree, referenceElement.NodesCoord);
[L,U,P] = lu(V');

f = zeros(size(X,1),1);
for elem = truelems
    %Local coordinates from vertex coordinates (linear mapping)
    vc = X( T(elem,1:3), : );
    xieta = inverseLinearMapping(vc, s);

    %Pointsource term: N_i(s) for i = 1:size(X,1)
    p = orthopoly2D(xieta,referenceElement.degree);
    N = U \ ( L \ (P*p) );
    f(T(elem,:)) = f(T(elem,:)) + N;
end

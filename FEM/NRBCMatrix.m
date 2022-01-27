function C = NRBCMatrix(X,TRobin,theReferenceElement,k)

% Input:
%  X: nodal coordinates
%  TRobin: connectivity matrix for 1D boundary elements
%  theReferenceElement: information of the reference element
%  k: nodal values of wave number
%  ccg: nodal values of the c*cg parameter
%  R: radius of circular non-reflecting boundary (radiation)
% Output:
%  C: non reflecting boundary condition matrix

%Number of elements and number of mesh nodes
nOfElements = size(TRobin,1);
nOfNodes = size(X,1);
nOfBoundaryNodes = length(unique(TRobin));
nOfLinearNodes = nOfElements + 1;
nOfElementNodes = theReferenceElement.degree + 1;
nOfInnerNodes = nOfBoundaryNodes - nOfLinearNodes;

%Memory allocation
NNZ = nOfElementNodes*(2*nOfLinearNodes + nOfInnerNodes);
% C = spalloc(nOfNodes,nOfNodes,NNZ);
aux_ones = ones(1,nOfElementNodes);
I = zeros(NNZ,1);
J = I;
Creal = I;
Cimag = I;
 
%Loop in 1D boundary elements
for ielem = 1:nOfElements
 Te = TRobin(ielem,:);
 Xe = X(Te,:); % Coordinates of the current element nodes
 ke = k(Te);
 Ce = ElementalNRBCMatrix(Xe,theReferenceElement,ke);
 
 % Assembling
 Te_transp = transpose(Te);
 aux_row = Te_transp(:,aux_ones);
 aux_col = Te(aux_ones,:);
 indK = (ielem-1)*nOfElementNodes^2+1:ielem*nOfElementNodes^2;
 I(indK) = aux_row(:);
 J(indK) = aux_col(:);
 Creal(indK) = real(Ce);
 Cimag(indK) = imag(Ce);
 
 clear Ce;
end

C = sparse(I(I~=0),J(I~=0),Creal(I~=0) + sqrt(-1)*Cimag(I~=0),nOfNodes,nOfNodes);

%______________________________________________________________
function Ce = ElementalNRBCMatrix(Xe,theReferenceElement,ke)

nOfElementNodes = size(Xe,1);
Ce = zeros(nOfElementNodes,nOfElementNodes);

%Information of the reference element
IPw = theReferenceElement.IPweights1d; 
N = theReferenceElement.N1d; 
Nxi = theReferenceElement.N1dxi;

%Number of Gauss points
ngauss = length(IPw); 

%Loop in integration points
for g = 1:ngauss
  %Shape functions and derivatives at the current integration point 
  N_g = N(g,:);
  Nxi_g = Nxi(g,:);
  %Integration weight
  xyDer_g = Nxi_g*Xe;
  xyDerNorm_g = norm(xyDer_g);
  dline=IPw(g)*xyDerNorm_g;
  %wavenumber at current gauss point
  k = N_g*ke;
  %Contribution of the current integration point to the elemental matrix
  Ce = Ce + k*(N_g')*N_g*dline;
end
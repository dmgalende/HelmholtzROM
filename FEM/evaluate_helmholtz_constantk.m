function sol = evaluate_helmholtz_constantk(data, mesh, sample)

nf = length(sample.freq);
ns = length(sample.source_pos);
ndof = size(mesh.X, 1);
sol = zeros(ndof, nf, ns);

% Volume matrices
onesv = ones(ndof, 1);
[K, varM] = PGDberkhoffVolumeMatrices(...
    mesh.X,...
    mesh.T,...
    mesh.referenceElement,...
    onesv,... %for Kint
    onesv,... %for Kxpml
    onesv,... %for Kypml
    onesv,... %for varMint, varMpml
    []);
K = K{1};
varM = varM{1};

% Boundary matrices
BC = NRBCMatrix(mesh.X, mesh.Tb, mesh.referenceElement, onesv);

for i = 1:nf
    
    if data.verbose, disp(['Freq = ', num2str(sample.freq(i))]), end
    
    % Wavenumber and signal amplitude
    wn = 2 * pi * sample.freq(i) ./ mesh.velocity(1);
    amp = getSignalAmplituteFreq(sample.freq(i));    

    % System matrix and decomposition
    Mat = (wn^2 * varM + sqrt(-1) * wn * BC - K);
    [L, U, P, Q] = lu(Mat);
    
    % Source vector
    for j = 1:ns
        
        if data.verbose, disp(['   Pos = ', num2str(sample.source_pos(j))]), end
        
        % Point source located at the input parameter
        ys = mesh.ysourcepos;
        xs = mesh.xlim(1) + sample.source_pos(j) * (mesh.xlim(2) - mesh.xlim(1));
        f = pointSource1D_TRI(mesh.X, mesh.T, [xs, ys], mesh.referenceElement);
        
        % Solution
        sol(:,i,j) = Q * (U \ (L \ (P * (amp * f))));
    end
end








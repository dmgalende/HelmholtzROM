function mesh = build_mesh(data)

% Test choice definition
switch data.test
    case 1 % homogeneous (using mean heterogeneous velocity)
        
        % Domain and velocity
        vel_data = load('vel_z6.25m_x12.5m_exact_CROP.mat');
        model = mean(vel_data.C0(:));
        lx = vel_data.lx;
        lz = vel_data.lz;
        ax = 0; bx = lx;
        ay = 0; by = lz;
        homogeneous_flag = true;
        ysourcepos = by;
        
    case 2 % heterogeneous
        
        % Domain and velocity
        vel_data = load('vel_z6.25m_x12.5m_exact_CROP.mat');
        model = flipud(vel_data.C0);
        nx0 = vel_data.nx0;
        nz0 = vel_data.nz0;
        dx0 = vel_data.dx0;
        lx = vel_data.lx;
        lz = vel_data.lz;
        ax = 0; bx = lx;
        ay = 0; by = lz;
        homogeneous_flag = false;
        ysourcepos = by;
end

% Mesh properties
wavelength_min = min(model(:)) / data.max_freq;
elemSize = wavelength_min * data.elem_deg / data.wave_res;
nOfrefNodes = (data.elem_deg + 1) * (data.elem_deg + 2) / 2;
nx = ceil((bx-ax) / elemSize);
ny = ceil((by-ay) / elemSize);

% Nodal positions and mesh connectivities
xx = linspace(ax, bx, nx);
yy = linspace(ay, by, ny);
T = createConnectivity_linTri(nx, ny);
X = mesh1Dto2D(xx', yy');
[X, T] = increaseOrderFromLinearMesh(X, T, nOfrefNodes);
referenceElement = createReferenceElement(1, nOfrefNodes, []);

% Boundaries
tol = elemSize / (2 * data.elem_deg);
numNodes = 1:size(X,1);
posl = find(X(:,1) < ax+tol & X(:,1) > ax-tol); [~, sl] = sort(X(posl,2), 'descend');
posr = find(X(:,1) < bx+tol & X(:,1) > bx-tol); [~, sr] = sort(X(posr,2), 'ascend');
posd = find(X(:,2) < ay+tol & X(:,2) > ay-tol); [~, sd] = sort(X(posd,1), 'ascend');
posu = find(X(:,2) < by+tol & X(:,2) > by-tol); [~, su] = sort(X(posu,1), 'descend');
nodesl = numNodes(posl); nodesls = nodesl(sl); Tl = create1Dconec(ny-1, data.elem_deg, nodesls);
nodesr = numNodes(posr); nodesrs = nodesr(sr); Tr = create1Dconec(ny-1, data.elem_deg, nodesrs);
nodesd = numNodes(posd); nodesds = nodesd(sd); Td = create1Dconec(nx-1, data.elem_deg, nodesds);
nodesu = numNodes(posu); nodesus = nodesu(su); Tu = create1Dconec(nx-1, data.elem_deg, nodesus);
Tb = [Tl ; Td ; Tr ; Tu];

% Interpolated model
if homogeneous_flag
    mesh.velocity = model(1) * ones(size(X, 1), 1);
else
    xx0 = (ax:nx0-1) * dx0;
    yy0 = (ay:nz0-1) * dx0;
    X0 = mesh1Dto2D(xx0', yy0');
    model_t = model.';
    scatterer = scatteredInterpolant(X0, model_t(:));
    mesh.velocity = scatterer(X);
end

% Mesh output
mesh.X = X;
mesh.T = T;
mesh.Tb = Tb;
mesh.referenceElement = referenceElement;
mesh.xlim = [ax, bx];
mesh.ylim = [ay, by];
mesh.homogeneous_flag = homogeneous_flag;
mesh.ysourcepos = ysourcepos;





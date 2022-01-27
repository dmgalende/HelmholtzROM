
% Continuous finite element solution u(x; w,s) of the Helmholtz equation in 2D at frequency w and relative 
% source position s, with first-order non-absorbing boundary conditions and point source.
%
%           Div(Grad(u)) + k^2 * u = f  in computational domain
%          n x Grad(u) - i * k * u = 0  on boundaries
%
% where k = k(x; w) and f = f(x; w,s).
%
% Problem outline:
%                                        point source
%                                             |             
%                                             |
%                     x = a + s * (b-a)       |
%              <--------------------------->  |
%            x = a                           \ /   x = b
%              -------------------------------------
%             |                                     |
%             |                                     |
%             |         computational domain        |
%             |                                     |
%             |                                     |
%              -------------------------------------

clearvars 
close all
home

setpath()

%% TEST CONFIG

% Model design (computational domain, wave velocity and mesh)
data.test = 1;      % 1: homogeneous, 2: heterogeneous
data.wave_res = 8;  % Minimum wave resolution: points per (minimum) wavelength
data.elem_deg = 3;  % Order of the finite element approximation
data.max_freq = -1; % Maximum linear frequency in Hz that determines the spatial discretization. 
                    % Use -1 to get it from frequency samples
data.verbose = 1;   % Boolean value that prints sample info
 
% Evaluation sample
sample.freq = [1, 3, 5];             % Linear frequencies in Hz
sample.source_pos = [0.1, 0.5, 0.7]; % Relative source positions in [0,1] along the 1D surface of the domain

%% BUILD MESH

if data.max_freq == -1, data.max_freq = max(sample.freq); end % get maximum frequency from sample input
mesh = build_mesh(data);

%% EVALUATE HELMHOLTZ SOLUTION AT SAMPLES

sol = evaluate_helmholtz(data, mesh, sample); % Stored as (number_nodes, number_freq, number_sources) array

%% PLOTS

maxfigs = 10;

% Model and mesh
figure, hold on
sol_handle = plotSolution(mesh.X, mesh.T, mesh.velocity(:), mesh.referenceElement);
mesh_handle = plotMesh(mesh.X, mesh.T, mesh.referenceElement.faceNodes);
set(mesh_handle, 'FaceAlpha', 0, 'EdgeAlpha', 0.1)
set(sol_handle(2), 'Location', 'eastoutside')
title('MODEL VELOCITY')

% Solutions
[~, nf, ns] = size(sol);
counter = 1;
for i = 1:nf
    for j = 1:ns
        if counter == maxfigs, break, end
        sampletext = ['(', num2str(sample.freq(i)) ', ', num2str(sample.source_pos(j)), ')'];
        figure
        subplot(1, 2, 1), title(['Real solution at sample (freq, pos) = ', sampletext])
        sol_handle = plotSolution(mesh.X, mesh.T, real(sol(:,i,j)), mesh.referenceElement);
        set(sol_handle(2), 'Location', 'eastoutside')
        subplot(1, 2, 2), title(['Wave height at sample (freq, pos) = ', sampletext])
        sol_handle = plotSolution(mesh.X, mesh.T, abs(sol(:,i,j)), mesh.referenceElement);
        set(sol_handle(2), 'Location', 'eastoutside')
        counter = counter + 1;
    end
    if counter == maxfigs, break, end
end


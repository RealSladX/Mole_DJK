% ================ Problem 3 (Traverso et al. 2013) =================
% Heterogeneous, anisotropic, full-tensor conductivity (Darcy flow)
%
% PDE (potential u):
%   -div( C(x,y) * grad(u) ) = f(x,y),   (x,y) in (0,1)^2
%
% Full-tensor conductivity (Eq. 55):
%   C11 = y^2 + alpha*x^2
%   C22 = x^2 + alpha*y^2
%   C12 = C21 = (alpha-1)*x*y
%
% Analytic solution (Eq. 56):
%   u(x,y) = exp( -20*pi * ((x-1/2)^2 + (y-1/2)^2) )
%
% Boundary conditions:
%   Dirichlet on all boundaries: u = u_exact
%
% Discretization:
%   Mimetic Finite Differences using the SDSU MOLE mimetic operators.
%   Assemble: A = - D * K * G, where
%     G : mimetic gradient, mapping cell-centers u -> faces
%     D : mimetic divergence, mapping faces -> centers
%     K : block tensor conductivity at faces
%
% =====================================================================
close all; clc;

addpath('../../src/matlab_octave');

% ---- grid ----
k = 4;            % mimetic order (typical: 2, 4)
m = 64;            % cells in x
n = 64;            % cells in y
dx = 1/m;
dy = 1/n;


% cell centers (include boundary nodes)
xc = [0 dx/2: dx :1-dx/2 1]';
yc = [0 dy/2: dy :1-dy/2 1]';
[Yc,Xc] = meshgrid(yc, xc);

% interior cell centers 
xc_int = [dx/2: dx :1-dx/2]';
yc_int = [dy/2: dy :1-dy/2]';

% faces
xf = (0:dx:1)';
yf = (0:dy:1)';
% [Yf,Xf] = meshgrid(yf, xf);
[Yf_x, Xf_x] = meshgrid(yc_int, xf);    % (m+1) x n  → x-faces
[Yf_y, Xf_y] = meshgrid(yf, xc_int);    % m x (n+1)  → y-faces


%Exact Solution (Calculated at the centers)
a = -20 * pi;
z = (Xc - 0.5).^2 + (Yc - 0.5).^2;
ue = exp(a * z);


% Boundary Conditions (Dirichlet)
dc = [1;1;1;1];
nc = [0;0;0;0];
bcl = squeeze(ue(1,:))'; % left bc (y increases)
bcr = squeeze(ue(end,:))'; % right bc (y increases)
bcb = squeeze(ue(:,1)); % bottom bc (x increases)
bct = squeeze(ue(:,end)); % top bc (x increases)
bcl = bcl(2:end-1,1);
bcr = bcr(2:end-1,1);
v = {bcl;bcr;bcb;bct};

alpha = 100; % 1, 10, or 100

% Gradient Operator
G = grad2D(k, m, dx, n, dy);
% Divergence Operator
D = div2D(k, m, dx, n, dy);

% Building tensor components at Faces
C11F = Yf_x.^2 + alpha*Xf_x.^2;
% C12F_x = (alpha - 1).*Xf.*Yf;
% C12F_y = (alpha - 1).*Xf.*Yf;
C12F_x = (alpha - 1).*Xf_x.*Yf_x;
C12F_y = (alpha - 1).*Xf_y.*Yf_y;
% C12F_x = zeros(size(C12F_x));
% C12F_y = zeros(size(C12F_y));
C22F = Xf_y.^2 + alpha*Yf_y.^2;

% Get dimension of Gradient Operator
rowsGx = (m+1)*n;
rowsGy = m*(n+1);

Nx = rowsGx;
Ny = rowsGy;
Icf = interpolCentersToFacesD2D(k, m, n);
Ifc = interpolFacesToCentersG2D(k, m, n);
N  = (m+2)*(n+2);
Cx2x = Icf(1:Nx, 1:N);                 % center scalar -> x-faces
Cy2y = Icf(Nx+1:Nx+Ny, N+1:2*N);       % center scalar -> y-faces
Fx2c = Ifc(1:N, 1:Nx);                 % x-faces -> center scalar
Fy2c = Ifc(N+1:2*N, Nx+1:Nx+Ny);       % y-faces -> center scalar
Iy2x = Cx2x * Fy2c;   % y-faces -> centers -> x-faces
Ix2y = Cy2y * Fx2c;   % x-faces -> centers -> y-faces

% Iy2x = buildYFacesToXFaces(m, n);
% Ix2y = buildXFacesToYFaces(m, n);
% Full Heterogeneous Anisotropic 2x2 Tensor at Faces
K11 = spdiags(reshape(C11F, [], 1), 0, rowsGx, rowsGx);
K12 = spdiags(reshape(C12F_x, [], 1), 0, rowsGx, rowsGx) * Iy2x;
K21 = spdiags(reshape(C12F_y, [], 1), 0, rowsGy, rowsGy) * Ix2y;
% K12 = spdiags(reshape(C12F_x, [], 1), 0, rowsGx, rowsGy);
% K21 = spdiags(reshape(C12F_y, [], 1), 0, rowsGy, rowsGx);
% K12 = sparse(rowsGx, rowsGy);
% K21 = sparse(rowsGy, rowsGx);
% K21 = K12';
K22 = spdiags(reshape(C22F, [], 1), 0, rowsGy, rowsGy);

K = [K11 K12;
     K21 K22];

% Analytical Derivatives (For RHS)
ux = a * (2*Xc -1) .* ue;
uy = a * (2*Yc -1) .* ue;
uxx = (2*a + a^2*(2*Xc - 1).^2) .* ue;
uyy = (2*a + a^2*(2*Yc - 1).^2) .* ue;
uxy = a^2 * (2*Xc - 1) .* (2*Yc - 1) .* ue;

% Building tensor components at Centers (For RHS)
C11 = Yc.^2 + alpha*Xc.^2;
C1221 = (alpha - 1).*Xc.*Yc;
C22 = Xc.^2 + alpha*Yc.^2;

% Tensor component derivatives at centers (For RHS)
C11_x = 2*alpha*Xc;
C11_y = 2*Yc;

C22_x = 2*Xc;
C22_y = 2*alpha*Yc;

C12_x = (alpha-1).*Yc;
C12_y = (alpha-1).*Xc;
% C12_x = zeros(size(C12_x));
% C12_y = zeros(size(C12_y));

% Build RHS
dxdot = C11_x.*ux + C11.*uxx + C12_x.*uy + C1221.*uxy;
dydot = C12_y.*ux + C1221.*uxy + C22_y.*uy + C22.*uyy;
RHS = -(dxdot + dydot);
F = reshape(RHS, [], 1); 

% Build System of Equations
L = -D*K*G;
[L0, F0] = addScalarBC2D(L, F, k, m, dx, n, dy, dc, nc, v);
ua = L0\F0;
ua = reshape(ua, m+2, n+2);

err = ua - ue;
fprintf('max abs err = %e\n', max(abs(err(:))));
fprintf('rel L2 err   = %e\n', norm(err(:))/norm(ue(:)));

% PLOTTING
figure(69);
contour3(Xc, Yc, ue);
title(sprintf("Exact Solution (alpha = %.2f)", alpha));
shading interp;
view([0 90]);
colorbar;


figure(420);
contour3(Xc, Yc, ua);
title(sprintf("Approximate Solution (alpha = %.2f)", alpha));
shading interp;
view([0 90]);
colorbar;


figure(1);
contour3(Xc, Yc, err);
title(sprintf('Error (alpha = %.2f)', alpha));
view([0 90]);
colorbar;


% function Ix2y = buildXFacesToYFaces(m, n)
% % Map x-face field ((m+1) x n) -> y-face field (m x (n+1))
% % Uniform-grid, node-averaging transfer.
% %
% % Vectorization convention:
% %   x-faces stored as reshape(Vx, [], 1) with size(Vx) = [m+1, n]
% %   y-faces stored as reshape(Vy, [], 1) with size(Vy) = [m, n+1]
%
%     rowsGx = (m+1)*n;
%     rowsGy = m*(n+1);
%
%     I = [];
%     J = [];
%     S = [];
%
%     for i = 1:m          % y-face x-index (centered in x)
%         for j = 1:(n+1)  % y-face y-index (nodal in y)
%
%             row = sub2ind([m, n+1], i, j);
%
%             if j == 1
%                 % bottom boundary: use first x-face column only
%                 col1 = sub2ind([m+1, n], i,   1);
%                 col2 = sub2ind([m+1, n], i+1, 1);
%
%                 I = [I; row; row];
%                 J = [J; col1; col2];
%                 S = [S; 0.5; 0.5];
%
%             elseif j == n+1
%                 % top boundary: use last x-face column only
%                 col1 = sub2ind([m+1, n], i,   n);
%                 col2 = sub2ind([m+1, n], i+1, n);
%
%                 I = [I; row; row];
%                 J = [J; col1; col2];
%                 S = [S; 0.5; 0.5];
%
%             else
%                 % interior: average four surrounding x-faces
%                 col1 = sub2ind([m+1, n], i,   j-1);
%                 col2 = sub2ind([m+1, n], i+1, j-1);
%                 col3 = sub2ind([m+1, n], i,   j);
%                 col4 = sub2ind([m+1, n], i+1, j);
%
%                 I = [I; row; row; row; row];
%                 J = [J; col1; col2; col3; col4];
%                 S = [S; 0.25; 0.25; 0.25; 0.25];
%             end
%         end
%     end
%
%     Ix2y = sparse(I, J, S, rowsGy, rowsGx);
% end
%
% function Iy2x = buildYFacesToXFaces(m, n)
% % Map y-face field (m x (n+1)) -> x-face field ((m+1) x n)
% % Uniform-grid, node-averaging transfer.
% %
% % Vectorization convention:
% %   y-faces stored as reshape(Vy, [], 1) with size(Vy) = [m, n+1]
% %   x-faces stored as reshape(Vx, [], 1) with size(Vx) = [m+1, n]
%
%     rowsGx = (m+1)*n;
%     rowsGy = m*(n+1);
%
%     I = [];
%     J = [];
%     S = [];
%
%     for i = 1:(m+1)      % x-face x-index (nodal in x)
%         for j = 1:n      % x-face y-index (centered in y)
%
%             row = sub2ind([m+1, n], i, j);
%
%             if i == 1
%                 % left boundary: use first y-face row only
%                 col1 = sub2ind([m, n+1], 1, j);
%                 col2 = sub2ind([m, n+1], 1, j+1);
%
%                 I = [I; row; row];
%                 J = [J; col1; col2];
%                 S = [S; 0.5; 0.5];
%
%             elseif i == m+1
%                 % right boundary: use last y-face row only
%                 col1 = sub2ind([m, n+1], m, j);
%                 col2 = sub2ind([m, n+1], m, j+1);
%
%                 I = [I; row; row];
%                 J = [J; col1; col2];
%                 S = [S; 0.5; 0.5];
%
%             else
%                 % interior: average four surrounding y-faces
%                 col1 = sub2ind([m, n+1], i-1, j);
%                 col2 = sub2ind([m, n+1], i-1, j+1);
%                 col3 = sub2ind([m, n+1], i,   j);
%                 col4 = sub2ind([m, n+1], i,   j+1);
%
%                 I = [I; row; row; row; row];
%                 J = [J; col1; col2; col3; col4];
%                 S = [S; 0.25; 0.25; 0.25; 0.25];
%             end
%         end
%     end
%
%     Iy2x = sparse(I, J, S, rowsGx, rowsGy);
% end

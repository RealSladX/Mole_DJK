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
m = 2*(k) + 2;            % cells in x
n = 2*(k) + 1;            % cells in y
dx = 1/m;
dy = 1/n;

% cell centers (include boundary nodes)
xc = [0 dx/2: dx :1-dx/2 1]';
yc = [0 dy/2: dy :1-dy/2 1]';
[Yc,Xc] = meshgrid(yc, xc);

% faces
xf = (0:dx:1)';
yf = (0:dy:1)';
[Yf,Xf] = meshgrid(yf, xf);

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

alpha = 1; % 1, 10, or 100

% Gradient Operator
G = grad2D(k, m, dx, n, dy);
% Divergence Operator
D = div2D(k, m, dx, n, dy);

% Building tensor components at Faces
C11F = Yf.^2 + alpha*Xf.^2;
C1221F = (alpha - 1)*Xf.*Yf;
C22F = Xf.^2 + alpha*Yf.^2;

% Get dimension of Gradient Operator
rowsGx = find(G(:, 2), 1)-1;
rowsGy = size(G, 1)-rowsGx;

% Full Heterogeneous Anisotropic 2x2 Tensor at Faces
K11 = spdiags(reshape(C11F, [], 1), 0, rowsGx, rowsGx);
K12 = spdiags(reshape(C1221F, [], 1), 0, rowsGx, rowsGy);
K21 = spdiags(reshape(C1221F, [], 1), 0, rowsGy, rowsGx);
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


% PLOTTING
figure(69);
contour3(Xc, Yc, ue);
title("Exact Solution");
shading interp;
view([0 90]);
colorbar;


figure(420);
contour3(Xc, Yc, ua);
title("Approximate Solution");
shading interp;
view([0 90]);
colorbar;

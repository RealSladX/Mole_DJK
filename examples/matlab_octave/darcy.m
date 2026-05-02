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
%   Assemble: L = - D * K * G, where
%     G : mimetic gradient, mapping cell-centers u -> faces
%     D : mimetic divergence, mapping faces -> centers
%     K : block tensor conductivity at faces
% =====================================================================

function [errMax, errL2, qxErrL2, qyErrL2, uminVal, umaxVal] = darcy(m, n, alpha, plot=false)
  addpath('../../src/matlab_octave');

  % ---- grid ----
  k = 4;            % mimetic order (typical: 2, 4)
  % m = 32;            % cells in x
  % n = 32;            % cells in y
  dx = 1/m;
  dy = 1/n;

  % Gradient Operator
  G = grad2D(k, m, dx, n, dy);
  % Divergence Operator
  D = div2D(k, m, dx, n, dy);

  % alpha = 100; % 1, 10, or 100

  % cell centers (include boundary nodes)
  xc = [0 dx/2: dx :1-dx/2 1]';
  yc = [0 dy/2: dy :1-dy/2 1]';
  [Yc,Xc] = meshgrid(yc, xc);

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


  % interior cell centers 
  xc_int = [dx/2: dx :1-dx/2]';
  yc_int = [dy/2: dy :1-dy/2]';

  % faces
  xf = (0:dx:1)';
  yf = (0:dy:1)';
  [Yf,Xf] = meshgrid(yf, xf);
  [Yf_x, Xf_x] = meshgrid(yc_int, xf);    % (m+1) x n  → x-faces
  [Yf_y, Xf_y] = meshgrid(yf, xc_int);    % m x (n+1)  → y-faces

  % Building tensor components at Faces
  C11F = Yf_x.^2 + alpha*Xf_x.^2;
  C12F_x = (alpha - 1).*Xf_x.*Yf_x;
  C12F_y = (alpha - 1).*Xf_y.*Yf_y;
  % C12F_x = (alpha - 1).*Xf.*Yf;
  % C12F_y = (alpha - 1).*Xf.*Yf;
  C22F = Xf_y.^2 + alpha*Yf_y.^2;

  Nx = (m+1)*n;
  Ny = m*(n+1);
  Icf = interpolCentersToFacesD2D(k, m, n);
  Ifc = interpolFacesToCentersG2D(k, m, n);
  N  = (m+2)*(n+2);
  Cx2x = Icf(1:Nx, 1:N);                 % center scalar -> x-faces
  Cy2y = Icf(Nx+1:Nx+Ny, N+1:2*N);       % center scalar -> y-faces
  Fx2c = Ifc(1:N, 1:Nx);                 % x-faces -> center scalar
  Fy2c = Ifc(N+1:2*N, Nx+1:Nx+Ny);       % y-faces -> center scalar
  Iy2x = Cx2x * Fy2c;   % y-faces -> centers -> x-faces
  Ix2y = Cy2y * Fx2c;   % x-faces -> centers -> y-faces

  % Full Heterogeneous Anisotropic 2x2 Tensor at Faces
  K11 = spdiags(reshape(C11F, [], 1), 0, Nx, Nx);
  K12 = spdiags(reshape(C12F_x, [], 1), 0, Nx, Nx) * Iy2x;
  K21 = spdiags(reshape(C12F_y, [], 1), 0, Ny, Ny) * Ix2y;
  K22 = spdiags(reshape(C22F, [], 1), 0, Ny, Ny);

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

  err = ua - ue;
  errMax = max(abs(err(:)));
  errL2  = norm(err(:)) / norm(ue(:));

  % numerical gradients on faces
  g = G * ua(:);
  uxh = reshape(g(1:Nx), m+1, n); % x-face gradient
  uyh = reshape(g(Nx+1:end), m, n+1); % y-face gradient

  % transported gradients
  uy_on_x = reshape(Iy2x * uyh(:), m+1, n);
  ux_on_y = reshape(Ix2y * uxh(:), m, n+1);

  % numerical fluxes
  qxh = -(C11F .* uxh + C12F_x .* uy_on_x);
  qyh = -(C12F_y .* ux_on_y + C22F .* uyh);
  qx_c = 0.5 * (qxh(1:m, :) + qxh(2:m+1, :));      % m x n
  qy_c = 0.5 * (qyh(:, 1:n) + qyh(:, 2:n+1));      % m x n

  [Xcc, Ycc] = meshgrid(xc_int, yc_int);
  Xcc = Xcc'; Ycc = Ycc';
  % exact solution on face grids
  ue_xf = exp(a * ((Xf_x - 0.5).^2 + (Yf_x - 0.5).^2));
  ue_yf = exp(a * ((Xf_y - 0.5).^2 + (Yf_y - 0.5).^2));

  ux_exact_xf = a * (2*Xf_x - 1) .* ue_xf;
  uy_exact_xf = a * (2*Yf_x - 1) .* ue_xf;
  ux_exact_yf = a * (2*Xf_y - 1) .* ue_yf;
  uy_exact_yf = a * (2*Yf_y - 1) .* ue_yf;

  % exact fluxes
  qx_exact = -(C11F .* ux_exact_xf + C12F_x .* uy_exact_xf);
  qy_exact = -(C12F_y .* ux_exact_yf + C22F .* uy_exact_yf);

  % flux errors
  qxErrL2 = sqrt(dx*dy) * norm(qxh(:) - qx_exact(:), 2);
  qyErrL2 = sqrt(dx*dy) * norm(qyh(:) - qy_exact(:), 2);

  % min/max
  uminVal = min(ua(:));
  umaxVal = max(ua(:));

  if plot
    % fprintf('max abs err = %e\n', max(abs(err(:))));
    % fprintf('rel L2 err   = %e\n', norm(err(:))/norm(ue(:)));
    % PLOTTING
    figure(1);
    contour3(Xc, Yc, ue);
    hold on;
    quiver(Xcc, Ycc, qx_c, qy_c, 1.0, 'b');
    title(sprintf("Exact Solution (alpha = %.2f)", alpha));
    shading interp;
    view([0 90]);
    colorbar;
    hold off;
    saveas(gcf,sprintf('/mnt/shared/COMP/670/Figures/Darcy_Exact_h%d_alpha%d.png', m, alpha))
    figure(2);
    contour3(Xc, Yc, ua);
    hold on;
    quiver(Xcc, Ycc, qx_c, qy_c, 1.0, 'b');
    title(sprintf("Approximate Solution (alpha = %.2f)", alpha));
    shading interp;
    view([0 90]);
    colorbar;
    hold off;
    saveas(gcf,sprintf('/mnt/shared/COMP/670/Figures/Darcy_Approx_h%d_alpha%d.png', m, alpha))


    figure(3);
    contour3(Xc, Yc, err);
    title(sprintf('Error (alpha = %.2f)', alpha));
    view([0 90]);
    colorbar;

    saveas(gcf,sprintf('/mnt/shared/COMP/670/Figures/Darcy_Error_h%d_alpha%d.png', m, alpha))

    figure(4);
    [C,h] = contour(Xc, Yc, RHS, 10, 'LineWidth', 0.8);
    clabel(C, h, 'FontSize', 8);
    title('(b) Source term f(x)');
  end
end



function heatDiffusion2D(alpha, Tmax, hotspotTemp)
% heatDiffusion2D - Simulates 2D heat diffusion using explicit method.
%
% Syntax:
%  heatDiffusion2D(alpha, Tmax, hotspotTemp)
%
% Inputs:
%   alpha       - Thermal diffusivity (e.g., 0.01)
%   Tmax        - Simulation time in seconds (e.g., 0.1)
%   hotspotTemp - Initial hotspot temperature (e.g., 100)
%
% Example:
%   heatDiffusion2D(0.01, 0.1, 100)

    % Grid and time setup
    Lx = 1; Ly = 1;
    Nx = 50; Ny = 50;
    dx = Lx / (Nx - 1); 
    dy = Ly / (Ny - 1);
    Nt = 1000;
    dt = Tmax / Nt;

    % Stability coefficients
    rx = alpha * dt / dx^2;
    ry = alpha * dt / dy^2;

    % Check for stability
    if rx + ry > 0.5
        warning('Unstable configuration: rx + ry = %.3f > 0.5', rx + ry);
    end

    % Grid coordinates
    x = linspace(0, Lx, Nx);
    y = linspace(0, Ly, Ny);
    [X, Y] = meshgrid(x, y);

    % Initial temperature matrix
    u = zeros(Ny, Nx);
    u(round(Ny/2), round(Nx/2)) = hotspotTemp;

    % Time-stepping loop
    for n = 1:Nt
        u_new = u;
        for i = 2:Ny-1
            for j = 2:Nx-1
                u_new(i,j) = u(i,j) + ...
                    rx * (u(i,j+1) - 2*u(i,j) + u(i,j-1)) + ...
                    ry * (u(i+1,j) - 2*u(i,j) + u(i-1,j));
            end
        end

        % Apply Dirichlet boundary conditions
        u_new(:, [1 end]) = 0;
        u_new([1 end], :) = 0;
        u = u_new;

        % Visualization
        if mod(n, 50) == 0
            surf(X, Y, u, 'EdgeColor', 'none');
            title(sprintf('Time = %.3f s', n*dt));
            xlabel('x');
            ylabel('y'); 
            zlabel('Temperature');
            axis([0 Lx 0 Ly 0 hotspotTemp]);
            view(2); 
            colorbar; 
            drawnow;
        end
    end
end
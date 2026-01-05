function heatDiffusion1D(alpha, Tmax, hotspotTemp)
% heatDiffusion1D - Simulates 1D heat diffusion using explicit method.
%
% Syntax:
%  heatDiffusion1D(alpha, Tmax, hotspotTemp)
%
% Inputs:
%   alpha       - Thermal diffusivity
%   Tmax        - Total simulation time (seconds)
%   hotspotTemp - Initial hotspot temperature in the center

    % Parameters
    L = 1;           % Length of the rod
    Nx = 50;         % Number of spatial points
    dx = L / (Nx-1); % Spatial step size
    Nt = 1000;       % Number of time steps
    dt = Tmax / Nt;

    r = alpha * dt / dx^2;

    if r > 0.5
        warning('Unstable configuration: r = %.3f > 0.5', r);
    end

    % Spatial grid
    x = linspace(0, L, Nx)';
    u = zeros(Nx, 1);
    u(round(Nx/2)) = hotspotTemp;

    % Time stepping
    for n = 1:Nt
        u_new = u;
        for i = 2:Nx-1
            u_new(i) = u(i) + r * (u(i+1) - 2*u(i) + u(i-1));
        end

        % Dirichlet boundary conditions
        u_new(1) = 0;
        u_new(end) = 0;

        u = u_new;

        % Visualization
        if mod(n, 50) == 0
            plot(x, u);
            title(sprintf('Time = %.3f s', n*dt));
            xlabel('x');
            ylabel('Temperature');
            axis([0 L 0 hotspotTemp]);
            drawnow;
        end
    end
end
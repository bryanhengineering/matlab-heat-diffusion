% Heat Diffusion Project

%% One Dimension

% Simulation Setup
L = 1; % L is the length of the rod
Nx = 50; % Nx is the number of grid points along the rod
dx = L/(Nx-1); % spatial step size
               % (distance between adjacent points on the rod)

% Time Setup
T = 0.5;  % Total time for simulation (set to 0.5 for test)
Nt = 1000; % Number of time steps to divide the simulation into
dt = T/Nt; % Time step size (how much time elapses per iteration)

% Physical Parameters
alpha = 0.01; % Thermal diffusivity of material (how fast heat spreads)
r = alpha*dt/dx^2; % stability ratio (must be less than or equal to 0.5 
                    % for simualtion to remain numerically stable)

% Grid and Initial Conditions
x = linspace(0,L,Nx)'; % Column vector of Nx equally spaced from 0 to L 
                       % (spatial coordinates)
u = zeros(Nx,1); % Column Vector storing temperature values at each point 
                 % on the rod (zeros as it starts at 0 everywhere)
u(round(Nx/2)) = 100; % initial heat peak (hot spot) at the middle of the 
                      % rod to 100
u(3)= 60; % Additional hot spot at point 3 of the rod 
          % (for testing arbitarily)

% Time Evolution Loop
for n = 1:Nt % Loops over each time step from 1 to Nt
    u_new = u; % Makes a copy of current temperature value to update 
               % them seperately

    % finite difference formula (Forward Euler)
    for i = 2:Nx-1 
        u_new(i) = u(i) + r*(u(i+1)-2*u(i)+u(i-1));  

        % to calculate the new temperature at each interior point based on 
        % its neighbors
    end

% Enforce boundary conditions:
    u_new(1) = 0; u_new(end) = 0;  
    % Dirichlet boundary conditions: both ends of the rod are kept at 0 
    % temperature (fixed)

    % Update and Plot
    u = u_new; % Updates the current temperature distribution for the next 
    % time step

    if mod(n,50)==0  % every 50 steps
        plot(x,u);   % Plots the current temperature along the rod.
        axis([0 L 0 100]);  % Keeps the y-axis fixed between 0 and 100.
        title(sprintf('Time = %.2f s', n*dt));
        % Shows the current time in seconds
        drawnow;
    end
end

%% 2 Dimension

% Simulation Setup
Lx = 1; Ly = 1; % length of the 2D domain (x and y directions)
Nx = 50; Ny = 50; % Number of grid points (x and y directions)[50x50 nodes]
dx = Lx / (Nx-1); % Spatial Step Size in the x direction
dy = Ly / (Ny-1); % Spatial Step Size in the y direction

% Time Setup
Tmax = 0.1;  % Total time to simulate (can be modified)
Nt = 1000;   % Number of time steps
dt = Tmax / Nt; % Time step duration 
% (how much simulated time passes per step)

% Physical Parameters
alpha = 0.01; % Thermal diffusivity of material 
              % (how fast heat spreads) (can be modified)
rx = alpha * dt / dx^2; 
ry = alpha * dt / dy^2;
% Precomputed constants for the finite difference formula 
% Must satisfy rx + ry ≤ 0.5 for numerical stability

% Grid
x = linspace(0, Lx, Nx); % Vectors representing spatial coordinates in x 
y = linspace(0, Ly, Ny); % Vectors representing spatial coordinates in y
[X, Y] = meshgrid(x, y); % uppercase X and Y for 2D coordinate grids 
                         % used for plotting and surface visualization

% Initial Condition
u = zeros(Ny, Nx);     % Note: rows = y, cols = x (creating a 2D matrix of 
                       % temperature values for the grid (Ny x Nx) )
u(round(Ny/2), round(Nx/2)) = 100; 
% Heat in center as starting point for diffusion (can be modified)

% Time stepping
for n = 1:Nt % loop that runs the simulation from time 0 to Tmax
    u_new = u; % u_new is a copy of the current state used to
               % update values without overwriting u

    %Explicit method in 2D
    for i = 2:Ny-1
        for j = 2:Nx-1
            u_new(i,j) = u(i,j) + ...
                rx * (u(i,j+1) - 2*u(i,j) + u(i,j-1)) + ...
                ry * (u(i+1,j) - 2*u(i,j) + u(i-1,j));
            % Each point is updated based on its four neighbors and itself
        end
    end

    % Boundary Conditions: Dirichlet (fixed at 0)
    u_new(:,[1 end]) = 0;     % left/right edges
    u_new([1 end],:) = 0;     % top/bottom edges

    % Sets all four edges to stay at 0 temperature 
    % (representing perfect heat sinks or fixed cold edges).

    u = u_new; % Updates the current temperature distribution for the next 
    % time step

    % Plot every 50 steps (shows a live heatmap of the temperature field)
    if mod(n,50)==0
        surf(X, Y, u, 'EdgeColor', 'none'); %surf creates a 3D surface plot
        title(sprintf('Time = %.3f s', n*dt));
        xlabel('x'); 
        ylabel('y'); 
        zlabel('Temperature');
        axis([0 Lx 0 Ly 0 100]);
        view(2); % Top-down view
        colorbar; % shows the temperature scale
        drawnow; % Tells MATLAB to update the figure 
    end
end


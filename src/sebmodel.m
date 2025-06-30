% Model ice melt below debris by solving simplified debris-surface energy
% balance
%
% Michael McCarthy 2025

function [T_s,T_d,melt] = sebmodel(S_in,S_out,T_a,SD,timestep,h,debk,...
    debrho,debc,c1,c2,constT_i,constL_f,constrho,nLayers)

% Get number of timesteps
nTimesteps = length(S_in);

% Get thickness of each debris sublayer
hLayer = h/nLayers;

% Preallocate space for debris surface and internal temperatures
T_s = nan(nTimesteps,1);
T_d = nan(nLayers+1,nTimesteps);

% Initial condition is linear temperature gradient between debris surface
% and ice surface
T_s(1) = T_a(1);
T_d(:,1) = linspace(T_a(1),constT_i,nLayers+1);

% Specify tolerances for Newton's method
rangeT_s = 0.5;
tolT_s = 0.01;
maxIter = 100;
maxStepT_s = 1;

% Loop through time steps solving energy balance for surface temperature. 
% Use air temperature or previous surface temperature as initial point for 
% solver. Aim is to find the surface temperature value for which the energy
% balance function f(Ts) = 0
for iTimestep = 2:nTimesteps

    % If snow depth is more than zero, surface temperature is zero
    if SD(iTimestep) > 0
        T_s(iTimestep) = constT_i;
        
    % Otherwise...
    else
        
        % Create a function handle for f(Ts)
        fun = @(Ts) energybalance(Ts,S_in(iTimestep),T_a(iTimestep),c1,...
            c2,S_out(iTimestep),debk,hLayer,T_d(2:end-1,iTimestep-1),...
            T_s(iTimestep-1),nLayers,constT_i,timestep,debrho,debc);
        
        % First guess of T_s for each timestep is T_a
        T_s0 = T_a(iTimestep);
        
        % Solve for surface temperature
        T_s(iTimestep) = newtonsmethod(fun,T_s0,tolT_s,maxIter,rangeT_s,...
            maxStepT_s);
    end
    
    % Calculate temperature profile, given surface temperature solution
    T_d(2:end-1,iTimestep) = tempprof(T_s(iTimestep),T_s(iTimestep-1),...
        T_d(2:end-1,iTimestep-1),nLayers,hLayer,constT_i,debk,timestep,...
        debrho,debc);
end

% Add surface temperature and ice temperature to temperature profiles
T_d(1,:) = T_s;
T_d(end,:) = constT_i;

% Calculate melt per timestep. Melt cannot be negative
G_i = debk*(T_d(end-1,:)-T_d(end))/hLayer;
melt = G_i*timestep/(constrho*constL_f);
melt(melt <= 0) = 0;

end

% Compute residual psi of debris-surface energy balance
%
% Michael McCarthy 2025

function psi = energybalance(Ts,S_in,T_a,c1,c2,S_out,debk,hLayer,...
    TdP,TsP,nLayers,constT_i,timestep,debrho,debc)

% Calculate temperature profile given air temperature or previous
% surface temperature and previous temperature profile
tempProf = tempprof(Ts,TsP,TdP,nLayers,hLayer,constT_i,debk,...
    timestep,debrho,debc);

% psi is the debris surface energy flux, which must equal zero for each
% time step
psi = S_in-S_out+c1*(T_a-Ts)+c2-debk*(Ts-tempProf(1))/hLayer;

end

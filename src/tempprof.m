% Calculate debris temperature profile according to Reid and Brock (2010)
% 
% Michael McCarthy 2025
%
% Notes
% - Implements Reid and Brock (2010) Equations A8 to A12
% - TsC is surface temperature of current time step T_s(t+1)
% - TsP is surface temperature of previous time step T_s(t)
% - Td is temperature profile of previous time step T_d(t,:)

function tempProf = tempprof(TsC,TsP,TdP,nLayers,hLayer,constT_i,debk,...
    timestep,debrho,debc)

    % Calculate parameter C
    C = debk*timestep/(2*debrho*debc*hLayer^2);
    
    % Calculate parameters a, b and c
    a = C;
    b = 2*C+1;
    c = C;
    
    % Calculate parameter d
    d(1) = C*TsC+C*TsP+(1-2*C)*TdP(1)+C*TdP(2);
    d(2:nLayers-2) = C*TdP(1:nLayers-3)+(1-2*C)*TdP(2:nLayers-2)+C*TdP...
        (3:nLayers-1);
    d(nLayers-1) = 2*C*constT_i+C*TdP(nLayers-2)+(1-2*C)*TdP(nLayers-1);
    
    % Calculate parameters A and S
    A = zeros(nLayers-1,1);
    S = zeros(nLayers-1,1);
    A(1) = b;
    S(1) = d(1);
    for iLayer = 2:nLayers-1
        A(iLayer) = b-a/A(iLayer-1)*c;
        S(iLayer) = d(iLayer)+a/A(iLayer-1)*S(iLayer-1);
    end
    
    % Calculate new temperature profile
    tempProf = zeros(nLayers-1,1);
    tempProf(nLayers-1) = S(nLayers-1)/A(nLayers-1);
    for iLayer = 1:nLayers-2
        tempProf(nLayers-1-iLayer) = 1/A(nLayers-1-iLayer)*...
            (S(nLayers-1-iLayer)+c*tempProf(nLayers-1-iLayer+1));
    end
end



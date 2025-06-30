% Run simplified energy-balance model for ice melt below debris
%
% Michael McCarthy 2025

% Get parent folder and add paths
foParent = fileparts(pwd);
addpath([foParent '/src'])
addpath([foParent '/inputs'])

% Specify meteorological and UDG data filenames
fnAwsData = 'inputs/arolla_aws.csv';
fnUdgData = 'inputs/arolla_udg.csv';

% Specify debris properties
debh = 0.06; % Debris thickness (m)
debk = 0.96; % Debris thermal conductivity (W m^-1 K^-1)
debc = 948; % Debris specific heat capacity (J kg^-1 K^-1)
debrho = 1496; % Debris density (kg m^-3)

% Specify number of debris sublayers
nLayers = max(10,round(debh/0.01));

% Specify physical constants
constL_f = 334000; % Latent heat of fusion of water (J kg^-1)
constT_i = 273.15; % Temperature of melting ice (K)
constrho_i = 915; % Ice density (kg m^-3)
constrho_w = 999.7; % Density of water (kg m^-3)

% Specify free parameters
c1 = 64.8; % (W m^-2 K^-1)
c2 = 88.1; % (W m^-2)

% Read meteorological and UDG data
awsData = readtable(fnAwsData);
udgData = readtable(fnUdgData);

% Get meteorological data dates/times and timestep
awsData.dateTime = datetime(awsData.year,awsData.month,awsData.day,...
    awsData.hour,awsData.minute,awsData.second);
timestep = seconds(duration(awsData.dateTime(2)-awsData.dateTime(1)));

% Run simplified energy-balance model
[T_s,T_d,melt] = sebmodel(awsData.S_in,awsData.S_out,awsData.T_a,...
    awsData.SD,timestep,debh,debk,debrho,debc,c1,c2,constT_i,constL_f,...
    constrho_i,nLayers);

% Get surface temperature observations for validation
epsilon = 0.95; % Debris emissivity ()
sigma = 5.67e-8; % Stefan-Boltzmann constant (W m^-2 K^-4)
T_sObs = (awsData.L_out/(sigma*epsilon)).^0.25;

% Make modelled melt rates daily to match observed melt rates
dailyMelt = melt(:);
dateTime = awsData.dateTime;
meltModTT = timetable(dateTime,dailyMelt);
meltModTT = retime(meltModTT,'daily','sum');
dailyMelt = meltModTT.dailyMelt;

% Calculate cumulative daily melt rate
cumMelt = cumsum(dailyMelt);

% Normalise melt rates by day number
nDays = length(cumMelt);
dayNum = 0:nDays-1;
dayNum = dayNum(:);
rCmMelt = cumMelt./dayNum;
rCmMeltObs = udgData.cmObs./dayNum;

% Calculate performance statistics
r2T_s = 1-nansum((T_sObs-T_s).^2)./(nansum((T_sObs-nanmean...
        (T_sObs)).^2));

% Calculate melt rate MAE
maeMelt = nanmean(abs(rCmMeltObs-rCmMelt));

% Plot surface temperature
figure(1)
plot(dateTime,T_sObs-273.15,'r-'); hold on
plot(dateTime,T_s-273.15,'b');
ylabel('Surface temperature (\circC)')
xlim([dateTime(1) dateTime(end)])
text(0.02,0.9,fnAwsData(8:10),'Units','Normalized')
text(0.98,0.1,['NSE = ' num2str(r2T_s,2)],'Units','Normalized',...
    'HorizontalAlignment','right')
formatfigure(15,5)

% Plot melt rate
figure(2)
plot(udgData.dateTime,udgData.meltObs,'r-'); hold on
plot(meltModTT.dateTime,dailyMelt,'b-');
ylabel('Melt (m i.e. hr^-^1)')
xlim([meltModTT.dateTime(1)-days(1) meltModTT.dateTime(end)+days(1)])
text(0.02,0.9,fnAwsData(8:10),'Units','Normalized')
formatfigure(15,5)

% Plot cumulative melt rate
figure(3)
plot(udgData.dateTime,udgData.cmObs,'r-'); hold on
plot(meltModTT.dateTime,cumMelt,'b-');
ylabel('Cumulative melt (m i.e.)')
xlim([meltModTT.dateTime(1)-days(1) meltModTT.dateTime(end)+days(1)])
text(0.02,0.9,fnAwsData(8:10),'Units','Normalized')
text(0.98,0.1,['MAE = ' num2str(maeMelt,2) ' m d^-^1'],'Units',...
    'Normalized','HorizontalAlignment','right')
formatfigure(15,5)

% Plot debris temperatures
figure(4)
imagesc('XData',0:length(dateTime),'YData',[0 -debh],...
    'CData',T_d-273.15)
cBar = colorbar;
ylabel('Depth (m)')
cBar.Label.String = 'Temperature (\circC)';
xlim([0 length(dateTime)]);
ylim([-debh 0])
text(0.02,0.9,fnAwsData(8:10),'Units','Normalized')
formatfigure(15,5)

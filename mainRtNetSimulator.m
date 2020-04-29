clear
close all
clc

rtFolder = "../qd-realization/src";

addpath("classes",...
    "functions",...
    fullfile("functions", "rayTracerUtils"),...
    rtFolder,...
    fullfile(rtFolder, "utils"))

%% setup
qdCampaignFolder = fullfile(rtFolder, 'custom_scenarios/Journal1Lroom');
scenario = 'refl2_qd0_relTh-40_floorMetal';

params.processRatios = true;
params.saveHref = false;

params.bsAnt = Antenna(8, 8, 0.5, 0.5, @(t,f) 1);
params.utAnt = Antenna(4, 4, 0.5, 0.5, @(t,f) 1);

params.bsIdxs = [2, 4];
params.utIdxs = [1, 3];

params.bsSectorDir = 0;
params.bsDowntilt = pi/2;
params.utSectorDir = 0;
params.utDowntilt = pi/2;

params.Ptx = 30;
params.fc = 60e9;
params.BW = 29 * 13.889e6;
params.F = 5;

params.dataDirection = "DL";
params.bfMode = "SVD";

params.txRefIdx = 2;
params.rxRefIdx = 1;
params.bsInterfIdxs = 4;
params.utInterfIdxs = 3;

%% simulation
params.rtFolder = rtFolder;
out = launchRtNetSimulation(fullfile(qdCampaignFolder, scenario), params);

%% Plot (sanity check)
figure
plot(out.SINR_db, 'DisplayName', 'SINR'); hold on
plot(out.SNR_db, 'DisplayName', 'SNR')
xlabel('t [step]')
ylabel('[dB]')
legend('show')
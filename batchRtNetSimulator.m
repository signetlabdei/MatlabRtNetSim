clear
close all
clc

rtFolder = "../qd-realization/src";

addpath("classes",...
    "functions",...
    rtFolder,...
    fullfile(rtFolder, "utils"),...
    fullfile("functions", "rayTracerUtils"),...
    "3GPP38900ChannelModel/scenarios",...
    "3GPP38900ChannelModel/utils")

%%
parpool(8)
campaigns = ["Journal1GroundOnly", "Journal1Indoor1", "Journal1Lroom", "Journal1ParkingLot"];
bfModes = ["SVD"];

antConfigs(1).bsAnt = Antenna(8, 8, 0.5, 0.5, @(t,f) 1);
antConfigs(1).utAnt = Antenna(4, 4, 0.5, 0.5, @(t,f) 1);

antConfigs(2).bsAnt = Antenna(2, 2, 0.5, 0.5, @(t,f) 1);
antConfigs(2).utAnt = Antenna(2, 2, 0.5, 0.5, @(t,f) 1);

%% setup
for campaign = campaigns
    for bfMode = bfModes
        for antConfig = antConfigs
            qdCampaignFolder = fullfile(rtFolder, 'custom_scenarios', campaign);
            
            params.processRatios = true;
            params.saveHref = false;
            
            params.bsAnt = antConfig.bsAnt;
            params.utAnt = antConfig.utAnt;
            
            switch(campaign)
                case {"Journal1Indoor1", "Journal1GroundOnly"}
                    params.bsIdxs = [2];
                    params.utIdxs = [1];
                    
                    params.txRefIdx = 2;
                    params.rxRefIdx = 1;
                    params.bsInterfIdxs = [];
                    params.utInterfIdxs = [];
                case {"Journal1Lroom", "Journal1ParkingLot"}
                    params.bsIdxs = [2,4];
                    params.utIdxs = [1,3];
                    
                    params.txRefIdx = 2;
                    params.rxRefIdx = 1;
                    params.bsInterfIdxs = [4];
                    params.utInterfIdxs = [3];
                otherwise
                    error()
            end
            
            params.bsSectorDir = 0;
            params.bsDowntilt = pi/2;
            params.utSectorDir = 0;
            params.utDowntilt = pi/2;
            
            params.Ptx = 30;
            params.fc = 60e9;
            params.BW = 29 * 13.889e6;
            params.F = 9;
            
            params.dataDirection = "DL";
            params.bfMode = bfMode;
            
            %% simulations
            matFileName = sprintf('matlab_stats_%d_%d_%s.mat',...
                params.bsAnt.getNumAnt(), params.utAnt.getNumAnt(), params.bfMode);
            scenarios = dir(qdCampaignFolder);
            
            i = 1;
            while i <= length(scenarios)
                if ~scenarios(i).isdir ||...
                        any(scenarios(i).name == [".", "..", "Input/"]) ||...
                        ~isValidScenarioPath(fullfile(scenarios(i).folder, scenarios(i).name))
                    
                    warning('''%s'' has been discarded as a scenario', scenarios(i).name)
                    scenarios(i) = [];
                    
                else
                    i = i+1;
                end
            end
            
            parfor i = 1:length(scenarios)
                fprintf('Processing %2d/%2d (%s)\n', i, length(scenarios), scenarios(i).name);
                
                scenario = scenarios(i).name;
                
                out = launchRtNetSimulation(fullfile(qdCampaignFolder, scenario), params);
                out.scenario = scenario;
                
                rtNetResults(i) = out;
                
            end
            save(fullfile(qdCampaignFolder, matFileName),...
                    'params', 'rtNetResults')
            
        end
    end
end
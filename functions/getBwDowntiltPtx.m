function [BW,downtilt,Ptx] = getBwDowntiltPtx(fc,scenario)
% Parameters from 38.901

% Bandwidth from 38.901, Table 7.8-1
switch(fc)
    case 6e9
        BW = 20e6;
    case {30e9,70e9}
        BW = 100e6;
    otherwise
        warning("fc=%.1f is not default carrier frequency",fc/1e9);
        BW = NaN;
end

% Downtilt, Ptx from 38.901, Table 7.8-1
switch(scenario)
    case "UMiSC"
        downtilt = deg2rad(12);
        switch(fc)
            case 6e9
                Ptx = 44;
            case {30e9,70e9}
                Ptx = 35;
            otherwise
                error("Frequency %.1f GHz not valid for '%s'",...
                    fc/1e9,scenario);
        end
    case "UMa"
        downtilt = deg2rad(12);
        switch(fc)
            case 6e9
                Ptx = 49;
            case {30e9,70e9}
                Ptx = 35;
            otherwise
                error("Frequency %.1f GHz not valid for '%s'",...
                    fc/1e9,scenario);
        end
    case {"InOo","InMo"}
        downtilt = deg2rad(20);
        Ptx = 24;
    otherwise
        error("Scenario '%s' not recognized",scenario);
end

end
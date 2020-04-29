classdef SimulationEnvironment < handle
    %SIMULATIONENVIRONMENT Class containing simulation-level variables,
    %such as pathloss, UT-BS attachment, as well as a list of all the BSs
    %and UTs in the simulation scenario. Function to obtained received
    %power/interference are also included.
    
    
    %% Properties
    properties
        bs % array of SimulationNode objects representing the BSs in the
        % simulation scenario. nBs := length(bs)
        ut % array of SimulationNode objects representing the UTs in the
        % simulation scenario. nUt := length(ut)
        
        bsPos % quick lookup for the position of BSs
        utPos % quick lookup for the position of UTs
        dist2dMat % (nUt x nBs) 2D distance matrix [m]
        
        o2iType % (nUt x 1) string vector containing the O2I type given
        % by getO2iType
        isLos % (nUt x nBs) boolean matrix
        sf % (nUt x nBs) matrix containing shadow fading. Positive values
        % correspond to an increased received power
        pl % (nUt x nBs) matrix containing pathloss (including O2I pathloss)
        
        utAttachedTo % column vector of length nUt. The value j at index i
        % means that ut(i) is attached to bs(j).
        
        bfMode % how to compute the beamforming vector (SVD or codebook). Currently supporting only SVD

    end
    
    
    %% Public methods
    methods
        % Constructor
        function E = SimulationEnvironment(bs, ut, attachMetric, bfMode)
            if ~exist('attachMetric','var')
                warning('Default attachment metric (PL-SF) will be used..')
                attachMetric = 'pl-sf';
            end
            if ~exist('bfMode','var')
                warning('Default beamforming (SVD) will be used..')
                bfMode = 'SVD';
            end                
            E.bs = bs(:)'; % row
            E.ut = ut(:); % column
            
            E.bfMode = bfMode;
            
            E.setLargeScaleParameters();
            E.attach(attachMetric);
            
        end
        
        
        % Methods
        function utIdx = getUtFromBs(E,bsIdx)
            % Finds the first UT attached to the desired BS
            utIdx = find(E.utAttachedTo == bsIdx, 1);
        end
        
        
        function [pRx,H,Htime,tau,H_params,bfGain] = getRefPrx(E,bsIdx,utRefIdx,dataDirection,largeScaleOnly)
            switch(dataDirection)
                case "DL"
                    % do nothing
                case "UL"
                    error("Uplink not supported")
                otherwise
                    error("Data direction '%s' not recognized",dataDirection)
            end
            
            % NaN initialization
            Htime = NaN;
            H = NaN;
            H_params = NaN;
            tau = NaN;
            utRefBf = NaN;
            bsBf = NaN;
            switch(largeScaleOnly)
                case "largeLoss"
                    % everything defaults to NaN
                    
                case {"largeLoss + bfGain", "largeLoss + bfGain no small fading"}
                    % Compute channel between BS and ref UT
                    [H,tau,~,Htime,H_params] = SimulationNode.getChannelMatrix(E.bs(bsIdx), E.ut(utRefIdx),...
                        E.bs(bsIdx).type,E.isLos(utRefIdx,bsIdx),E.o2iType(utRefIdx),largeScaleOnly);
                    % Compute BF
                    [bsBf,utRefBf] = SimulationNode.getBfVectors(...
                        E.bfMode,H,E.bs(bsIdx),E.ut(utRefIdx));
                    % Store UT ref BF vector towards BS for interference
                    % computation
                    E.ut(utRefIdx).storeBfVector(utRefBf);
                    
                otherwise
                    error("largeScaleOnly '%s' specification not recognized.", largeScaleOnly)
                    
            end
            
            [pRx,bfGain] = E.getPrx(bsIdx,utRefIdx,H,tau,utRefBf,bsBf,dataDirection,largeScaleOnly);
        end
        
        
        function [I_dbm, I_perSec_lin,bfGain] = getInterf(E,utRefIdx,dataDirection,largeScaleOnly)
            I_perSec_lin = nan(length(E.bs), E.bs(1).getNumSectors());
            bfGain = nan(length(E.bs), E.bs(1).getNumSectors());
            I_mw = 0;
            switch(dataDirection)
                case "DL"
                    for bsIdx = 1:length(E.bs)
                        utList = E.getInterferers(bsIdx,utRefIdx);
                        
                        for sectorIdx = 1:length(utList)
                            rxUtIdx = utList(sectorIdx);
                            
                            if rxUtIdx ~= -1 && ~isnan(rxUtIdx)
                                [interRx, bfGain(bsIdx,sectorIdx)] = E.getInterfPrx(bsIdx,...
                                    rxUtIdx,utRefIdx,dataDirection,largeScaleOnly);
                                I_mw = I_mw + interRx;
                                                                
                                I_perSec_lin(bsIdx,sectorIdx) = interRx;
                            end
                        end
                    end
                    
                case "UL"
                    error("Data direction '%s' not supported yet",dataDirection)
                    
                otherwise
                    error("Data direction '%s' not recognized",dataDirection)
            end
            
            I_dbm = 10*log10(I_mw);
        end
        
        
        % Visualize
        function plotEnvironment(E)
            scatter3(E.utPos(:,1),E.utPos(:,2),E.utPos(:,3),"filled"); hold on
            scatter3(E.bsPos(:,1),E.bsPos(:,2),E.bsPos(:,3),"filled"); hold off
            xlabel("x [m]")
            ylabel("y [m]")
            zlabel("z [m]")
            
            view(0,90)
            legend("UT","BS")
        end
    end
    
    
    %% Private methods
    methods (Access = private)
        function setLargeScaleParameters(E)
            E.getDist2dMat();
            
            % consider fc from a BS and assume they all use the same
            fc = E.bs(1).fc;
            scenarios = ["UMa","UMi","RMa","InOo","InMo"];
            
            % Determine o2iType
            E.o2iType = E.getO2iType();
            
            % Determine LoS probability
            pLos = NaN(length(E.ut),length(E.bs));
            for scenario = scenarios % Handle multiple propagation type per
                % scenario (e.g., Dense Urban from 38.913)
                scenarioMask = [E.bs.type] == scenario;
                pLos(:,scenarioMask) = getProbLos(scenario,...
                    E.bsPos(scenarioMask,:),E.utPos);
            end
            assert( ~any(isnan(pLos),"all"),...
                "Not all scenarios are being considered")
            
            % Determine isLos
            E.isLos = rand(size(pLos)) < pLos;
            E.isLos(E.o2iType ~= "None",:) = false; % force NLoS if O2I
            
            % outdoor PL computation
            outdoorPl = NaN(length(E.ut),length(E.bs));
            sfStd = NaN(length(E.ut),length(E.bs));
            for scenario = scenarios % Handle multiple propagation type per
                % scenario (e.g., Dense Urban from 38.913)
                scenarioMask = [E.bs.type] == scenario;
                
                [outdoorPl(:,scenarioMask),sfStd(:,scenarioMask)] =...
                    pathloss(fc,scenario,E.bsPos(scenarioMask,:),...
                    E.utPos,E.isLos(:,scenarioMask)); % TODO: Check TR
            end
            assert( ~any(isnan(outdoorPl),"all"),...
                "Not all scenarios are being considered")
            assert( ~any(isnan(sfStd),"all"),...
                "Not all scenarios are being considered")
            
            % Pathloss and Shadow Fading
            o2iPl = E.getO2iPathloss(fc);
            E.pl = outdoorPl + o2iPl;
            E.sf = sfStd.*randn(size(E.pl));
            
        end
        
        
        function attach(E,metric)
            % Associate each UT with a single BS
            % Pick the least pathloss value (max SNR, given equal BS tx power)
            % Positive shadow fading (sf) means higher received power (3gpp)
            switch(metric)
                case 'pl'
                    [~,E.utAttachedTo] = min(E.pl,[],2,"omitnan");
                case 'pl-sf'
                    [~,E.utAttachedTo] = min(E.pl-E.sf,[],2,"omitnan");
                otherwise
                    error('Attachment metric %s is not recognized!',metric)
            end
            
%             % Unattach UTs too far from their BS
%             idx = sub2ind(size(E.dist2dMat),(1:length(E.ut))',E.utAttachedTo);
%             dMax = max(vecnorm(E.bs(1).pos - E.bsPos, 2, 2));
%             E.utAttachedTo(E.dist2dMat(idx) > dMax) = -1;
        end
        
        
        function getDist2dMat(E)
            E.bsPos = SimulationNode.getPosMatrix(E.bs);
            E.utPos = SimulationNode.getPosMatrix(E.ut);
            
            E.dist2dMat = getDistanceMatrix(E.utPos(:,1:2), E.bsPos(:,1:2));
        end
        
        
        function utList = getInterferers(E,bsIdx,utRef)
            % Trivial scheduler: simply chooses one UT at random per sector
            % NaN means that no UT is in that sector ID
            % -1 means that the sector is already communicating with the
            % reference UT
            
            utsAttachedToBs = find(E.utAttachedTo == bsIdx)';
            
            utList = nan(1,E.bs(bsIdx).getNumSectors());
            for utIdx = utsAttachedToBs
                sector = E.bs(bsIdx).getSector(E.ut(utIdx));
                
                if isnan(utList(sector))
                    utList(sector) = utIdx;
                end
                
                if nnz(isnan(utList)) == length(utList)
                    % all sectors have been covered
                    break
                end
            end
            
            % if bsIdx is the reference one
            if any(utsAttachedToBs == utRef)
                % do not transmit twice in the utRef sector
                refSector = E.bs(bsIdx).getSector(E.ut(utRef));
                utList(refSector) = -1; % set sector containing ref to -1
            end
        end
        
        function [pRx_mw, bfGain] = getInterfPrx(E,bsIdx,rxUtIdx,utRefIdx,dataDirection,largeScaleOnly)
            switch(dataDirection)
                case "DL"
                    % do nothing
                case "UL"
                    error("Uplink not supported")
                otherwise
                    error("Data direction '%s' not recognized",dataDirection)
            end
            
            HToRef = NaN;
            auToRef = NaN;
            utRefBf = NaN;
            bsBf = NaN;
            bfGain = NaN;
            switch(largeScaleOnly)
                case {"largeLoss + bfGain","largeLoss + bfGain no small fading"}
                    % Compute channel between bs and rxUt
                    [HToRx,~,posInfo] = SimulationNode.getChannelMatrix(E.bs(bsIdx), E.ut(rxUtIdx),...
                        E.bs(bsIdx).type,E.isLos(rxUtIdx,bsIdx),E.o2iType(rxUtIdx),largeScaleOnly);
                    % Compute channel between bs and utRef with currently active bs
                    % sector
                    [HToRef,tauToRef] = SimulationNode.getChannelMatrix(E.bs(bsIdx), E.ut(utRefIdx),...
                        E.bs(bsIdx).type,E.isLos(utRefIdx,bsIdx),E.o2iType(utRefIdx),...
                        largeScaleOnly,posInfo.bsSector);
                    % Compute bs BF towards rxUt
                    bsBf = SimulationNode.getBfVectors(E.bfMode,HToRx,E.bs(bsIdx),E.ut(rxUtIdx));
                    % Obtain utRef BF towards btRef
                    utRefBf = E.ut(utRefIdx).getStoredBfVector();
                    
                case "largeLoss"
                    % do nothing
                otherwise
                    error("largeScaleOnly '%s' specification not recognized.", largeScaleOnly)
                    
            end
            
            [pRx_dbm,bfGain] = E.getPrx(bsIdx,utRefIdx,HToRef,tauToRef,utRefBf,bsBf,dataDirection,largeScaleOnly);
            pRx_mw = 10^(pRx_dbm/10);
        end
        
        
        function o2iType = getO2iType(E)
            o2iType = strings(length(E.ut),1);
            
            utType = [E.ut.type]';
            bsType = [E.bs.type];
            
            o2iType(utType == "Outdoor") = "None";
            o2iType(utType == "Car") = "Car";
            
            % Split indoor in high and low O2I loss
            indoorBs = any(bsType == ["InOo","InMo"]',1);
            if ~all(indoorBs)
                indoorIdxMask = find(utType == "Indoor");
                highLossMask = rand(length(indoorIdxMask),1) < 0.5;
                o2iType(indoorIdxMask(highLossMask)) = "IndoorHighLoss";
                o2iType(indoorIdxMask(~highLossMask)) = "IndoorLowLoss";
            end
            
            assert( ~any(o2iType == "","all"),...
                "Not all scenarios are being considered")
        end
        
        
        function o2iPl = getO2iPathloss(E,freq)
            o2iPl = zeros(size(E.o2iType));
            outdoorBsTypes = ["UMi","UMa","RMa"];
            
            % NOTE: The O2I penetration is UT-specifically generated
            for scenario = outdoorBsTypes
                scenarioMask = [E.bs.type] == scenario;
                % Car
                o2iTypeMask = E.o2iType == "Car";
                utO2iPl = additionalO2IPathloss(nnz(o2iTypeMask),freq,scenario,"Car");
                o2iPl(o2iTypeMask, scenarioMask) = repmat(utO2iPl,1,nnz(scenarioMask));
                
                % Indoor High loss
                o2iTypeMask = E.o2iType == "IndoorHighLoss";
                utO2iPl = additionalO2IPathloss(nnz(o2iTypeMask),freq,scenario,"High");
                o2iPl(o2iTypeMask, scenarioMask) = repmat(utO2iPl,1,nnz(scenarioMask));
                
                % Indoor Low loss
                o2iTypeMask = E.o2iType == "IndoorLowLoss";
                utO2iPl = additionalO2IPathloss(nnz(o2iTypeMask),freq,scenario,"Low");
                o2iPl(o2iTypeMask, scenarioMask) = repmat(utO2iPl,1,nnz(scenarioMask));
                
            end
            
            
            function o2iPl = applyO2iLoss(o2iPl,utO2iPl,mask,utMask)
                utIdx = find(utMask);
                for i = 1:length(utIdx)
                    idx = utIdx(i);
                    o2iPl(idx,mask(idx,:)) = utO2iPl(i);
                end
            end
        end
        
        
        function [Prx_dbm,bfGain] = getPrx(E,bsIdx,utIdx,H,tau,utBf,bsBf,dataDirection,largeScaleOnly)
            % Ptx + bf gain
            bfGain = nan;
            switch(dataDirection)
                case "DL"
                    Ptx = E.bs(bsIdx).Ptx;
                    bfRx = utBf;
                    bfTx = bsBf;
                case "UL"
                    Ptx = E.ut(utIdx).Ptx;
                    bfRx = bsBf;
                    bfTx = utBf;
                otherwise
                    error("Data direction '%s' not recognized",dataDirection)
            end
            % large scale
            % TR 38.901 (v15.0.0) Table 7.5-6 Part-1 NOTE 2: The sign of
            % the shadow fading is defined so that positive SF means more
            % received power at UT than predicted by the path loss model
            largeLoss = E.pl(utIdx,bsIdx) - E.sf(utIdx,bsIdx);
            switch(largeScaleOnly)
                case "largeLoss"
                    Prx_dbm = Ptx - largeLoss;
                    
                case {"largeLoss + bfGain", "largeLoss + bfGain no small fading"}
                    bfGain = getBfGain(H,tau,bfRx,bfTx);
                    Prx_dbm = Ptx - largeLoss + bfGain;
                    
                otherwise
                    error("largeScaleOnly '%s' specification not recognized.", largeScaleOnly)
            end
        end
        
    end
    
end
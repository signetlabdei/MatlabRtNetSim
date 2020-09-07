classdef SimulationNode < handle
    %SIMULATIONNODE Generic node in a network simulation. It can be a UT or
    %a BS (gNB). It is equipped with an antenna array (Antenna object) that
    %is used in all sectors. Other positional and TX/RX information are
    %also contained in its properties
    
    
    %% Properties
    properties
        nodeIdx % Node index. Useful for some scenarios

        ant % Antenna object used by the node. The same antenna will be
        % used for all sectors, if more than one are required
        sectorDir % broadside directions of the sectors [rad]. The number
        % of sectors is identified as the number of directions
        downtilt % mechanical downtilt in [rad] with respect to the zenith
        % (upward vertical direction)
        
        pos % [x,y,z] position of the node in [m]
        vel % [vx,vy,vz] velocity of the node in [m]
        type % BS: {"UMa","UMi","RMa","Indoor"}, UT: {"Indoor","Outdoor","Car"}
        % for propagation model and O2I
        
        F % RX noise figure [dB]
        Ptx % TX power [dBm]
        fc % carrier frequency [Hz]
        BW % bandwidth [Hz]
        
        storedBfVec % See: storeBfVector, getStoredBfVector
    end
    
    
    %% Public methods
    methods
        function N = SimulationNode(ant,sectorDir,downtilt,pos,vel,type,F,Ptx,fc,BW)
            N.ant = ant;
            N.downtilt = downtilt;
            N.pos = pos;
            N.vel = vel;
            N.type = type;
            N.F = F;
            N.Ptx = Ptx;
            N.fc = fc;
            N.BW = BW;
            
            N.sectorDir = fastWrapToPi(sectorDir);
            
        end
        
        
        % Methods
        function n = getNumSectors(N)
            %GETNUMSECTORS Returns number of sectors for the given node
            n = length(N.sectorDir);
        end
        
        
        function sector = getSector(N,node)
            %GETSECTOR Returns the sector of "N" is contained "node"
            
            sectorVecs = [cos(N.sectorDir); sin(N.sectorDir)];
            diffPos = node.pos(1:2) - N.pos(1:2);
            % max cosine similarity
            % NOTE: all sectorVecs have norm=1, thus invariant to diffPos
            % norm
            [~,sector] = max(diffPos * sectorVecs);
            
        end
        
        
        function n = getNoise(N)
            %GETNOISE Returns noise power based on BW and F
            n = -174 + N.F + 10*log10(N.BW);
        end
        
        
        function storeBfVector(N,bfVec)
            %STOREBFVECTOR Stores a predefined BF vector
            %It should be used to store the utRef's BF towards its BS
            %SEE ALSO: GETSTOREDBFVECTOR
            N.storedBfVec = bfVec;
        end
        
        
        function bfVec = getStoredBfVector(N)
            %GETSTOREDBFVECTOR Returns the previously stored BF vector
            %SEE ALSO: STOREBFVECTOR
            bfVec = N.storedBfVec;
        end
    end
    
    
    %% Static methods
    methods(Static)
        function [H,tau,posInfo,Htime,H_params] = getChannelMatrix(bs,ut,...
                scenario,isLos,o2iType,largeScaleOnly,bsActiveSector)
            %GETCHANNELMATRIX Method to obtain H matrix between a given BS
            %and UT
            %Note: by deafult, downlink communication is assumed
            
            % Input check
            if ~exist("bsActiveSector","var")
                bsActiveSector = NaN;
            end
            posInfo = SimulationNode.getPosInfo(bs,ut,bsActiveSector);
            
            switch(largeScaleOnly)
                case "largeLoss + bfGain no small fading"
                    t0 = tic;
                    bsAng = posInfo.bsLosAngles;
                    utAng = posInfo.utLosAngles;
                    
                    bsSteerVec = bs.ant.getSteeringVector(bsAng.el,bsAng.az);
                    bsElPatt = bs.ant.ElementPattern_lin(bsAng.el,bsAng.az);
                    bsVec = bsElPatt*bsSteerVec;
                    
                    utSteerVec = ut.ant.getSteeringVector(utAng.el,utAng.az);
                    utElPatt = ut.ant.ElementPattern_lin(utAng.el,utAng.az);
                    utVec = utElPatt*utSteerVec;
                    
                    H = utVec*bsVec.';
                    H_params = nan;
                    tau = nan;
                    Htime = toc(t0);
                case "largeLoss + bfGain"
                    t0 = tic;
                    [H,tau,H_params] = compute3gppHmatrix(struct(...
                        "scenario", scenario,...
                        "distance2d", posInfo.dist2d,...
                        "hBs", bs.pos(3),...
                        "hUt", ut.pos(3),...
                        "bsLosAngles", posInfo.bsLosAngles,...
                        "utLosAngles", posInfo.utLosAngles,...
                        "bsNAnt", struct("nv", bs.ant.M, "nh", bs.ant.N),...
                        "utNAnt", struct("nv", ut.ant.M, "nh", ut.ant.N),...
                        "bsAntSpacing", struct("dv", bs.ant.dv, "dh", bs.ant.dh),...
                        "utAntSpacing", struct("dv", ut.ant.dv, "dh", ut.ant.dh),...
                        "bsElemPattern", bs.ant.ElementPattern_lin,...
                        "utElemPattern", ut.ant.ElementPattern_lin,...
                        "utVelocity", ut.vel,...
                        "isO2i", o2iType ~= "None",...
                        "isLos", isLos,...
                        "freq", bs.fc,...
                        "tr", "38.901"));
                    Htime = toc(t0);
            end
            
        end
        
        
        function angle = getAbsoluteAngle(node1,node2)
            %GETABSOLUTEANGLE Returns [elevation,azimuth] from the given
            %SimulationNode N to node
            directionVec = node2.pos - node1.pos;
            linearDistance = norm(directionVec);
            
            azimuth = atan2(directionVec(2),directionVec(1));
            elevation = acos(directionVec(:,3)./linearDistance);
            
            angle = struct("el",elevation,"az",azimuth);
        end
        
        
        function posInfo = getPosInfo(bs,ut,bsActiveSector)
            %GETPOSINFO Returns positional information between the given BS
            %and UT, i.e., 2D and 3D distance, BS LoS angles [el,az]
            %towards UT, UT LoS angles [el,az] towards BS
            
            % Input check
            if ~exist("bsActiveSector","var")
                bsActiveSector = NaN;
            end
            
            % distance
            posInfo.dist2d = norm(bs.pos([1,2]) - ut.pos([1,2]));
            posInfo.dist3d = norm(bs.pos - ut.pos);
            
            % BS sector angles
            posInfo.bsAbsAngle = SimulationNode.getAbsoluteAngle(bs,ut);
            if isnan(bsActiveSector)
                posInfo.bsSector = bs.getSector(ut);
            else
                posInfo.bsSector = bsActiveSector;
            end
            posInfo.bsLosAngles = struct("el", posInfo.bsAbsAngle.el-bs.downtilt,...
                "az", posInfo.bsAbsAngle.az-bs.sectorDir(posInfo.bsSector)); % TODO: check angle wrapping
            
            % UTsector angles
            posInfo.utAbsAngle = SimulationNode.getAbsoluteAngle(ut,bs);
            posInfo.utSector = ut.getSector(bs);
            posInfo.utLosAngles = struct("el", posInfo.utAbsAngle.el-ut.downtilt,...
                "az", posInfo.utAbsAngle.az-ut.sectorDir(posInfo.utSector)); % TODO: check angle wrapping
        end
        
        
        function [txBf,rxBf] = getBfVectors(bfMode,H,tx,rx)
            %GETBFVECTOR Returns BF vectors for both BS and UT given the
            %channel matrix H. Implementation may not consider the channel.
            
            % SVD-based BF. Optimal for analog BF
            % Considering narrowband channel
            switch(lower(bfMode)) % case insensitive
                case 'svd'
                    [U,~,V] = svd(sum(H,3), "econ");
                    rxBf = conj(U(:,1));
                    txBf = V(:,1);
                case 'geometric'
                    angleTxRx = SimulationNode.getAbsoluteAngle(tx,rx);
                    txBf = tx.ant.getSteeringVector(angleTxRx.el, angleTxRx.az, tx.downtilt, tx.sectorDir);
                    txBf = txBf / norm(txBf);

                    angleRxTx = SimulationNode.getAbsoluteAngle(rx,tx);
                    rxBf = rx.ant.getSteeringVector(angleRxTx.el, angleRxTx.az, rx.downtilt, rx.sectorDir);
                    rxBf = rxBf / norm(rxBf);
                case 'codebook'
                    txCodebook = tx.ant.codebook;
                    rxCodebook = rx.ant.codebook;
                    
                    % Compute BF gain for each tx/rx codeword pair
                    Hnarrow = sum(H,3);
                    bfGainMatrix = 20*log10(abs(rxCodebook.' * Hnarrow * txCodebook));
                    
                    % Find best pair
                    [~, idx] = max(bfGainMatrix, [], 'all', 'linear');
                    [rxIdx, txIdx] = ind2sub(size(bfGainMatrix), idx);
                    
                    txBf = txCodebook(:, txIdx);
                    rxBf = rxCodebook(:, rxIdx);
                otherwise
                    error("Mode '%s' not yet supported",bfMode)
            end
        end
        
        
        function pos = getPosMatrix(nodes)
            %GETPOSMATRIX Utility function to obtain the position matrix
            %given an array of SimulationNodes
            pos = reshape([nodes.pos],3,length(nodes))';
        end
    end
end
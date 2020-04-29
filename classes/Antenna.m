classdef Antenna < handle
    %ANTENNA class. Contains all the parameters of an antenna array as
    %defined in 3GPP TR 38.901, using the same nomenclature.
    %Antenna-related methods are also included here (e.g., steering
    %vectors), since they depend on the antenna parameters.
    
    
    %% Properties
    properties
        % Array level properties
        M % number of vertical elements
        N % number of horizontal elements
        dv % spacing between vertical elements normalized in lambda units
        dh % spacing between horizontal elements normalized in lambda units
        
        ElementPattern_lin % Single antenna element pattern. Function accepting
        % two inputs (theta,phi) and return a gain in [linear] units
        P % number of polarizations
        
        % Panel level properties
        Mg % number of vertical panels
        Ng % number of horizontal panels
        dgv % spacing between vertical panels normalized in lambda units
        dgh % spacing between horizontal panels normalized in lambda units
    end
    
    
    %% Public methods
    methods
        % Constructor
        % TODO: multiple panels not supported. Might not be needed.
        % TODO: add quantized phases?
        function A = Antenna(M,N,dv,dh,ElementPattern_lin,P,Mg,Ng,dgv,dgh)
            % Array configuration
            if exist("M","var")
                A.M = M;
            else
                A.M = 1;
            end
            
            if exist("N","var")
                A.N = N;
            else
                A.N = 1;
            end
            
            if exist("dv","var")
                A.dv = dv;
            else
                A.dv = 1;
            end
            
            if exist("dh","var")
                A.dh = dh;
            else
                A.dh = 1;
            end
            
            % Pattern
            if exist("ElementPattern_lin","var")
                A.ElementPattern_lin = ElementPattern_lin;
            else
                A.ElementPattern_lin = @(t,f) 1; % isotropic
            end
            
            % Polarization
            if exist("P","var")
                assert(P==1, "No dual polarization is supported");
                A.P = P;
            else
                A.P = 1;
            end
            
            % Panels
            if exist("Mg","var")
                assert(Mg==1, "Multi-panels antennas are not supported")
                A.Mg = Mg;
            else
                A.Mg = 1;
            end
            if exist("Ng","var")
                assert(Ng==1, "Multi-panels antennas are not supported")
                A.Ng = Ng;
            else
                A.Ng = 1;
            end
            if exist("dgv","var")
                A.dgv = dgv;
            else
                A.dgv = 1;
            end
            if exist("dgh","var")
                A.dgh = dgh;
            else
                A.dgh = 1;
            end
        end
        
        function pos = getElementPositions(A)
            vIdx = repelem(0:A.M-1, A.N).'; % [0,0,0,...,M-1,M-1,M-1]
            hIdx = repmat(0:A.N-1, 1, A.M).'; % [0,1,2,...,N-1,0,1,2,...,N-1]
            
            pos = [zeros(size(vIdx)), hIdx*A.dh, vIdx*A.dv];
        end
        
        
        % Methods
        function a = getSteeringVector(A,theta,phi,downtilt,sectorDir)
            if ~isrow(theta) || ~isrow(phi)
                error("theta and phi should be scalar or row vectors")
            end
            if length(theta) ~= length(phi)
                error("theta and phi should have the same length")
            end
            
            [theta1, phi1] = transformAngles(theta, phi, downtilt, sectorDir);
            rho1 = getWaveVector(theta1(:), phi1(:), 1);
            pos = A.getElementPositions();
            
            a = exp(-1j * pos * rho1);
        end
        
        
        function gain = getElementPattern(A, theta, phi, downtilt, sectorDir)
            [theta1, phi1] = transformAngles(theta, phi, downtilt, sectorDir);
            gain = A.ElementPattern_lin(theta1, phi1);
        end
            
        
        function N = getNumAnt(A)
            if A.Mg == 1 && A.Ng == 1
                N = A.M * A.N;
            else
                error("Mg or Ng ~= 1 not supported")
            end
        end
        
        
        % Visualize
        function plotAntenna(A)
            pos = A.getElementPositions();
            scatter3(pos(:,1), pos(:,2), pos(:,3) ,"filled")
            
            axis equal
            xlabel("horiz [$\lambda$]")
            ylabel("vert [$\lambda$]")
        end
    end
end
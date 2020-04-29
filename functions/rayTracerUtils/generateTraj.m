function generateTraj(wayPointList,varargin)
%GENERATETRAJ Function that generates a trajectory for each node, as
%specified by wayPointsList and save them according to the /Input directory
%struct of the qd-realization.
%
% INPUTS:
% - wayPointList: cell array. The number of entries in wayPointList
% corresponds to the number of nodes in the scenario. Each entry contains the Nx3 list of
% waypoints of the corresponding node. Currently, a single moving node is
% supported, therefore all entries must contain 1x3 arrays (static nodes)
% except one (moving node).
% - velocities: double, matrix (Optional). Each row specifies the velocity
% of the corresponding node of the scenaro. If velocities is given, also
% samplingTime must be specified.
% - samplingTime: double, scalar (Optional). Trajectory sampling time. If
% samplingTime is given, also velocities must be specified.
% - nPoints: double, scalar (Optional). Number of points in the
% trajectory.
% - inputPath: string (Optional). Path where to save the generated files.
% If not specified, files will be saved in the /Input folder in the current
% directory.
%
% Either velocities and samplingTime are given, or numPoints.
%
% OUTPUTS:
% - None. All files are created in the Input folder in the specified path.
% If no path is given, the Input folder is created in the current
% directory.


numNodes = numel(wayPointList); % wayPointList is a numNodes x 1 cell array

% input parsing
p = inputParser;
validCellPos = @(x) isCellPos(x);
addRequired(p,'wayPointList',validCellPos);
addParameter(p,'velocities',@(x) size(x,1)==numNodes);
addParameter(p,'samplingTime',@(x) x>0);
addParameter(p,'nPoints',@(x) x>0);
addParameter(p,'inputPath',@(x) isstring(x))
parse(p,wayPointList,varargin{:});

% Custom validation functions
    function isPos(pos)
        % Test if pos has three columns (x,y,z)
        if size(pos,2)~=3
            error('You must specify coordinates as (x,y,z).')
        end
    end
    function isCellPos(cellpos)
        for x = cellpos
            isPos(x{1})
        end
    end

notGiven = p.UsingDefaults;
fields = fieldnames(p.Results);
args = fields(~contains(fields,notGiven));

if xor(any(contains(args,'velocities')),any(contains(args,'samplingTime'))) % input must contain both velocities and sampling time, or neither
    error('You have to specify both velocities and samplingTime');
else
    if any(contains(args,'velocities')) && any(contains(args,'samplingTime')) % if both velocities and sampling time are given
        if any(contains(args,'nPoints'))
            error('You can either specify (velocities,samplingTime) or nPoints, not both.')
        end
        velocities = p.Results.velocities;
        samplingTime = p.Results.samplingTime;
        isMoving = ~(abs(velocities)<eps(velocities)); % boolean mask selecting moving nodes (velocity~=0)
        % currently, a single moving node is supported
    elseif any(contains(args,'nPoints')) % if nPoints is given
        nPoints = p.Results.nPoints;
    end
end

% inputPath can
if any(contains(args,'inputPath'))
    inputPath = p.Results.inputPath;
else
    inputPath = fullfile(pwd,'Input');
    warning('Saving files in %s.',inputPath);
end
if ~exist(inputPath,'dir')
    mkdir(inputPath);
    warning('Folder %s created.',inputPath);
end
nodeLocPath = fullfile(inputPath,"nodes.dat");
nodeVelocitiesPath = fullfile(inputPath,"nodeVelocities.dat");


initPos = [];
for nodeId = 1:numNodes
    wayPoints = wayPointList{nodeId};
    dists = sqrt(sum((wayPoints(2:end,:)-wayPoints(1:end-1,:)).^2,2));
    distTot = sum(sqrt(sum((wayPoints(2:end,:)-wayPoints(1:end-1,:)).^2,2)));
    
    if (exist('velocities','var') && exist('samplingTime','var'))
        if isMoving(nodeId)
            vel = velocities(nodeId);
            timeTot = distTot/vel;
            samplingDist = vel*samplingTime;
        elseif ~isMoving(nodeId)
            wayPoints = wayPointList{nodeId};
            if size(wayPoints,1)~=1
                error("You specified velocity=0 but multiple waypoints for node %d",nodeId);
            end
            nodePos=repmat(wayPoints,totNumPoints,1);
            initPos = [initPos;nodePos(1,:)];
            filename = sprintf('NodePosition%d.dat',nodeId);
            csvwrite(fullfile(inputPath, filename), nodePos) % tx position
            continue
        end
    elseif (exist('nPoints','var'))
        samplingDist = distTot/nPoints;
    end
    
    allSamples = 0:samplingDist:distTot;
    totNumPoints = length(allSamples);
    
    prevDist = 0;
    currDist = 0;
    nPointsPerSegm = nan(size(wayPoints,1)-1,1);
    nodePos = wayPoints(1,:);
    if size(wayPoints,1)==1
        nodePos = repmat(wayPoints,[nPoints,1]);
        initPos = [initPos;nodePos(1,:)];
        filename = sprintf('NodePosition%d.dat',nodeId);
        csvwrite(fullfile(inputPath, filename), nodePos) % tx position
        continue
    else
        for j = 1:length(dists)
            currDist = currDist+dists(j);
            currNPoints = sum(prevDist<=allSamples & allSamples<=currDist);
            nPointsPerSegm(j) = currNPoints;
            prevDist = currDist;
            currTraj = linearTraj(wayPoints(j,:),wayPoints(j+1,:),currNPoints);
            nodePos=[nodePos;currTraj(2:end,:)];
        end
    end
    initPos = [initPos;nodePos(1,:)];
    filename = sprintf('NodePosition%d.dat', nodeId);
    csvwrite(fullfile(inputPath, filename), nodePos) % tx position
    
    scatter(nodePos(:,1), nodePos(:,2) ,10, flip(1:length(nodePos)))
    axis equal
    hold on
end
csvwrite(nodeLocPath,initPos) % initial positions
initVel = initPos; % TODO: check if this is ok
csvwrite(nodeVelocitiesPath,initVel) % initial velocities

end


%% utils
function positions=linearTraj(wayPoint1,wayPoint2,nPoints)
x = computePath(wayPoint1(1),wayPoint2(1),nPoints);
y = computePath(wayPoint1(2),wayPoint2(2),nPoints);
z = computePath(wayPoint1(3),wayPoint2(3),nPoints);
positions = [x,y,z];
end

function path = computePath(wayPoint1,wayPoint2,nPoints)
if wayPoint1==wayPoint2
    path = wayPoint1*ones(nPoints,1);
else
    path = linspace(min(wayPoint1,wayPoint2),max(wayPoint1,wayPoint2),nPoints).';
    if wayPoint1>wayPoint2
        path = flip(path);
    end
end
end
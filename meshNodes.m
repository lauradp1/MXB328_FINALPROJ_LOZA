function [nodes] = meshNodes(meshLims,r,meshLines,noNodeZones)
%MESHNODES Generates x and z vectors corresponding to nodes of mesh with
%each node geometrically distributed based on corresponding growth factor.
% Inputs:
%   meshLims: d by 2 matrix containing min and max values for each
%               dimension
%   r: d sized vector of growth factors
%   meshLines: d sized vector with how many nodes to generate for each
%               dimension (output will be +1 as it includes boundary nodes)
%   noNodeZones: d sized cell array with 2 by n matrices containing ranges
%                   in which no nodes should be placed
% Outputs:
%   nodes: d sized cell array containing node vector for each dimension d

% Obtain mesh information and initialise output
[dimension,~] = size(meshLims);
for d = 1:dimension
    nodes{d} = zeros(1,meshLines(d));
    delta{d} = (meshLims(2,d)-meshLims(1,d))/(sum(r(d).^(0:1:meshLines(d)-1)));
end

% Generate nodes for each dimension
for d = 1:dimension
    dVal = meshLims(1,d);
    for i = 1:meshLines(d) + 1
        if isempty(noNodeZones{d})
            nodes{d}(i) = dVal;
        else
            if ~ismember(1, noNodeZones{d}(1,:) <= dVal & dVal <= noNodeZones{d}(2,:)); nodes{d}(i) = dVal; end
        end
        dVal = dVal + r(d)^(i-1)*delta{d};
    end
end

end


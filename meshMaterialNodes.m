function [nodes,refinements,meshLims,materialsPlot] = meshMaterialNodes(materials,c_r)
%MESHMATERIALNODES Generates material distribution and populates mesh
%with base nodes that account for every corner of each material. Also adds
%refining nodes around these base nodes.
% Inputs:
%   materials: structure containing d by 2 by n arrays with each matrix
%               consisting of the d1 limits (1,:), d2 limits (2,:) etc. of
%               the d-D material for n d-D boxes of the material
%   c_r: vector of length d containing refinement constants of each
%           dimension
% Outputs:
%   nodes: d sized cell array containing the nodes for each dimension
%   refinements: d sized cell array containing 2 by n array with 
%                   corresponding refinement ranges for the relevant nodes
%   meshLims: 2 by d matrix containing limits of each dimension (min vals
%               on row 1 and max vals on row 2)
%   meshPlot: figure of material distribution

% Define variables
materialsPlot = figure('visible', 'off'); hold on;
matNames = fieldnames(materials);
colours = jet(length(matNames));
[dimension,~] = size(materials.(matNames{1}));

% Initialise outputs
for d = 1:dimension
    nodes{d} = [];
    refinements{d} = [];
end

% Plot material distribution and get base nodes
for k = 1:numel(matNames)
    % Obtain material blocks and no. of blocks
    matBlocks = materials.(matNames{k});
    [~,~,blocks] = size(matBlocks);
    
    for block = 1:blocks
        % Draw a 2D box to display the defined material distribution if 2D mesh
        if dimension == 2
            x = [matBlocks(1,1,block) matBlocks(1,:,block) matBlocks(1,end:-1:1,block)];
            y = [matBlocks(2,:,block) matBlocks(2,end:-1:1,block) matBlocks(2,1,block)];
            plot(x, y, 'Color', colours(k,:));
        end
        
        % Store unique values - these are the base nodes
        for d_1 = 1:dimension
            for d_2 = 1:dimension
                if ~ismember(matBlocks(d_1,d_2,block),nodes{d_1}); nodes{d_1}(end+1) = matBlocks(d_1,d_2,block); end
            end
        end
    end
end

% Get mesh limits
meshLims = zeros(2,dimension);
for d = 1:dimension
    meshLims(:,d) = [min(nodes{d});max(nodes{d})];
end

% Refine base nodes based on c_r and collate the 'no node' zones (any node
% that would lie within a refinement)
for d = 1:dimension
    dNodes = nodes{d};
    for dNode = dNodes
        if ~ismember(dNode - c_r(d), dNodes)
            if dNode - c_r(d) > meshLims(1,d)
                nodes{d}(end+1) = dNode - c_r(d);
                refinements{d}(1,end+1) = nodes{d}(end);
            else
                refinements{d}(1,end+1) = dNode;
            end
        end
        if ~ismember(dNode + c_r(d), dNodes)
            if dNode + c_r(d) < meshLims(2,d)
                nodes{d}(end+1) = dNode + c_r(d);
                refinements{d}(2,end) = nodes{d}(end);
            else
                refinements{d}(2,end) = dNode;
            end
        end
    end
end

end


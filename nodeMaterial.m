function material = nodeMaterial(materials,nodeCoordinates)
%NODEMATERIAL Determines what material the given node falls within. Will
%return false if node is on the boundary of two or more materials.
% Inputs:
%   materials: structure containing d by 2 by n arrays with each matrix
%               consisting of the d1 limits (1,:), d2 limits (2,:) etc. of
%               the d-D material for n d-D boxes of the material
%   nodeCoordinates: 1 by d sized vector containing coordinates of node
% Outputs:
%   material: the structure key of whichever material the node falls within

% Obtain material names/keys and dimension of mesh
matNames = fieldnames(materials);
dimension = length(nodeCoordinates);
material = false;

% Iterate through each material region until node falls within one
for k = 1:numel(matNames)
    % Obtain material blocks (and no. of blocks)
    matBlocks = materials.(matNames{k});
    [~,~,blocks] = size(matBlocks);
    
    % Iterate through each block of material
    for block = 1:blocks
        if sum(matBlocks(:,1,block) < nodeCoordinates' & nodeCoordinates' < matBlocks(:,2,block)) == dimension
            material = matNames{k};
            break;
        end
    end
    
    % Break if material has been found
    if material
        break;
    end
end

end


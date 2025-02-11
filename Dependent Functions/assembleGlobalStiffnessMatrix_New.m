function [S,T] = assembleGlobalStiffnessMatrix_New(N_Coord, elements, elementProps, E, G)

    % Number of elements
    numElements = size(elements, 1);
    
    % Initialize global stiffness matrix
    S = cell(numElements, 1);
    T=S;
    % Loop over each element to compute local and global stiffness matrices
    for i = 1:numElements
        % Element start and end nodes
        startNode = elements(i, 1);
        endNode = elements(i, 2);
        
        % Element properties: Area, Iy, Iz, J
        A = elementProps(i, 1);
        Iy = elementProps(i, 2);
        Iz = elementProps(i, 3);
        J = elementProps(i, 4);
        
        % Node coordinates
        startCoord = N_Coord(startNode, :);
        endCoord = N_Coord(endNode, :);
        
        % Element length and direction cosines
        L = norm(endCoord - startCoord);
        directionCosines = (endCoord - startCoord) / L;
        lx = directionCosines(1);
        ly = directionCosines(2);
        lz = directionCosines(3);
        
        
        % Compute the local stiffness matrix (12x12)
        kLocal = computeLocalStiffnessMatrix(E, G, A, Iy, Iz, J, L);
        % Compute the transformation matrix (12x12)
        T{i,1} = computeTransformationMatrix(lx, ly, lz);
        
        % Transform the local stiffness matrix to global coordinates
        S{i,1} = kLocal;
    end
end
% Function to compute the local stiffness matrix (12x12)
function kLocal = computeLocalStiffnessMatrix(E, G, A, Iy, Iz, J, L)
    % Stiffness coefficients
    EA_L = E * A / L;
    EIz_L3 = E * Iz / L^3;
    EIy_L3 = E * Iy / L^3;
    GJ_L = G * J / L;
    
    % Local stiffness matrix
    kLocal = zeros(12, 12);
    
    % Axial
    kLocal(1,1) = EA_L; kLocal(1,7) = -EA_L;
    kLocal(7,1) = -EA_L; kLocal(7,7) = EA_L;
    
    % Bending about z-axis
    kLocal(2,2) = 12*EIz_L3; kLocal(2,8) = -12*EIz_L3;
    kLocal(2,6) = 6*EIz_L3*L; kLocal(2,12) = 6*EIz_L3*L;
    kLocal(6,2) = 6*EIz_L3*L; kLocal(6,6) = 4*EIz_L3*L^2;
    kLocal(6,8) = -6*EIz_L3*L; kLocal(6,12) = 2*EIz_L3*L^2;
    kLocal(8,2) = -12*EIz_L3; kLocal(8,6) = -6*EIz_L3*L;
    kLocal(8,8) = 12*EIz_L3; kLocal(8,12) = -6*EIz_L3*L;
    kLocal(12,2) = 6*EIz_L3*L; kLocal(12,6) = 2*EIz_L3*L^2;
    kLocal(12,8) = -6*EIz_L3*L; kLocal(12,12) = 4*EIz_L3*L^2;
    
    % Bending about y-axis
    kLocal(3,3) = 12*EIy_L3; kLocal(3,9) = -12*EIy_L3;
    kLocal(3,5) = -6*EIy_L3*L; kLocal(3,11) = -6*EIy_L3*L;
    kLocal(5,3) = -6*EIy_L3*L; kLocal(5,5) = 4*EIy_L3*L^2;
    kLocal(5,9) = 6*EIy_L3*L; kLocal(5,11) = 2*EIy_L3*L^2;
    kLocal(9,3) = -12*EIy_L3; kLocal(9,5) = 6*EIy_L3*L;
    kLocal(9,9) = 12*EIy_L3; kLocal(9,11) = 6*EIy_L3*L;
    kLocal(11,3) = -6*EIy_L3*L; kLocal(11,5) = 2*EIy_L3*L^2;
    kLocal(11,9) = 6*EIy_L3*L; kLocal(11,11) = 4*EIy_L3*L^2;
    
    % Torsion
    kLocal(4,4) = GJ_L; kLocal(4,10) = -GJ_L;
    kLocal(10,4) = -GJ_L; kLocal(10,10) = GJ_L;
end

function T = computeTransformationMatrix(lx, ly, lz)
    % Direction cosine matrix for the local x-axis
    l = [lx, ly, lz];
    
    % Handling the case where the element is nearly aligned with the global z-axis
    if abs(lx) < abs(lz) && abs(ly) < abs(lz)
        % Choose a vector that is not parallel to l (use global x-axis)
        v = [1, 0, 0];
    else
        % Otherwise, use global z-axis
        v = [0, 0, 1];
    end
    
    % Compute local z-axis (perpendicular to l and v)
    z = cross(l, v);
    z = z / norm(z); % Normalize
    
    % Compute local y-axis (perpendicular to l and z)
    y = cross(z, l);
    y = y / norm(y); % Normalize
    
    % Assemble the rotation matrix
    R = [l; y; z]; % 3x3 rotation matrix, where:
    % - First row is the direction cosines of the local x-axis.
    % - Second row is the direction cosines of the local y-axis.
    % - Third row is the direction cosines of the local z-axis.
    
    % Initialize the full transformation matrix (12x12)
    T = zeros(12, 12);
    
    % Populate the transformation matrix for both nodes (start and end)
    for i = 1:4
        T(3*i-2:3*i, 3*i-2:3*i) = R;
    end
end

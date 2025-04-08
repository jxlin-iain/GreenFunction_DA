function hex_coords = generate_hexagonal_structure(n)
    % Generates hexagonal lattice centered around (0,0) up to layer n.
    
    % Initialize coordinates with the center point
    hex_coords = [0, 0];
    
    % Define the six primary direction vectors
    directions = [0.5, -sqrt(3)/2;
                   1, 0; 
                  0.5, sqrt(3)/2; 
                  -0.5, sqrt(3)/2; 
                  -1, 0; 
                  -0.5, -sqrt(3)/2; 
                  ];
    
    % Loop through each layer (starting from layer 2 to n)
    for layer = 1:n-1
        % Start at the "leftmost" point of the layer
        pos = [-layer, 0]; 
        
        % Loop through each of the six directions to form a hexagon
        for dir = 1:6
            for step = 1:layer
                hex_coords = [hex_coords; pos]; % Add current position
                pos = pos + directions(dir, :); % Move to next position
            end
        end
    end
    
end
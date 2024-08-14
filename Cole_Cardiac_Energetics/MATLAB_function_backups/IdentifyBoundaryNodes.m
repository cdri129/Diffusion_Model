% Read in .neigh file into Neighbors array
clear all;
close all;
meshdir = ['/Users/vrajagopal/Documents/heart/meshing/cellSpecificPLC/schneider_tomo/Cell_schneider_tomo/lower-res/']
fileroot = 'combined_tet_mesh_wryrgap.1';
filename = [meshdir fileroot '.neigh'];
delimiterIn = ' ';
headerlinesIn = 1;
NeighborsFile = importdata(filename, delimiterIn, headerlinesIn);

% Read in .ele file into Elements array
filename = [meshdir fileroot '.ele'];;
delimiterIn = ' ';
headerlinesIn = 1;
ElementsFile = importdata(filename, delimiterIn, headerlinesIn);

Neighbors = NeighborsFile.data;
 [rows, cols] = size(Neighbors);
 
 Elements = ElementsFile.data;
 count = 0;
 
 %Initializing result array 
 BoundaryNodes = zeros(1,1);
 
while rows > 0,
 while cols > 1,
     
     % Checking to see if there are -1 neighbors and how many are -1
     if Neighbors(rows,cols) == -1
         pos = cols;
         count = count + 1;
     else
      end
     cols = cols - 1;
     
 end
     cols = 5;
     
% If count > 1 then all nodes are on the boundary. Just return the entire
% row from Elements array. Else return all elements except the one at pos.

    
    if count > 0
        BNodes = Elements(rows,:);
        if count > 1
        BNodes
        else
        %resetting element at pos to 0
        BNodes(:, pos) = [0]
        end 
        for nodeidx = 2:(size(BNodes,2) - 1)
        if BNodes(1,nodeidx) ~= 0
        BoundaryNodes(end + 1,1) = BNodes(1,nodeidx); 
        end
    end
    
    
    end
    
 rows = rows - 1;
 count = 0;
end
 
Boundary = unique(BoundaryNodes( [2:end] , : ));
Boundary = [size(Boundary,1);Boundary]
dlmwrite([meshdir fileroot '.bdnode'], Boundary,'precision',10);
      


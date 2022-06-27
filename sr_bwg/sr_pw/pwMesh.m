% This function aims to generate an unstructured 3D Mesh in pointwise having
% already defined an empty volume mesh
%       Author  : Gabriel Buendia
%       Version : 1
% Inputs:
%      ID         -> ID of the already existing bulk data entity
%      BLoptions  -> Structure containing options for the boundary layer
%                    Fields:
%                    -maxLayers       -> This default is the maximum number of 
%                    T-Rex layers of an unstructured block when it is created : 
%                    Number (Default 0)
%                    -fullLayers       -> This default is the minimum number of 
%                    fully structured T-Rex layers of an unstructured block when 
%                    it is created : 
%                    Number (Default 0)
%                    -growthRate       -> This default is the growth rate of 
%                    T-Rex layers of an unstructured block when it is created : 
%                    Number (Default 1.2) 
%                    -push             -> This default is the flag for pushing 
%                    T-Rex attributes onto the connectors and domains of an 
%                    unstructured block when a it is created : 
%                    'true' or 'false' (Default false)
%                    -solverAttribute  -> Sets the named unstructured solver 
%                    attribute : 
%                    'TetPyramidPrismHex'||'AllAndReducePyramids'||
%                    'AllAndConvertWallDoms'||'LegacyTetPyramidPrismHex'||
%                    'TetPyramid'
%                    -firstLayerHeight -> Specifies the BL thickness of the
%                    first cell : 
%                    Number (User specified)
%      ISOoptions -> -minLength -> Minimum length of the isotropic tets :
%                    Number||'Boundary' (Default Boundary)
%                   -maxLength -> Maximum length of the isotropic tets :
%                    Number||'Boundary' (Default Boundary)
% Outputs:
%      volumeMeshCommand -> String containing the command to write the 3D
%                           mesh
function volumeMeshCommand = pwMesh(ID, BLoptions, ISOoptions)
volumeMeshCommand = append(['set _TMP(mode_1) [pw::Application begin UnstructuredSolver [list $_BL(',num2str(ID),')]]\n' ...
      '  $_BL(',num2str(ID),') setUnstructuredSolverAttribute TRexMaximumLayers ', num2str(BLoptions.maxLayers), '\n' ...
      '  $_BL(',num2str(ID),') setUnstructuredSolverAttribute TRexFullLayers ', num2str(BLoptions.fullLayers), '\n' ...
      '  $_BL(',num2str(ID),') setUnstructuredSolverAttribute TRexGrowthRate ',  num2str(BLoptions.growthRate), '\n' ...
      '  $_BL(',num2str(ID),') setUnstructuredSolverAttribute TRexPushAttributes ', num2str(BLoptions.push), '\n' ...
      '  $_BL(',num2str(ID),') setUnstructuredSolverAttribute TRexCellType ', num2str(BLoptions.solverAttribute), '\n' ...
      '  set _TMP(PW_1) [pw::TRexCondition getByName {Boundary Layer}]\n' ...
      '  $_TMP(PW_1) setValue ', num2str(BLoptions.firstLayerHeight), '\n' ...
      '  $_BL(',num2str(ID),') setUnstructuredSolverAttribute EdgeMinimumLength ', num2str(ISOoptions.minLength), '\n' ...
      '  $_BL(',num2str(ID),') setUnstructuredSolverAttribute EdgeMaximumLength ', num2str(ISOoptions.maxLength), '\n' ...
      '  $_BL(',num2str(ID),') setUnstructuredSolverAttribute EdgeMaximumGrowthRate ', num2str(ISOoptions.maxGrowth), '\n' ...
      '  $_TMP(mode_1) setStopWhenFullLayersNotMet true\n' ...
      '  $_TMP(mode_1) setAllowIncomplete true\n' ...
      '  $_TMP(mode_1) run Initialize\n' ...
      '$_TMP(mode_1) end\n' ...
      'unset _TMP(mode_1)\n']);
end
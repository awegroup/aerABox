% This function aims to generate an empty farfield (3D mesh) from a list of 
% domains for Pointwise language
%       Author  : Gabriel Buendia
%       Version : 1
% Inputs:
%       ID          -> Pointwise ID desired farfield number
%       list        -> List of domains IDs enclosing the geometry of interest
%                      You may also specify 'all' by typing 0
%       optionsFF   -> 'Box' : Box shape farfield
%                      'Cylinder' : Cylindrical shape farfield
%                      'Sphere' : Spherical shape farfield
%       optionsCell -> 'Unstructured' : Unstructured elements in the FF
%                      'Voxel' : Voxel elements in the FF
%       dimensions  -> They depend on the 'optionsFF' specified
%                      if 'Box'      -> 3x2 Matrix : [x_positive x_negative;
%                                                     y_positive y_negative;
%                                                     z_positive z_negative]
%                      if 'Cylinder' -> 3x2 Matrix : [L_positive L_negative;
%                                                     Radius_1   Radius_2;
%                                                     Axis_alignment    0]
%                      if 'Sphere'   -> 1x1 Variable : Radius
%       arrangement -> 1x3 Vector specifying the largest and the shortest
%                      directions e.g. [1 2 3] where 1 is the largest and 3
%                      the lowest (only useful for 'Box' case)
% Outputs:
%       farfieldCommand -> String containing the command to write the
%                          farfield
function farfieldCommand = pwFarfield(ID,list,optionsFF,optionsCell,dimensions,arrangement)
size = length(list);
string = '[list ';
if size ~= 1    
    for i = 1:size
        if i < size
            auxstr = append('$_DM(',num2str(list(i)),') '); 
        else
            auxstr = append('$_DM(',num2str(list(i)),')]]');
        end
        string = append(string,auxstr);
    end
    allBoundary = '';
else
    allBoundary = 'set _TMP(mode_1) [pw::Grid getAll -type pw::Domain]\n';
    string = append(string,'$_TMP(mode_1)]]');
end
switch optionsFF
    case 'Box'
        x_forw = dimensions(arrangement(1),1); x_back = dimensions(arrangement(1),2);
        y_forw = dimensions(arrangement(2),1); y_back = dimensions(arrangement(2),2);
        z_forw = dimensions(arrangement(3),1); z_back = dimensions(arrangement(3),2);
        farfieltype = append(['  $_TMP(mode_2) setFarfieldShapeType Box\n'...
            '  $_TMP(mode_2) setFarfieldLength {',num2str(x_forw),' ',num2str(x_back),'}\n'...
            '  $_TMP(mode_2) setFarfieldWidth {',num2str(y_forw),' ',num2str(y_back),'}\n'...
            '  $_TMP(mode_2) setFarfieldHeight {',num2str(z_forw),' ',num2str(z_back),'}\n']);
    case 'Cylinder'
        L_forw = dimensions(1,1); L_back = dimensions(1,2);
        r_forw = dimensions(2,1); r_back = dimensions(2,2);
        axis = dimensions(3,1);
        switch axis
            case 1
                alignment = 'X';
            case 2
                alignment = 'Y';
            case 3
                alignment = 'Z';
            otherwise
                alignment = 'Automatic';
        end
        farfieltype = append(['  $_TMP(mode_2) setFarfieldShapeType Cylinder\n'...
            '  $_TMP(mode_2) setFarfieldLength {',num2str(L_forw),' ',num2str(L_back),'}\n'...
            '  $_TMP(mode_2) setFarfieldCapRadii {',num2str(r_forw),' ',num2str(r_back),'}\n'...
            '  $_TMP(mode_2) setShapeAlignment ',alignment,'\n']);
    case 'Sphere'
        R = dimensions;
        farfieltype = append(['  $_TMP(mode_2) setFarfieldShapeType Sphere\n'...
            '  $_TMP(mode_2) setFarfieldRadius ',num2str(R),'\n']);
    otherwise
        x_forw = dimensions(1,1); x_back = dimensions(1,2);
        y_forw = dimensions(2,1); y_back = dimensions(2,2);
        z_forw = dimensions(3,1); z_back = dimensions(3,2);
        farfieltype = append(['  $_TMP(mode_2) setFarfieldShapeType Box\n'...
            '  $_TMP(mode_2) setFarfieldHeight {',num2str(x_forw),' ',snum2str(x_back),'}\n'...
            '  $_TMP(mode_2) setFarfieldLength {',num2str(y_forw),' ',snum2str(y_back),'}\n'...
            '  $_TMP(mode_2) setFarfieldWidth {',num2str(z_forw),' ',snum2str(z_back),'}\n']);
end
switch optionsCell
    case 'Unstructured'
        meshtype = append('  $_TMP(mode_2) setMeshType Unstructured\n');
    case 'Voxel'
        meshtype = append('  $_TMP(mode_2) setMeshType Voxel\n');
    otherwise
        meshtype = append('  $_TMP(mode_2) setMeshType Unstructured\n');
end

farfieldCommand = append([allBoundary, ...
    'set _TMP(mode_2) [pw::Application begin VolumeMesher ',string,'\n'...    
    farfieltype,...
    meshtype,...
    '  $_TMP(mode_2) createGridEntities\n'...
    '  $_TMP(mode_2) end\n'...
    'unset _TMP(mode_1)\n'...
    'unset _TMP(mode_2)\n'...
    'set _BL(',num2str(ID),') [pw::Grid getAll -type pw::Block]\n'...
    'pw::Application markUndoLevel {Automatic Volume Mesh}\n']);
end

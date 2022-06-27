% This function aims to generate an empty source (3D refinement region) from
% in any location and shape specified by the user for Pointwise language
%       Author  : Gabriel Buendia
%       Version : 1
% Inputs:
%       ID           -> Pointwise ID desired source number
%       optionsSrc   -> 'Box' : Box shape source
%                       'Cylinder' : Cylindrical shape source
%                       'Sphere' : Spherical shape source
%       startface    -> 'x' : Source start and end faces are defined
%                       following the unitary vector i (x axis)
%                       'y' : Source start and end faces are defined
%                       following the unitary vector j (y axis)
%                       'z' : Source start and end faces are defined
%                       following the unitary vector k (z axis)
%       center       -> Specify the COG of the source in cartesian coords
%                       1x3 Vector : [x_center y_center z_center]
%       dimensions   -> They depend on the 'optionsSrc' specified
%                       if 'Box'      -> 1x3 Vector : [length height width]
%                       if 'Cylinder' -> 1x3 Vector : [radius top_Radius length]
%                                        Note that radius and top_Radius are 
%                                        the same for cylinder and different
%                                        if defining a cone
%                       if 'Sphere'   -> 1x3 Vector : [radius base_Angle top_Angle]
%                                        Note that base_Angle -> [0,90] deg and
%                                        top_Angle -> [90,180] deg are used
%                                        to define cut sphere sources
%       distribution -> 'Constant': The begin values will be used 
%                        throughout the source
%                       'Parametric' : The begin values will be used at the 
%                        minimum parametric limits and the end values will be
%                        used at the maximum parametric limits of the source
%                       'AxisToPerimeter : The begin values will be used along 
%                       the axis and the end values will be used at the perimeter
%                       of the source
%                       'CenterToPerimeter : The begin values will be used at 
%                       the center and the end values will be used at the perimeter
%                       of the source 
%       spec         -> 2x2 Matrix: [spacing_initial spacing_final;
%                                    decay_initial   decay_final]
%                       Note that for the 'constant' case only the initial
%                       values are read
% Outputs:
%       sourceCommand -> String containing the command to write the
%                        source
function sourceCommand = pwSource(ID,optionsSrc,startface,center,dimensions,distribution,spec)
switch optionsSrc
    case 'Box'
        switch startface
            case 'x'
                R = rotationMatrix(0,270,0);
                T = translationMatrix(-center(3),center(2),center(1));
                len = dimensions(1,1); height = dimensions(1,2); width = dimensions(1,3);
            case 'y'
                R = rotationMatrix(90,0,0);
                T = translationMatrix(center(1),-center(3),center(2));
                len = dimensions(1,2); height = dimensions(1,3); width = dimensions(1,1);
            case 'z'
                R = rotationMatrix(0,0,0);
                T = translationMatrix(center(1),center(2),center(3));
                len = dimensions(1,3); height = dimensions(1,2); width = dimensions(1,1);
        end
        typeDim = append('  $_SR(',num2str(ID),') box  -length ',num2str(len),' -height ',num2str(height),' -width ',num2str(width));
    case 'Cylinder'
        radius = dimensions(1,1); topRadius = dimensions(1,2); len = dimensions(1,3);
        switch startface
            case 'x'
                R = rotationMatrix(0,270,0);
                T = translationMatrix(-center(3),center(2),center(1));
                
            case 'y'
                R = rotationMatrix(90,0,0);
                T = translationMatrix(center(1),-center(3),center(2));
            case 'z'
                R = rotationMatrix(0,0,0);
                T = translationMatrix(center(1),center(2),center(3));
        end
        typeDim = append('  $_SR(',num2str(ID),') cylinder -radius ',num2str(radius),' -topRadius ',num2str(topRadius),' -length ',num2str(len));
    case 'Sphere'
        radius = dimensions(1,1); baseAngle = dimensions(1,2); topAngle = dimensions(1,3);
        switch startface
            case 'x'
                R = rotationMatrix(0,270,0);
                T = translationMatrix(-center(3),center(2),center(1));
            case 'y'
                R = rotationMatrix(90,0,0);
                T = translationMatrix(center(1),-center(3),center(2));
            case 'z'
                R = rotationMatrix(0,0,0);
                T = translationMatrix(center(1),center(2),center(3));
        end
        typeDim = append('  $_SR(',num2str(ID),') sphere -radius ',num2str(radius),' -baseAngle ',num2str(baseAngle),' -topAngle ',num2str(topAngle));
end
Transform = R*T;
Transform = reshape(Transform,[1 16]);
size = length(Transform);
string = '[list ';
for i = 1:size
    if i < size
        auxstr = append(num2str(Transform(i)),' '); 
    else
        auxstr = append(num2str(Transform(i)),']');
    end
    string = append(string,auxstr);
end
sourceCommand = append(['set _TMP(mode_1) [pw::Application begin Create]\n'...
    '  set _SR(',num2str(ID),') [pw::SourceShape create]\n'...
    '  $_SR(',num2str(ID),') setPivot Center\n'...
    typeDim,'\n'...
    '  $_SR(',num2str(ID),') setTransform ',string,'\n' ...
    '$_TMP(mode_1) end\n' ...
    'unset _TMP(mode_1)\n']);
switch distribution 
    case 'Constant'
        spacingB = spec(1,1); spacingE = spacingB;
        decayB = spec(2,1); decayE = 0.5;
        typeDist = append('  $_SR(',num2str(ID),') setSpecificationType Constant');
    case 'Parametric'
        spacingB = spec(1,1); spacingE = spec(1,2);
        decayB = spec(2,1); decayE = spec(2,2);
        typeDist = append('  $_SR(',num2str(ID),') setSpecificationType Parametric');
    case 'AxisToPerimeter'
        spacingB = spec(1,1); spacingE = spec(1,2);
        decayB = spec(2,1); decayE = spec(2,2);
        typeDist = append('  $_SR(',num2str(ID),') setSpecificationType AxisToPerimeter');
    case 'CenterToPerimeter'
        spacingB = spec(1,1); spacingE = spec(1,2);
        decayB = spec(2,1); decayE = spec(2,2);
        typeDist = append('  $_SR(',num2str(ID),') setSpecificationType CenterToPerimeter');
end
sourceCommand = append([sourceCommand,...
    'set _TMP(mode_1) [pw::Application begin Modify [list $_SR(',num2str(ID),')]]\n'...
    typeDist,'\n'...
    '  $_SR(',num2str(ID),') setBeginSpacing ',num2str(spacingB),'\n'...
    '  $_SR(',num2str(ID),') setEndSpacing ',num2str(spacingE),'\n'...
    '  $_SR(',num2str(ID),') setBeginDecay ',num2str(decayB),'\n'...
    '  $_SR(',num2str(ID),') setEndDecay ',num2str(decayE),'\n'...
    '$_TMP(mode_1) end\n'...
    'unset _TMP(mode_1)\n']);
end

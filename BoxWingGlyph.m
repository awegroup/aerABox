% This function generates the glyph file to be run by Pointwise to generate
% the CFD mesh
%       Author  : Gabriel Buendia
%       Version : 1
% Inputs:
%       name         -> Name of the input file
%       Slices_Wing  -> Matrix [N_points_airfoil x 3 x N_slices] containing 
%                         all the wing cross sections [m]
%       b            -> Wingspan [m]
%       R            -> Cant radius [m]
%       h            -> Height between wings [m]
%       s            -> Wing stagger [m]
%       c_r1         -> Chord at the root of the lower wing [m]
%       c_t1         -> Chord at the tip of the lower wing [m]
%       c_r2         -> Chord at the root of the upper wing [m]
%       c_t2         -> Chord at the tip of the upper wing [m]
%       N            -> Number of sections to discretize the cant region
%       N_bu         -> Number of spanwise sections upper wing
%       N_bl         -> Number of spanwise sections lower wing
%       N_bw         -> Number of spanwise sections winglet
%       el_size      -> Size of the element [m]
%       tolBound     -> This attribute is the tolerance to use when splitting 
%                      and joining the curves specified as the boundaries, for
%                      detecting end to end connections of the boundaries, and 
%                      for creating the surface from the boundaries
%       tolThres     -> This attribute specifies the percentage of length of 
%                      the boundary curves of a created surface that must be 
%                      database constrained to automatically set the fitting 
%                      entities
%       AR           -> Aspect ratio of the desired cells
%       nodesTEz     -> Nodes to discretize finite trailing edge
%       advanced     -> Structure containing options for unstructured mesh only
%                      -Fields:
%                      --IsoCellType              -> 'Triangle'||'TriangleQuad'
%                      --Algorithm                -> 'Delaunay (Default)'||
%                                                   'AdvancingFront'||
%                                                   'AdvancingFrontOrtho'||
%                                                   'ThinSurfaceInterpolation'
%                      --EdgeMinimumLength        -> Number||
%                                                   'Automatic (Default)'||
%                                                   'Boundary'||
%                                                   'TRexBoundary'||
%                                                   'NotApplied'
%                      --EdgeMaximumLength        -> Number||
%                                                   'Automatic (Default)'||
%                                                   'Boundary'||
%                                                   'TRexBoundary'||
%                                                   'NotApplied'
%                      --NormalMaximumDeviation   -> Number (Default 0)
%                      --NormalSurfaceDeviation   -> Number (Default 0)
%                      --QuadMaximumIncludedAngle -> Number (Default 150)
%                      --QuadMaximumWarpAngle     -> Number (Default 30)
%       FF           -> Structure containing options for defining the farfield
%                      -Fields:
%                      --Type -> 'Box' : Box shape farfield
%                                'Cylinder' : Cylindrical shape farfield
%                                'Sphere' : Spherical shape farfield
%                      --Cell -> 'Unstructured' : Unstructured elements in the FF
%                                'Voxel' : Voxel elements in the FF
%                      --Dims -> They depend on the 'optionsFF' specified
%                                if 'Box'      -> 3x2 Matrix : [x_positive x_negative;
%                                                               y_positive y_negative;
%                                                               z_positive z_negative]
%                                if 'Cylinder' -> 3x2 Matrix : [L_positive L_negative;
%                                                               Radius_1   Radius_2;
%                                                               Axis_alignment    0]
%                                if 'Sphere'   -> 1x1 Variable : Radius
%      BLoptions     -> Structure containing options for the boundary layer
%                      -Fields:
%                      --maxLayers        -> This default is the maximum number of 
%                      T-Rex layers of an unstructured block when it is created : 
%                      Number (Default 0)
%                      --fullLayers       -> This default is the minimum number of 
%                      fully structured T-Rex layers of an unstructured block when 
%                      it is created : Number (Default 0)
%                      --growthRate       -> This default is the growth rate of 
%                      T-Rex layers of an unstructured block when it is created : 
%                      Number (Default 1.2) 
%                      --push             -> This default is the flag for pushing 
%                      T-Rex attributes onto the connectors and domains of an 
%                      unstructured block when a it is created : 'true' or 'false' 
%                      (Default false)
%                      --solverAttribute  -> Sets the named unstructured solver 
%                      attribute : 'TetPyramidPrismHex'||'AllAndReducePyramids'||
%                                  'AllAndConvertWallDoms'||'LegacyTetPyramidPrismHex'||
%                                  'TetPyramid'
%                      --firstLayerHeight  -> Specifies the BL thickness of the
%                      first cell : Number (User specified)
%      ISOoptions    -> -minLength -> Minimum length of the isotropic tets :
%                       Number||'Boundary' (Default Boundary)
%                       -maxLength -> Maximum length of the isotropic tets :
%                       Number||'Boundary' (Default Boundary)
%      sources       -> Cell array {1xN_source}
%                      Each of the elements of the cell array is a structure
%                      defining the source.
%                      -Fields:
%                      --Type        -> 'Box' : Box shape source
%                                       'Cylinder' : Cylindrical shape source
%                                       'Sphere' : Spherical shape source
%                      --Startface   -> 'x' : Source start and end faces are defined
%                                      following the unitary vector i (x axis)
%                                       'y' : Source start and end faces are defined
%                                      following the unitary vector j (y axis)
%                                       'z' : Source start and end faces are defined
%                                      following the unitary vector k (z axis)
%                      --Center      -> Specify the COG of the source in cartesian coords
%                                       1x3 Vector : [x_center y_center z_center]
%                       --Dims        -> They depend on the 'optionsSrc' specified
%                                       if 'Box'      -> 1x3 Vector : 
%                                                       [length height width]
%                                       if 'Cylinder' -> 1x3 Vector : 
%                                                       [radius top_Radius length]
%                                                       Note that radius and top_Radius
%                                                       are the same for cylinder and 
%                                                       different if defining a cone
%                                       if 'Sphere'   -> 1x3 Vector : 
%                                                       [radius base_Angle top_Angle]
%                                                       Note that base_Angle -> 
%                                                       [0,90] deg and top_Angle -> 
%                                                       [90,180] deg and are used to 
%                                                       define cut sphere sources
%                      --Distr       -> 'Constant': The begin values will be used 
%                                      throughout the source
%                                       'Parametric' : The begin values will be used at
%                                      the minimum parametric limits and the end values 
%                                      will be used at the maximum parametric limits of
%                                      the source
%                                       'AxisToPerimeter : The begin values will be 
%                                      used along the axis and the end values will be 
%                                      used at the perimeter of the source
%                                       'CenterToPerimeter : The begin values will be
%                                      used at the center and the end values will be 
%                                      used at the perimeter of the source 
%                     --Distrval     -> 2x2 Matrix: [spacing_initial spacing_final;
%                                                    decay_initial   decay_final]
%                                      Note that for the 'constant' case only the 
%                                      initial values are read
%      surfDom       -> Structure containing specification for surface mesh in terms of
%                      'structured' or 'unstructured'.
%                      Fields:
%                      -wings        -> Specify the structure on the upper and lower
%                                       wings
%                      -cant         -> Specify the structure on the cant regions
%                      -winglets     -> Specify the structure on the winglets
%       distribution -> Connector grid point distribution function is
%                      specified : 'MRQS' | 'Tanh'
%       refineVal    -> Minimum node spacing in the connector
%       N_cant       -> Number of nodes in the cant region

%       fixedPoint   -> Vector (1x3) defining the fixed point for the rotation 
%                      [xRot yRot zRot] 
%       rotAxis      -> if 'x' : The rotation axis is x
%                       if 'y' : The rotation axis is y
%                       if 'z' : The rotation axis is z
%       rotAngle     -> Angle to be rotated [deg]
%       exportPath   -> String that specifies where the mesh is exported to
%       savePath     -> String that specifies the path for the file to be saved
%       saveName     -> String that specifies the name of the file to be saved
function BoxWingGlyph(name,Slices_Wing,b,R,h,s,c_r1,c_t1,c_r2,c_t2,N,N_bl,N_bu,N_bw,el_size,tolBound,tolThres,AR,nodesTEz,advanced,FF,symmetry,BLoptions,ISOoptions,sources,surfDom,distribution,refineVal,N_cant,fixedPoint,rotAxis,rotAngle,exportPath,savePath,saveName)
addpath("sr_bwg\sr_aux\")
addpath("sr_bwg\sr_pw\")
% File
fileName = fopen(name,'w');
fprintf(fileName,'package require PWI_Glyph 5.18.5\n');

% Lines and sections
[N_lines,~,N_slices] = size(Slices_Wing);

% IDs
point_ID = zeros(N_slices,N_lines-1);
line_ID = zeros(N_slices,6);
surface_ID = zeros(N_slices-1,3);

% Counter Limits
counter_points = 0;
LE = N_lines/2;
counter_lines = 10^(ceil(log10(N_slices))*ceil(log10(N_lines-1)));
counter_surfaces = 10^(ceil(log10(N_slices))*ceil(log10(N_lines-1))+1);
counter_farfield = 10^(ceil(log10(N_slices))*ceil(log10(N_lines-1))+2);
counter_source = 10^(ceil(log10(N_slices))*ceil(log10(N_lines-1))+3);

% Surface Mesh Options
[N_LEls, N_LEus, N_AF, N_R, N_WL] = computeNumberOfElements(el_size,b,R,h,s,c_r1,c_t1,c_r2,c_t2,N,N_bl,N_bu,N_bw,AR);
nodesAirfoil = N_AF;

% Farfield
if isequal(symmetry,'off')
    dim = [s b h];
    arrangement = 2*ones(1,3);
    [~,posMax] = max(dim); arrangement(posMax) = 1;
    [~,posMin] = min(dim); arrangement(posMin) = 3;
else
    dim = [s b/2 h];
    arrangement = 2*ones(1,3);
    [~,posMax] = max(dim); arrangement(posMax) = 1;
    [~,posMin] = min(dim); arrangement(posMin) = 3;
end

% Sources
sources_size = length(sources);

%% Computations
for i = 1:N_slices
    for j = 1:N_lines-1
        counter_points = counter_points+1;
        point_ID(i,j) = counter_points;
        pointCommand = pwPoints(counter_points,Slices_Wing(j,1,i),Slices_Wing(j,2,i),Slices_Wing(j,3,i));
        fprintf(fileName,pointCommand);
    end
    if (i < N_bl+3)||(i > 2*N+N_bl+N_bw+1) %((i > 16)&&(i < 21))
            domainType = surfDom.wings;
            nodesLE = N_LEls;
            nodesTE1y = N_LEls;
            nodesTE2y = N_LEls;
    elseif ((i > N+N_bl+1)&&(i < N+N_bl+N_bw+3))
            domainType = surfDom.winglet;
            nodesLE = N_WL;
            nodesTE1y = N_WL;
            nodesTE2y = N_WL;
    else
            domainType = surfDom.cant;
            nodesLE = N_cant;
            nodesTE1y = N_cant;
            nodesTE2y = N_cant;
    end
    % Exterior curve
    counter_lines = counter_lines+1;
    line_ID(i,1) = counter_lines;
    curveCommand = pwCurve(counter_lines,point_ID(i,1:LE));
    connectorCommand = pwConnector(line_ID(i,1),nodesAirfoil,line_ID(i,1),'dimension',distribution,refineVal);
    fprintf(fileName,curveCommand);
    fprintf(fileName,connectorCommand);
    % Interior curve
    counter_lines = counter_lines+1;
    line_ID(i,2) = counter_lines;
    curveCommand = pwCurve(counter_lines,point_ID(i,LE:end));
    connectorCommand = pwConnector(line_ID(i,2),nodesAirfoil,line_ID(i,2),'dimension',distribution,refineVal);
    fprintf(fileName,curveCommand);
    fprintf(fileName,connectorCommand);
    % TE in z direction
    counter_lines = counter_lines+1;
    line_ID(i,3) = counter_lines;
    curveCommand = pwCurve(counter_lines,[point_ID(i,end) point_ID(i,1)]);
    connectorCommand = pwConnector(line_ID(i,3),nodesTEz,line_ID(i,3),'dimension');
    fprintf(fileName,curveCommand);
    fprintf(fileName,connectorCommand);
    
    if i <= 1
        line_ID(i,4:6) = NaN;
    else
        % LE
        counter_lines = counter_lines+1;
        line_ID(i,4) = counter_lines;
        curveCommand = pwCurve(counter_lines,[point_ID(i-1,LE) point_ID(i,LE)]);
        connectorCommand = pwConnector(line_ID(i,4),nodesLE,line_ID(i,4),'dimension');
        fprintf(fileName,curveCommand);
        fprintf(fileName,connectorCommand);
        % TE1 in y direction (lower)
        counter_lines = counter_lines+1;
        line_ID(i,5) = counter_lines;
        curveCommand = pwCurve(counter_lines,[point_ID(i-1,1) point_ID(i,1)]);
        connectorCommand = pwConnector(line_ID(i,5),nodesTE1y,line_ID(i,5),'dimension');
        fprintf(fileName,curveCommand);
        fprintf(fileName,connectorCommand);
        % TE2 in y direction (upper)
        counter_lines = counter_lines+1;
        line_ID(i,6) = counter_lines;
        curveCommand = pwCurve(counter_lines,[point_ID(i-1,end) point_ID(i,end)]);
        connectorCommand = pwConnector(line_ID(i,6),nodesTE2y,line_ID(i,6),'dimension');
        fprintf(fileName,curveCommand);
        fprintf(fileName,connectorCommand);
        %Surfaces
        % Exterior surface
        counter_surfaces = counter_surfaces+1;
        surface_ID(i-1,1) = counter_surfaces;
        surfaceCommand = pwSurface(counter_surfaces,[line_ID(i-1,1) line_ID(i,4) line_ID(i,1) line_ID(i,5)],tolBound,tolThres);
        domainCommand = pwDomain(counter_surfaces,[line_ID(i-1,1) line_ID(i,4) line_ID(i,1) line_ID(i,5)],domainType,advanced);
        fprintf(fileName,surfaceCommand);
        fprintf(fileName,domainCommand);
        % Interior surface
        counter_surfaces = counter_surfaces+1;
        surfaceCommand = pwSurface(counter_surfaces,[line_ID(i-1,2) line_ID(i,4) line_ID(i,2) line_ID(i,6)],tolBound,tolThres);
        domainCommand = pwDomain(counter_surfaces,[line_ID(i-1,2) line_ID(i,4) line_ID(i,2) line_ID(i,6)],domainType,advanced);
        fprintf(fileName,surfaceCommand);
        fprintf(fileName,domainCommand);
        surface_ID(i-1,2) = counter_surfaces;
        % TE surface
        domainType = 'structured';
        counter_surfaces = counter_surfaces+1;
        surfaceCommand = pwSurface(counter_surfaces,[line_ID(i-1,3) line_ID(i,6) line_ID(i,3) line_ID(i,5)],tolBound,tolThres);
        domainCommand = pwDomain(counter_surfaces,[line_ID(i-1,3) line_ID(i,6) line_ID(i,3) line_ID(i,5)],domainType);
        fprintf(fileName,surfaceCommand);
        fprintf(fileName,domainCommand);
        surface_ID(i-1,3) = counter_surfaces;
    end
    counter_points = counterUpdate(N_lines-1,counter_points);
    counter_lines = counterUpdate(6,counter_lines);
    counter_surfaces = counterUpdate(3,counter_surfaces);
end

% Angle of attack rotation
rotationCommand = pwRotation(fixedPoint,rotAxis,rotAngle,'allDataBaseandGrid',NaN,[point_ID(:);line_ID(:);surface_ID(:)],line_ID(:),surface_ID(:));
fprintf(fileName,rotationCommand);

% Symmetry copy
if isequal(symmetry,'off')
    mirrorCommand = pwMirror([point_ID(:);line_ID(:);surface_ID(:)],line_ID(:),surface_ID(:),'y');
    fprintf(fileName,mirrorCommand);
    domainNumber = 2*length(surface_ID(:));
else
    domainNumber = length(surface_ID(:));
end

% Source in between wings
sourceCommand = pwSource(counter_source,sources{1}.Type,sources{1}.StartFace,sources{1}.Center,sources{1}.Dims,sources{1}.Distr,sources{1}.DistrVal);
fprintf(fileName,sourceCommand);
rotationCommand = pwRotation(fixedPoint,rotAxis,rotAngle,'allSource',NaN,NaN,NaN,NaN);
fprintf(fileName,rotationCommand);

% Farfield generation
farfieldCommand = pwFarfield(counter_farfield,0,FF.Type,FF.TypeCell,FF.Dims,arrangement);
fprintf(fileName,farfieldCommand);

% Rest of sources (refinements)
if sources_size > 1
    for i = 2:sources_size
        counter_source = counter_source+1;
        sourceCommand = pwSource(counter_source,sources{i}.Type,sources{i}.StartFace,sources{i}.Center,sources{i}.Dims,sources{i}.Distr,sources{i}.DistrVal);
        fprintf(fileName,sourceCommand);
        rotationCommand = pwRotation(fixedPoint,rotAxis,rotAngle/2,'otherSource',counter_source,NaN,NaN,NaN);
        fprintf(fileName,rotationCommand);
    end
end

% Volumetric mesh generation
volumeMeshCommand = pwMesh(counter_farfield, BLoptions, ISOoptions);
fprintf(fileName,volumeMeshCommand);

% Export and save commands
exportCommand = pwExportOPF(counter_farfield,domainNumber,'3D',exportPath);
fprintf(fileName,exportCommand);
saveandfinishCommand = pwSaveAndFinish(savePath,saveName);
fprintf(fileName,saveandfinishCommand);
fclose(fileName);
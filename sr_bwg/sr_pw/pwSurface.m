% This function aims to generate single patched surface from a list of 
% curves for Pointwise language. Note that curves do not need to intersect,
% and the intersect tolerance will be specified by tolThres.
%       Author  : Gabriel Buendia
%       Version : 1
% Inputs:
%       ID       -> Pointwise ID desired database number
%       list     -> List of curves IDs that enclose the surface
%       tolBound -> This attribute is the tolerance to use when splitting 
%                   and joining the curves specified as the boundaries, for
%                   detecting end to end connections of the boundaries, and 
%                   for creating the surface from the boundaries
%       tolThres -> This attribute specifies the percentage of length of 
%                   the boundary curves of a created surface that must be 
%                   database constrained to automatically set the fitting 
%                   entities
% Outputs:
%       curveCommand -> String containing the command to write the surface
function surfaceCommand = pwSurface(ID,list,tolBound,tolThres)
size = length(list);
string = append('  set _DB(',num2str(ID),') [$fitter createPatch ');
for i = 1:size
    if i < size
        auxstr = append('[list $_DB(',num2str(list(i)),')] '); 
    else
        auxstr = append('[list $_DB(',num2str(list(i)),')]]\n');
    end
    string = append(string,auxstr);
end
surfaceCommand = append(['set fitter [pw::Application begin SurfaceFit]\n' ...
    '  $fitter setBoundaryTolerance ',num2str(tolBound),'\n' ...
    '  $fitter setFitEntitiesThreshold ',num2str(tolThres),'\n' ...
    string, ...
    '  $fitter run 0\n' ...
    '$fitter end\n']);
end
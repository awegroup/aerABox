% This function aims to generate a domain (2D mesh) from a list of 
% connectors for Pointwise language
%       Author  : Gabriel Buendia
%       Version : 1
% Inputs:
%       ID       -> Pointwise ID desired domain number
%       list     -> List of connector IDs enclosing the domain
%       options  -> 'structured' : Structured domain
%                   'unstructured' : Unstructured domain
%       advanced -> Structure containing options for unstructured mesh only
%                   Fields:
%                   -IsoCellType              -> 'Triangle'||'TriangleQuad'
%                   -Algorithm                -> 'Delaunay (Default)'||
%                                                'AdvancingFront'||
%                                                'AdvancingFrontOrtho'||
%                                                'ThinSurfaceInterpolation'
%                   -EdgeMinimumLength        -> Number||
%                                                'Automatic (Default)'||
%                                                'Boundary'||
%                                                'TRexBoundary'||
%                                                'NotApplied'
%                   -EdgeMaximumLength        -> Number||
%                                                'Automatic (Default)'||
%                                                'Boundary'||
%                                                'TRexBoundary'||
%                                                'NotApplied'
%                   -NormalMaximumDeviation   -> Number (Default 0)
%                   -NormalSurfaceDeviation   -> Number (Default 0)
%                   -QuadMaximumIncludedAngle -> Number (Default 150)
%                   -QuadMaximumWarpAngle     -> Number (Default 30)
% Outputs:
%       domainCommand -> String containing the command to write the
%                        domain
function domainCommand = pwDomain(ID,list,options,advanced)
size = length(list);
string = '[list ';
for i = 1:size
    if i < size
        auxstr = append('$_CN(',num2str(list(i)),') '); 
    else
        auxstr = append('$_CN(',num2str(list(i)),')]');
    end
    string = append(string,auxstr);
end
switch options
    case 'structured'
        domainCommand = append('set _DM(',num2str(ID),') [pw::DomainStructured createFromConnectors ',string,']\n');
    case 'unstructured'
        domainCommand = append('set _DM(',num2str(ID),') [pw::DomainUnstructured createFromConnectors ',string,']\n');
        if nargin > 3
            IsoCellType = advanced.IsoCellType;

            if sum(strcmp(fieldnames(advanced), 'Algorithm')) == 1
                Algorithm = advanced.Algorithm;
            else
                Algorithm = 'Delaunay';
            end

            if sum(strcmp(fieldnames(advanced), 'EdgeMinimumLength')) == 1
                EdgeMinimumLength = advanced.EdgeMinimumLength;
            else
                EdgeMinimumLength = 'Automatic';
            end

            if sum(strcmp(fieldnames(advanced), 'EdgeMaximumLength')) == 1 
                EdgeMaximumLength = advanced.EdgeMaximumLength;
            else
                EdgeMaximumLength = 'Automatic';
            end

            if sum(strcmp(fieldnames(advanced), 'NormalMaximumDeviation')) == 1 
                NormalMaximumDeviation = advanced.NormalMaximumDeviation;
            else
                NormalMaximumDeviation = 0;
            end

            if sum(strcmp(fieldnames(advanced), 'SurfaceMaximumDeviation')) == 1 
                SurfaceMaximumDeviation = advanced.NormalMaximumDeviation;
            else
                SurfaceMaximumDeviation = 0;
            end
                        
            if isequal(IsoCellType,'TriangleQuad')
                if sum(strcmp(fieldnames(advanced), 'QuadMaximumWarpAngle')) == 1 
                    QuadMaximumIncludedAngle = advanced.QuadMaximumIncludedAngle;
                else
                    QuadMaximumIncludedAngle = 150;
                end

                if sum(strcmp(fieldnames(advanced), 'QuadMaximumWarpAngle')) == 1 
                    QuadMaximumWarpAngle = advanced.QuadMaximumWarpAngle;
                else
                    QuadMaximumWarpAngle = 30;
                end
                domainCommand = append([domainCommand,...
                    'set _TMP(mode_1) [pw::Application begin UnstructuredSolver [list $_DM(',num2str(ID),')]]\n'...
                    '   $_DM(',num2str(ID),') setUnstructuredSolverAttribute IsoCellType ',IsoCellType,'\n'...
                    '   $_DM(',num2str(ID),') setUnstructuredSolverAttribute Algorithm ',Algorithm,'\n'...
                    '   $_DM(',num2str(ID),') setUnstructuredSolverAttribute EdgeMinimumLength ',num2str(EdgeMinimumLength),'\n'...
                    '   $_DM(',num2str(ID),') setUnstructuredSolverAttribute EdgeMaximumLength ',num2str(EdgeMaximumLength),'\n'...
                    '   $_DM(',num2str(ID),') setUnstructuredSolverAttribute NormalMaximumDeviation ',num2str(NormalMaximumDeviation),'\n'...
                    '   $_DM(',num2str(ID),') setUnstructuredSolverAttribute SurfaceMaximumDeviation ',num2str(SurfaceMaximumDeviation),'\n'...
                    '   $_DM(',num2str(ID),') setUnstructuredSolverAttribute QuadMaximumIncludedAngle ',num2str(QuadMaximumIncludedAngle),'\n'...
                    '   $_DM(',num2str(ID),') setUnstructuredSolverAttribute QuadMaximumWarpAngle ',num2str(QuadMaximumWarpAngle),'\n'...
                    '   $_TMP(mode_1) run Initialize\n' ...
                    '$_TMP(mode_1) end\n']);
            else
                domainCommand = append([domainCommand,...
                    'set _TMP(mode_1) [pw::Application begin UnstructuredSolver [list $_DM(',num2str(ID),')]]\n'...
                    '   $_DM(',num2str(ID),') setUnstructuredSolverAttribute IsoCellType ',IsoCellType,'\n'...
                    '   $_DM(',num2str(ID),') setUnstructuredSolverAttribute Algorithm ',Algorithm,'\n'...
                    '   $_DM(',num2str(ID),') setUnstructuredSolverAttribute EdgeMinimumLength ',num2str(EdgeMinimumLength),'\n'...
                    '   $_DM(',num2str(ID),') setUnstructuredSolverAttribute EdgeMaximumLength ',num2str(EdgeMaximumLength),'\n'...
                    '   $_DM(',num2str(ID),') setUnstructuredSolverAttribute NormalMaximumDeviation ',num2str(NormalMaximumDeviation),'\n'...
                    '   $_DM(',num2str(ID),') setUnstructuredSolverAttribute SurfaceMaximumDeviation ',num2str(SurfaceMaximumDeviation),'\n'...
                    '   $_TMP(mode_1) run Initialize\n' ...
                    '$_TMP(mode_1) end\n']);
            end
        end
end
end

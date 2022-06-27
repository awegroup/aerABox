% This function aims to generate a connector with nodes on a curve (1D mesh)
% for Pointwise language
%       Author  : Gabriel Buendia
%       Version : 1
% Inputs:
%       ID           -> Pointwise ID desired connector number
%       dim          -> Expected nodes on the connector or spacing
%       ID_DB        -> Pointwise ID of the database curve
%       spacing      -> 'dimension' : Nodes are created specifying a dimension
%                       'separation' : Nodes are created specifying an average spacing
%       distribution -> Connector grid point distribution function is
%                       specified : 'MRQS' | 'Tanh'
%       refineVal    -> Minimum node spacing in the connector
% Outputs:
%       connectorCommand -> String containing the command to write the
%                           connector
function connectorCommand = pwConnector(ID,dim,ID_DB,spacing,distribution,refineVal)
switch spacing
    case 'dimension'
    dimension = append('pw::Connector setDefault Dimension ',num2str(dim),'\n  ');
    case 'separation'
    dimension = append(['pw::Connector setCalculateDimensionMethod Spacing\n' ...
        'pw::Connector setCalculateDimensionSpacing '],num2str(dim),'\n  ');
end

connectorCommand = append([dimension, ...
    'set _CN(',num2str(ID),') [pw::Connector createOnDatabase $_DB(',num2str(ID_DB),')]\n', ...
    '$_CN(',num2str(ID),') replaceDistribution 1 [pw::DistributionTanh create]\n']);
if nargin > 4
    distributionCommand = append(['set _TMP(mode_1) [pw::Application begin Modify [list $_CN(',num2str(ID),')]]\n', ...
                               '  pw::Connector swapDistribution ',distribution,' [list [list $_CN(',num2str(ID),') 1]]\n', ...
                               '  [[$_CN(',num2str(ID),') getDistribution 1] getEndSpacing] setValue ',num2str(refineVal),'\n', ... 
                               '  [[$_CN(',num2str(ID),') getDistribution 1] getBeginSpacing] setValue ',num2str(refineVal),'\n', ... '
                               '$_TMP(mode_1) end\n', ...
                               'unset _TMP(mode_1)\n']);
    connectorCommand = append(connectorCommand,distributionCommand);
end
end
% This function aims to generate a mirror copy of all database, connectors 
% and domains selected with respect to certain cartesian axis for
% Pointwise language
%       Author  : Gabriel Buendia
%       Version : 1
% Inputs:
%       listDB -> Defines the list of database entities (IDs) to be copied
%       listCN -> Defines the list of conectors (IDs) to be copied
%       listDM -> Defines the list of domains (IDs) to be copied
%       axis   -> if 'x' : The reflection plane is yz
%                 if 'y' : The reflection plane is xz
%                 if 'z' : The reflection plane is xy
% Outputs:
%       mirrorCommand -> String containing the command to write the mirror copy
function mirrorCommand = pwMirror(listDB,listCN,listDM,axis)
switch axis
    case 'x'
        axis = [1 0 0];
    case 'y'
        axis = [0 1 0];
    case 'z'
        axis = [0 0 1];
end
sizeDB = length(listDB);
sizeCN = length(listCN);
sizeDM = length(listDM);
string = 'pw::Application setClipboard [list ';
for i = 1:sizeDB
    if ~isnan(listDB(i))
        string = append(string,'$_DB(',num2str(listDB(i)),') ');
    end
end
for i = 1:sizeCN
    if ~isnan(listCN(i))
        string = append(string,'$_CN(',num2str(listCN(i)),') ');
    end
end
for i = 1:sizeDM
    if ~isnan(listDM(i))
        if i < sizeDM
            string = append(string,'$_DM(',num2str(listDM(i)),') ');
        else
            string = append(string,'$_DM(',num2str(listDM(i)),')]\n');
        end
    end
end
mirrorCommand = append(['pw::Application clearClipboard\n', ...
    string, ...
    'pw::Application markUndoLevel Copy\n', ...
    'set _TMP(mode_1) [pw::Application begin Paste]\n', ...
    '  set _TMP(PW_1) [$_TMP(mode_1) getEntities]\n', ...
    '  set _TMP(mode_2) [pw::Application begin Modify $_TMP(PW_1)]\n', ...
    '    pw::Entity transform [pwu::Transform mirroring {',num2str(axis(1)),' ',num2str(axis(2)),' ',num2str(axis(3)),'} 0] [$_TMP(mode_2) getEntities]\n', ...
    '  $_TMP(mode_2) end\n', ...
    '  unset _TMP(mode_2)\n', ...
    '$_TMP(mode_1) end\n', ...
    'unset _TMP(mode_1)\n']);
end
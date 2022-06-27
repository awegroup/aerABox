% This function aims to generate a rotation of all entities selected with 
% respect to certain fixed point and cartesian axis for Pointwise language
%       Author  : Gabriel Buendia
%       Version : 1
% Inputs:
%       fixedPoint -> Vector (1x3) defining the fixed point for the rotation 
%                     [xRot yRot zRot] 
%       rotAxis    -> if 'x' : The rotation axis is x
%                     if 'y' : The rotation axis is y
%                     if 'z' : The rotation axis is z
%       rotAngle   -> Angle to be rotated [deg]
%       option     -> String defining the entities to be rotated
%                     -'allDataBaseandGrid': Rotates all database and grid entities
%                     (connectors and domains) specified in listDB, listCN and list DM
%                     -'allDataBase': Rotates all database entities created
%                     so far withouth list specification
%                     -'allGrid': Rotates all grid entities created so far withouth 
%                     list specification
%                     -'allSource': Rotates all source entities created so far withouth 
%                     list specification
%                     -'otherDataBase': Rotates all database entities
%                     spcified in the list listO
%                     -'otherConnector': Rotates all connector entities
%                     spcified in the list listO
%                     -'otherDomain': Rotates all domain entities
%                     spcified in the list listO
%                     -'otherSource': Rotates all source entities
%                     spcified in the list listO
%       listO      -> Defines the list of entities (IDs) specified in option 'other...'
%                     to be rotated
%       listDB     -> Defines the list of database entities (IDs) to be rotated
%       listCN     -> Defines the list of conectors (IDs) to be rotated
%       listDM     -> Defines the list of domains (IDs) to be rotated

% Outputs:
%       rotationCommand -> String containing the command to write the rotation
function rotationCommand = pwRotation(fixedPoint,rotAxis,rotAngle,option,listO,listDB,listCN,listDM)
switch rotAxis
    case 'x'
        rotAxis = [1 0 0];
    case 'y'
        rotAxis = [0 1 0];
    case 'z'
        rotAxis = [0 0 1];
end
switch option 
    case 'allDataBaseandGrid'
    sizeDB = length(listDB);
    sizeCN = length(listCN);
    sizeDM = length(listDM);
    string = '  set _TMP(mode_1) [pw::Application begin Modify [list ';
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
                string = append(string,'$_DM(',num2str(listDM(i)),')]]\n');
            end
        end
    end
    string2 = '';
    case 'allDataBase'
        string = append(['set _TMP(PW_1) [pw::Database getAll]\n' ...
                 '  set _TMP(mode_1) [pw::Application begin Modify $_TMP(PW_1)]\n']);
        string2 = 'unset _TMP(PW_1)\n';
    case 'allGrid'
        string = append(['set _TMP(PW_1) [pw::Grid getAll]\n' ...
                 '  set _TMP(mode_1) [pw::Application begin Modify $_TMP(PW_1)]\n']);
        string2 = 'unset _TMP(PW_1)\n';
    case 'allSource'
        string = append(['set _TMP(PW_1) [pw::Source getAll]\n' ...
                 '  set _TMP(mode_1) [pw::Application begin Modify $_TMP(PW_1)]\n']);
        string2 = 'unset _TMP(PW_1)\n';
    case 'otherDataBase'
        sizeO = length(listO);
        string = '  set _TMP(mode_1) [pw::Application begin Modify [list ';
        for i = 1:sizeO
            if ~isnan(listO(i))
                if i < sizeO
                    string = append(string,'$_DB(',num2str(listO(i)),') ');
                else
                    string = append(string,'$_DB(',num2str(listO(i)),')]]\n');
                end
            end
        end
        string2 = '';
    case 'otherConnector'
        sizeO = length(listO);
        string = '  set _TMP(mode_1) [pw::Application begin Modify [list ';
        for i = 1:sizeO
            if ~isnan(listO(i))
                if i < sizeO
                    string = append(string,'$_CN(',num2str(listO(i)),') ');
                else
                    string = append(string,'$_CN(',num2str(listO(i)),')]]\n');
                end
            end
        end
        string2 = '';
    case 'otherDomain'
        sizeO = length(listO);
        string = '  set _TMP(mode_1) [pw::Application begin Modify [list ';
        for i = 1:sizeO
            if ~isnan(listO(i))
                if i < sizeO
                    string = append(string,'$_DM(',num2str(listO(i)),') ');
                else
                    string = append(string,'$_DM(',num2str(listO(i)),')]]\n');
                end
            end
        end
        string2 = '';
    case 'otherSource'
        sizeO = length(listO);
        string = '  set _TMP(mode_1) [pw::Application begin Modify [list ';
        for i = 1:sizeO
            if ~isnan(listO(i))
                if i < sizeO
                    string = append(string,'$_SR(',num2str(listO(i)),') ');
                else
                    string = append(string,'$_SR(',num2str(listO(i)),')]]\n');
                end
            end
        end
        string2 = '';
end
rotationCommand = append([string, ...
    '    pw::Entity transform [pwu::Transform rotation -anchor {',num2str(fixedPoint(1)),' ',num2str(fixedPoint(2)),' ',num2str(fixedPoint(3)),'} {',num2str(rotAxis(1)),' ',num2str(rotAxis(2)),' ',num2str(rotAxis(3)),'} ',num2str(rotAngle),'] [$_TMP(mode_1) getEntities]\n' ...
    '  $_TMP(mode_1) end\n' ...
    '  unset _TMP(mode_1)\n' ...
    string2]);
end
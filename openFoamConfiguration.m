clear; close all; clc;
addpath("sr_of\")
vec_alpha = 4;
for i = vec_alpha
folder = num2str(i);
if i < 0
    folder(1) = 'n';
else
    for j = 1:length(folder)
        if folder(j) == '.'
            folder(j) = '_';
        end
    end
end
% Change respectively
nameFolder = append('C_BoxWing_',folder,'_deg');
localFolder = append('F:\Final_Code\',nameFolder,'\');
meshFolder = append('F:\Final_Code\BoxWingPaper\',folder);
templateFolder = 'F:\Final_Code\EXAMPLE\';
mkdir(localFolder)

localSys = append(localFolder,'\system\');
localCon = append(localFolder,'\constant\');
local_0 = append(localFolder,'\0\');
mkdir(localFolder,'0')
mkdir(localFolder,'constant')
mkdir(localFolder,'system')
copyfile(append(templateFolder,'0\k'),local_0)
copyfile(append(templateFolder,'0\nut'),local_0)
copyfile(append(templateFolder,'0\omega'),local_0)
copyfile(append(templateFolder,'0\p'),local_0)
copyfile(append(templateFolder,'0\U'),local_0)
mkdir(localCon,'polyMesh')
copyfile(append(templateFolder,'constant\transportProperties'),localCon)
copyfile(append(templateFolder,'constant\turbulenceProperties'),localCon)
copyfile(append(templateFolder,'system\fvSchemes'),localSys)
copyfile(append(templateFolder,'system\fvSolution'),localSys)
copyfile(append(templateFolder,'system\yPlus'),localSys)
copyfile(append(templateFolder,'system\wallShearStress'),localSys)
copyfile(meshFolder,append(localCon,'polyMesh\'))

load('Flight_Conditions.mat')
solver = 'simpleFoam';
startTime = 0;
endTime = 10001;
deltaT = 1;
writeInterval = 500;
writePrecision = 16;
timePrecision = 12;
AoA = 0; % Even if the actual angle of attack is non zero, the geometry is already rotated
nut = 0;
I = 2/100;
e = 1000;
wallname = 'cylinder';
CofR = [0 0 0];
writeIntervalForces = 1;
cores = 48;

job.FolderPath = localFolder;
job.Name = nameFolder;
job.Queue = 'awep-small';
job.Hours = 72;
job.Nodes = 1;
job.Cores = cores;
job.NodeType = 'g';
job.Solver = solver;
job.FinalTime = endTime;
job.Path = append('/home/gbuendiavela/',nameFolder);

ofControlDict(localSys,solver,startTime,endTime,deltaT,writeInterval,writePrecision,timePrecision) 
ofIncludeDict(localSys,Re,AoA,nu,nut,p,rho,I,e,MAC,airspeed) 
ofWingForces(localSys,wallname,CofR,writeIntervalForces) 
ofDecomposeParDict(localSys,job.Cores*job.Nodes)
ofClusterJob(job)
end
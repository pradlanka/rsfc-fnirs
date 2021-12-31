% This code creates the fig 1 in the manuscript which comapres AR partial
% correlation with AR correlation and Pearson's correlation in correcting
% for global spatial signal and temporal autocorrelation.

Curr_Folder = pwd; % Get path of working directory

% This code requires NIRS Brain AnalyzIR toolbox to run. The toolbox can be
% downloaded from  https://github.com/huppertt/nirs-toolbox. 


% Please download the toolbox and add to path

addpath(genpath('nirs-toolbox')); % Add toolbox to path

addpath(genpath('functions')); % Add to path folder containing important functions]


% Load the probe file & set the parameters for probe to draw
load([Curr_Folder filesep 'config' filesep 'probe.mat']);
probe.defaultdrawfcn = '10-20';
mesh=probe.getmesh;
mesh(1).transparency=0.1;
mesh(1).fiducials.Draw(:)=false;
mesh(1).fiducials.Draw([1,2,3,45, 85]) = true;
[raw , truth] = simData_connectivity_shortsep([],[1],0, probe);


%Preprocessing pipeline 
jobpre = nirs.modules.OpticalDensity;
jobpre = nirs.modules.BeerLambertLaw(jobpre);
jobpre=nirs.modules.LabelShortSeperation(jobpre);
jobpre = nirs.modules.KeepTypes(jobpre);
jobpre.types = 'hbo'; %Just HbO should be sufficient since, 
% HBR should give similar results as they were generated from distribution
% from same statisical properties.
% Remove short channels from the data and keep just HBO
jobpre_lg=nirs.modules.RemoveShortSeperations(jobpre);


% Generate the simulated ground truth. This ground truth will be used to
% evaluate the performance of several analysis pipelines

% With spatially correlated physiological noise & temporally autocorrelated
% noise
[raw1 , truth] = simData_connectivity_shortsep(truth,[1],150,probe); 
% truth contains the generated covariance matrix and raw the corresponding data

hbORest = jobpre.run(raw1);
hbORest_lg= jobpre_lg.run(raw1);
% Processed data with just long channels. Currently functions in nirs.sFC
% cannot process data with short sepration channels.

% Run connectivity  models

job = nirs.modules.Connectivity;


% Standard Pearsons' correlation without robust
job.corrfcn=@(data)nirs.sFC.corr(data,false); 
ConnStats11 = job.run(hbORest_lg);
ConnStats11(1).probe.defaultdrawfcn='10-20';
figure; ConnStats11(1).draw('R',[-1 1],'q<0.05');


 % AR correlation without robust
job.corrfcn=@(data)nirs.sFC.ar_corr(data,'10x',false);
ConnStats21 = job.run(hbORest_lg);
ConnStats21(1).probe.defaultdrawfcn='10-20';
ConnStats21(1).draw('R',[-1 1],'q<0.05')


% Partial correlation with PCA filtering of short channels with robust
% regression
job.corrfcn=@(data)partial_corr(data,'10xFs','short',false, true);
ConnStats31 = job.run(hbORest);
ConnStats31(1).probe = ConnStats21(1).probe;
ConnStats31(1).draw('R',[-1 1],'q<0.05')

    
% Plotting the ground truth on the 3D model of the probe
ConnStats0 = ConnStats11;
ConnStats0.R =  truth;
ConnStats0.draw('R',[-1 1],'q<0.05')


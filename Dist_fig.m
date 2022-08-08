% This code creates the figures 6-8 as well as fig 13 in the manuscript

Curr_Folder = pwd; % Get path of working directory

% This code requires NIRS Brain AnalyzIR toolbox to run. The toolbox can be
% downloaded from  https://github.com/huppertt/nirs-toolbox. 


% Please download the toolbox and add to path

addpath(genpath('nirs-toolbox')); % Add toolbox to path

addpath(genpath('functions')); % Add to path folder containing important functions

%Creating folder to save results
fig_dir  = [Curr_Folder filesep 'Figures'];
mkdir(fig_dir)



% Preprocessing pipeline to preprocess resting-state fNIRS data

job = nirs.modules.OpticalDensity;
job = nirs.modules.BeerLambertLaw(job);
job = nirs.modules.KeepTypes(job);
job.types = 'hbo'; % Keep only HBO as similar results should apply to HbR
job = nirs.modules.LabelShortSeperation(job);
job_lg=nirs.modules.RemoveShortSeperations(job); %Currently  sFC in nirs toolbox
% toolbox does not support short seperation channels so this step removes
% short channels


% Generate the simulated ground truth. This ground truth will be used to
% evaluate the performance of several analysis pipelines
[raw , truth] = simData_connectivity_shortsep([],[1],0); % truth contains 
% the generated covariance matrix and raw the corresponding data
hboRest = job.run(raw);
hboRest_lg= job_lg.run(raw); % Processed data with just long channels


% Setting simulation Parameters
nIters =100; % No of iterations, no of datasets to be generated
Fs = raw.Fs; %Sampling rate
Pmax = round(10*Fs); % Max model order for autoregressive (AR) models
criterion = 'BIC'; %Informtation criteria to compare model fits to find optimal model orders in AR and MVAR models 

%Extract info from data
nchan = sum(~raw.probe.link.ShortSeperation)/2; %No of long channels
nshort= sum(raw.probe.link.ShortSeperation)/2; %No of short channels
ntp = size(raw.data, 1); %No of timepoints in the data
zeroind = triu(true(size(truth)),1)&(~truth); %Indices of elements with true negatives in the data covariance  matrix
dataind = triu(true(nchan),1); % getting the upper right hand corner of the data covariance matrix.

pathsa = nnz(zeroind); %True negative paths
pathsb = 0.5*nchan*(nchan-1); % All paths



% Truth distribution
% This will later be used to compare the performance of simulation results
t = sqrt(ntp-Pmax).*((truth)./(sqrt(1-truth.^2)));
pvals = 2*nirs.math.tpvalue(-abs(t),ntp-Pmax);
truth_binp  = pvals<0.05; %Binarizing the "ground truth".
% ones in the truth_binp denote a connection between channels and 0 denotes no connection
nnzind = triu(true(size(truth_binp)),1)&(truth_binp); %Indices of elements with true positives
%in the data covariance  matrix
tp = nnz(nnzind); % No of true positives
tn = nnz(zeroind); % No of true negatives
nvals = min(tp,tn); %Often are data has less true positives than true negatives
truth_all_p = repmat(vertcat(ones(nvals,1),zeros(nvals,1)),nIters,1); % Will later be used to generate ROC
% curves which require equal no of true positives and true negatives.
truth_all_p = truth_all_p(:); % convert to a single long column

% Null distribution. This will later be used to compare the expected false
% positives to actual false positives for various analyis methods
 pvals= (1/(pathsa*nIters))*(1:(pathsa*nIters));
 int = 0.01:0.001:0.99;
 pvals_int = quantile(pvals, int);

% Data labels indicating the statisitical properties of each generated
% dataset
datalabels = {'random data','w temp corr noise','w temp & spatial corr noise',....
    'w temp & spatial corr noise & motion',};

for type=1:length(datalabels)
  mkdir([Curr_Folder filesep 'Figures' filesep datalabels{type}])
  cd([Curr_Folder filesep 'Figures' filesep datalabels{type}])
   
%% Running ths simulation
    
% Initializing variables to store values obtained from the analysis pipeline with zeros
Rdist_1a = zeros(pathsa,nIters); Rdist_2a = zeros(pathsa,nIters);
pdist_1a = zeros(pathsa,nIters); pdist_2a = zeros(pathsa,nIters); 
Rdist_11a = zeros(pathsa,nIters); Rdist_21a = zeros(pathsa,nIters);
pdist_11a = zeros(pathsa,nIters); pdist_21a = zeros(pathsa,nIters);
Rdist_3a = zeros(pathsa,nIters);  Rdist_31a = zeros(pathsa,nIters);
pdist_3a = zeros(pathsa,nIters); pdist_31a = zeros(pathsa,nIters);

Fdist_6a = zeros(pathsa,nIters); Fdist_61a = zeros(pathsa,nIters);
pdist_6a = zeros(pathsa,nIters); pdist_61a = zeros(pathsa,nIters);
Gdist_6a = zeros(pathsa,nIters); Gdist_61a = zeros(pathsa,nIters);


Rdist_1b = zeros(pathsb,nIters); Rdist_2b = zeros(pathsb,nIters);
pdist_1b = zeros(pathsb,nIters); pdist_2b = zeros(pathsb,nIters); 
Rdist_11b = zeros(pathsb,nIters); Rdist_21b = zeros(pathsb,nIters);
pdist_11b = zeros(pathsb,nIters); pdist_21b = zeros(pathsb,nIters);
Rdist_3b = zeros(pathsb,nIters); Rdist_31b = zeros(pathsb,nIters); 
pdist_3b = zeros(pathsb,nIters); pdist_31b = zeros(pathsb,nIters); 


Fdist_6b = zeros(pathsb,nIters); Fdist_61b = zeros(pathsb,nIters);
pdist_6b = zeros(pathsb,nIters); pdist_61b = zeros(pathsb,nIters);
Gdist_6b = zeros(pathsb,nIters); Gdist_61b = zeros(pathsb,nIters);


% Running several iterations
 for iter=1:nIters
   
   % Simulate the data based on statisical properites

    switch type-1
        case 1
          % With temporal autocorrelated noise 
          [raw1 , ~] = simData_connectivity_shortsep(truth,[1,0],0);
          data = job.run(raw1);
          data_lg= job_lg.run(raw1);

        case 2
           % With termporal autocorrelation & spatially correlated physiological noise
          [raw1 , ~] = simData_connectivity_shortsep(truth,[1,0],150);
          data = job.run(raw1);
          data_lg = job_lg.run(raw1);  

        case 3
          % WIth motion artifacts & spatially correlated physiological noise
          [raw1 , ~] = simData_connectivity_shortsep(truth,[1,0],150);
          raw1 = nirs.testing.simMotionArtifact(raw1); % adding motion artifacts
          data = job.run(raw1);
          data_lg = job_lg.run(raw1);

        otherwise
          % Random noise
          [raw1 , ~] = simData_connectivity_shortsep(truth,[1,0],0);      
          data = job.run(raw1);
          data_lg = job_lg.run(raw1);
          data.data = rand(ntp,nchan+nshort);
          data_lg.data = rand(ntp,nchan) ;
    end

    % Standard Pearsons' correlation with robust
    [R1, p1, dfe1] = nirs.sFC.corr(data_lg,true);
    Rdist_1a(:, iter) = R1(zeroind);
    pdist_1a(:, iter) = p1(zeroind);
    Rdist_1b(:, iter) = cat(1,R1(nnzind),R1(zeroind));
    pdist_1b(:, iter) = cat(1,p1(nnzind),p1(zeroind));

    % Standard Pearsons' correlation without robust
    [R11, p11, dfe11] = nirs.sFC.corr(data_lg,false);
    Rdist_11a(:, iter) = R11(zeroind);
    pdist_11a(:, iter) = p11(zeroind);
    Rdist_11b(:, iter) = cat(1,R11(nnzind),R11(zeroind));
    pdist_11b(:, iter) = cat(1,p11(nnzind),p11(zeroind));
    

    % AR correlation with robust
    [R2, p2,dfe2] = nirs.sFC.ar_corr(data_lg,'10xFs',true);
    Rdist_2a(:, iter) = R2(zeroind);
    pdist_2a(:, iter) = p2(zeroind);
    Rdist_2b(:, iter) = cat(1,R2(nnzind),R2(zeroind));
    pdist_2b(:, iter) = cat(1,p2(nnzind),p2(zeroind));

    % AR correlation without robust
    [R21, p21,dfe21] = nirs.sFC.ar_corr(data_lg,'10xFs',false);
    Rdist_21a(:, iter) = R21(zeroind);
    pdist_21a(:, iter) = p21(zeroind);
    Rdist_21b(:, iter) = cat(1,R21(nnzind),R21(zeroind));
    pdist_21b(:, iter) = cat(1,p21(nnzind),p21(zeroind));

    % AR Partial correlation with PCA filtering of short channels with robust
    % regression
    [R3, p3, dfe3] = partial_corr(data,'10xFs','short',true, true);
    Rdist_3a(:, iter) = R3(zeroind);
    pdist_3a(:, iter) = p3(zeroind);
    Rdist_3b(:, iter) = cat(1,R3(nnzind),R3(zeroind));
    pdist_3b(:, iter) = cat(1,p3(nnzind),p3(zeroind));
    
    % AR Partial correlation with PCA filtering of short channels without
    % robust regression
    [R31, p31, dfe31] = partial_corr(data,'10xFs','short',false, true);
    Rdist_31a(:, iter) = R31(zeroind);
    pdist_31a(:, iter) = p31(zeroind);
    Rdist_31b(:, iter) = cat(1,R31(nnzind),R31(zeroind));
    pdist_31b(:, iter) = cat(1,p31(nnzind),p31(zeroind));
    
    % Modified MVGC with zero-lag after PCA filtering of short channels
    % with robust regression
     [G6, F6, df16, df26, p6] = grangercausality(data,Pmax,'multivariate','short',true,true,true,criterion);
     Gdist_6a(:, iter) = G6(zeroind);
     Fdist_6a(:, iter) = F6(zeroind);
     pdist_6a(:, iter) = p6(zeroind);
     Gdist_6b(:, iter) = cat(1,G6(nnzind),G6(zeroind));
     Fdist_6b(:, iter) = cat(1,F6(nnzind),F6(zeroind));
     pdist_6b(:, iter) = cat(1,p6(nnzind),p6(zeroind));
     
    % Modified MVGC with zero-lag after PCA filtering of short channels
    % without robust regression
     [G61, F61, df161, df261, p61] = grangercausality(data,Pmax,'multivariate','short',true,false,true,criterion);
     Gdist_61a(:, iter) = G61(zeroind);
     Fdist_61a(:, iter) = F61(zeroind);
     pdist_61a(:, iter) = p61(zeroind);
     Gdist_61b(:, iter) = cat(1,G61(nnzind),G61(zeroind));
     Fdist_61b(:, iter) = cat(1,F61(nnzind),F61(zeroind));
     pdist_61b(:, iter) = cat(1,p61(nnzind),p61(zeroind));
  
   fprintf("Iteration completed: %d of %d\n", iter, nIters);   
 end
 
 %% Calculating the necessary metrics
 
 %  Converting the r values obtrained from anlaysis pipelines into column vector for later plotting of
 % bootstrap -p vs actual p
 Rdist_1a = Rdist_1a(:); pdist_1a = pdist_1a(:); 
 Rdist_11a = Rdist_11a(:); pdist_11a = pdist_11a(:);
 Rdist_2a = Rdist_2a(:); pdist_2a = pdist_2a(:); 
 Rdist_21a = Rdist_21a(:); pdist_21a = pdist_21a(:);
 Rdist_3a = Rdist_3a(:); pdist_3a = pdist_3a(:); 
 Rdist_31a = Rdist_31a(:); pdist_31a = pdist_31a(:);
 Gdist_6a = Gdist_6a(:); Fdist_6a = Fdist_6a(:); pdist_6a = pdist_6a(:); 
 Gdist_61a = Gdist_61a(:); Fdist_61a = Fdist_61a(:); pdist_61a = pdist_61a(:);  

 % Sorting the R values obtrained in ascending order 
 [R_sorted_1, indxR1]= sort(Rdist_1a, 'ascend'); [R_sorted_11, indxR11]= sort(Rdist_11a, 'ascend'); 
 [R_sorted_2, indxR2]= sort(Rdist_2a, 'ascend'); [R_sorted_21, indxR21]= sort(Rdist_21a, 'ascend'); 
 [R_sorted_3, indxR3]= sort(Rdist_3a, 'ascend'); [R_sorted_31, indxR31]= sort(Rdist_31a, 'ascend'); 
 [G_sorted_6, indxG6]= sort(Gdist_6a, 'ascend'); [G_sorted_61, indxG61]= sort(Gdist_61a, 'ascend'); 
 
 
 % Calculations for generating ROC curves and calculating AUC
  indx1 =  randperm(tp);
  indx2 =  randperm(tn);
 
% Since we have less true positives, we randomly sample a subset of true
% negatives to plot ROC curves which require equal number of true positives
% & negatives
 pdist_1b = cat(1,pdist_1b(indx1(1:nvals),:),pdist_1b(indx2(1:nvals),:));
 pdist_11b = cat(1,pdist_11b(indx1(1:nvals),:),pdist_11b(indx2(1:nvals),:)); 
 pdist_2b = cat(1,pdist_2b(indx1(1:nvals),:),pdist_2b(indx2(1:nvals),:));
 pdist_21b = cat(1,pdist_21b(indx1(1:nvals),:),pdist_21b(indx2(1:nvals),:));   
 pdist_3b = cat(1,pdist_3b(indx1(1:nvals),:),pdist_3b(indx2(1:nvals),:));
 pdist_31b = cat(1,pdist_31b(indx1(1:nvals),:),pdist_31b(indx2(1:nvals),:)); 
 pdist_6b = cat(1,pdist_6b(indx1(1:nvals),:),pdist_6b(indx2(1:nvals),:));
 pdist_61b = cat(1,pdist_61b(indx1(1:nvals),:),pdist_61b(indx2(1:nvals),:)); 
 
 pdist_1b = pdist_1b(:); pdist_11b = pdist_11b(:); 
 pdist_2b = pdist_2b(:); pdist_21b = pdist_21b(:);    
 pdist_3b = pdist_3b(:); pdist_31b = pdist_31b(:); 
 pdist_6b = pdist_6b(:); pdist_61b = pdist_61b(:); 
 
 % ROC curve analysis
[tp_p1,fp_p1,th_p1] = nirs.testing.roc(truth_all_p, pdist_1b); 
[tp_p11,fp_p11,th_p11] = nirs.testing.roc(truth_all_p, pdist_11b); 
[tp_p2,fp_p2,th_p2] = nirs.testing.roc(truth_all_p, pdist_2b);
[tp_p21,fp_p21,th_p21] = nirs.testing.roc(truth_all_p, pdist_21b);
[tp_p3,fp_p3,th_p3] = nirs.testing.roc(truth_all_p, pdist_3b);
[tp_p31,fp_p31,th_p31] = nirs.testing.roc(truth_all_p, pdist_31b); 
[tp_p6,fp_p6,th_p6] = nirs.testing.roc(truth_all_p, pdist_6b);
[tp_p61,fp_p61,th_p61] = nirs.testing.roc(truth_all_p, pdist_61b); 

% AUC metrics 
AUC_p1 = trapz(fp_p1, tp_p1); AUC_p11 = trapz(fp_p11, tp_p11); 
AUC_p2 = trapz(fp_p2, tp_p2); AUC_p21 = trapz(fp_p21, tp_p21); 
AUC_p3 = trapz(fp_p3, tp_p3); AUC_p31 = trapz(fp_p31, tp_p31); 
AUC_p6 = trapz(fp_p6, tp_p6); AUC_p61 = trapz(fp_p61, tp_p61); 
 
%% Plotting figures
% Figure 11
%Distribution of null values
 
title1= sprintf("Dist of R values using Pearsons corr w robust for %s",datalabels{type});
figure('Name',title1,'NumberTitle','off')
histogram(Rdist_1a,40,'Normalization','pdf'); hold on; 
plot(-1.1:0.01:1.1, normpdf(-1.1:0.01:1.1,0,1/sqrt(dfe1)), 'r--','LineWidth', 3);
ax =gca;
ax.FontSize = 24;
ax.FontWeight = 'bold';
ax.XLim = [-1.1 1.1];
ax.XLabel.String = 'Correlation (r)';
saveas(gcf,title1,'tiffn'); 
saveas(gcf,title1,'mfig');

title11= sprintf("Dist of R values using Pearsons corr for %s",datalabels{type});
figure('Name',title11,'NumberTitle','off')
histogram(Rdist_11a,40,'Normalization','pdf' ); hold on; 
plot(-1.1:0.01:1.1, normpdf(-1.1:0.01:1.1,0,1/sqrt(dfe11)), 'r--','LineWidth', 3);
ax =gca;
ax.FontSize = 24;
ax.FontWeight = 'bold';
ax.XLim = [-1.1 1.1];
ax.XLabel.String = 'Correlation (r)';
saveas(gcf,title11,'tiffn'); 
saveas(gcf,title11,'mfig');

title2= sprintf("Dist of R values using AR corr w robust for %s",datalabels{type});
figure('Name',title2,'NumberTitle','off')
histogram(Rdist_2a,40,'Normalization','pdf'); hold on;
plot(-1.1:0.01:1.1, normpdf(-1.1:0.01:1.1,0,1/sqrt(dfe2)), 'r--','LineWidth', 3);
ax =gca;
ax.FontSize = 24;
ax.FontWeight = 'bold';
ax.XLim = [-1.1 1.1];
ax.XLabel.String = 'Correlation (r)';
saveas(gcf,title2,'tiffn'); 
saveas(gcf,title2,'mfig');


title21= sprintf("Dist of R values using AR corr for %s",datalabels{type});
figure('Name',title21,'NumberTitle','off')
histogram(Rdist_21a,40,'Normalization','pdf'); hold on;
plot(-1.1:0.01:1.1, normpdf(-1.1:0.01:1.1,0,1/sqrt(dfe21)), 'r--','LineWidth', 3);
ax =gca;
ax.FontSize = 24;
ax.FontWeight = 'bold';
ax.XLim = [-1.1 1.1];
ax.XLabel.String = 'Correlation (r)';
saveas(gcf,title21,'tiffn'); 
saveas(gcf,title21,'mfig');

title3= sprintf("Dist of R values using Par corr (Short) w robust for %s",datalabels{type});
figure('Name',title3,'NumberTitle','off');
histogram(Rdist_3a,40,'Normalization','pdf'); hold on; 
plot(-1.1:0.01:1.1, normpdf(-1.1:0.01:1.1,0,1/sqrt(mean(dfe3,'all'))), 'r--','LineWidth', 3);
ax =gca;
ax.FontSize = 24;
ax.FontWeight = 'bold';
ax.XLim = [-1.1 1.1];
ax.XLabel.String = 'Correlation (r)';
saveas(gcf,title3,'tiffn'); 
saveas(gcf,title3,'mfig');


title31= sprintf("Dist of R values using Par corr (Short) for %s",datalabels{type});
figure('Name',title31,'NumberTitle','off');
histogram(Rdist_31a,40,'Normalization','pdf'); hold on; 
plot(-1.1:0.01:1.1, normpdf(-1.1:0.01:1.1,0,1/sqrt(mean(dfe31,'all'))), 'r--','LineWidth', 3);
ax =gca;
ax.FontSize = 24;
ax.FontWeight = 'bold';
ax.XLim = [-1.1 1.1];
ax.XLabel.String = 'Correlation (r)';
saveas(gcf,title31,'tiffn'); 
saveas(gcf,title31,'mfig');


% Type-1 Error control plot (Actual p vs bootstraped-p)
titlea= sprintf("FP plot for %s",datalabels{type});
figure('Name',titlea,'NumberTitle','off','Position',[100 100 800 600]);
plot(quantile(pdist_1a(indxR1),int), pvals_int, 'r','LineWidth',2); hold on; 
plot(quantile(pdist_11a(indxR11),int), pvals_int, 'r--','LineWidth',2); hold on; 
plot(quantile(pdist_2a(indxR2),int), pvals_int, 'g','LineWidth',2); hold on; 
plot(quantile(pdist_21a(indxR21),int), pvals_int, 'g--','LineWidth',2); hold on; 
plot(quantile(pdist_3a(indxR3),int), pvals_int, 'b','LineWidth',2); hold on;
plot(quantile(pdist_31a(indxR31),int), pvals_int, 'b--','LineWidth',2);
plot(quantile(pdist_6a(indxG6),int), pvals_int, 'm','LineWidth',2); hold on;
plot(quantile(pdist_61a(indxG61),int), pvals_int, 'm--','LineWidth',2);

xlabel("p-hat"); ylabel(" False positive rate"); pbaspect([1 1 1]);
l = legend(["Pearson's corr w/ robust","Pearson's corr", "AR corr w/ robust", ...
    "AR corr", "AR Par-corr w/ robust", "AR Par-corr", ...
    "Modified MVGC w/ robust", "Modified MVGC"], 'Location','eastoutside');
l.FontSize = 15;
% Plotting a x=y ref line to compare the false positive(FP) rate vs expected FP rate 
hline = refline([1,0]);
hline.LineStyle = ':';
hline.Color = 'k';
hline.LineWidth = 2;
hline.HandleVisibility = 'off';
ax =gca;
ax.FontSize = 18;
ax.FontWeight = 'bold';
saveas(gcf,titlea,'tiffn'); 
saveas(gcf,titlea,'mfig');

% Plotting the ROC curves to assess performance
titleb= sprintf("ROC curves for %s",datalabels{type});
figure('Name',titleb,'NumberTitle','off','Position',[100 100 800 600]);
plot(fp_p1,tp_p1,'r','LineWidth',2); hold on; 
plot(fp_p11,tp_p11,'r--','LineWidth',2); hold on; 
plot(fp_p2,tp_p2, 'g','LineWidth',2); hold on; 
plot(fp_p21,tp_p21, 'g--','LineWidth',2); hold on; 
plot(fp_p3,tp_p3, 'b','LineWidth',2); hold on;
plot(fp_p31,tp_p31, 'b--','LineWidth',2);
plot(fp_p6,tp_p6, 'm','LineWidth',2); hold on;
plot(fp_p61,tp_p61, 'm--','LineWidth',2);

xlabel('False positive rate'); ylabel('True positive rate'); pbaspect([1 1 1]);
l = legend(["Pearson's corr w/ robust","Pearson's corr", "AR corr w/ robust", ...
    "AR corr", "AR Par-corr w/ robust", "AR Par-corr", ...
    "Modified MVGC w/ robust", "Modified MVGC"], 'Location','eastoutside');

l.FontSize = 15;
% Plotting a x=y ref line to show performance obtained by random chance
hline = refline([1,0]);
hline.LineStyle = ':';
hline.Color = 'k';
hline.LineWidth = 2;
hline.HandleVisibility = 'off';
ax =gca;
ax.FontSize = 18;
ax.FontWeight = 'bold';
saveas(gcf,titleb,'tiffn'); 
saveas(gcf,titleb,'mfig');
cd(Curr_Folder)
end

%close all;

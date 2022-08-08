% This code creates the figure 9 in the manuscript

% Creating folder to save results

Curr_Folder = pwd; % Get path of working directory
% This code requires NIRS Brain AnalyzIR toolbox to run. The toolbox can be
% downloaded from  https://github.com/huppertt/nirs-toolbox. 
% Please download the toolbox and add to path

addpath(genpath('nirs-toolbox')); % Add toolbox to path

addpath(genpath('functions')); % Add to path folder containing important functions

%Creating folder to save results
fig_dir  = [Curr_Folder filesep 'Parcorr'];
mkdir(fig_dir)


% Preprocessing pipeline to preprocess resting-state fNIRS data

job = nirs.modules.OpticalDensity;
job = nirs.modules.BeerLambertLaw(job);
job = nirs.modules.KeepTypes(job);
job.types = 'hbo'; % Keep only HBO as similar results should apply to HbR
job = nirs.modules.LabelShortSeperation(job);



% Generate the simulated ground truth. This ground truth will be used to
% evaluate the performance of partial correlation pipelines
[raw , truth] = simData_connectivity_shortsep([],[1],0); % truth contains 
% the generated covariance matrix and raw the corresponding data
hboRest = job.run(raw); 

% Setting simulation Parameters
nIters =100; % No of iterations, no of datasets to be generated
Fs = raw.Fs;  %Sampling rate
Pmax = round(10*Fs);% Max model order for autoregressive (AR) models
criterion = 'BIC'; %Informtation criteria to compare model fits to find optimal model orders in AR and MVAR models 

%Extract info from data
nchan = sum(~raw.probe.link.ShortSeperation)/2; %No of long channels
nshort= sum(raw.probe.link.ShortSeperation)/2;  %No of short channels
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
datalabels = {'random data','w temp corr noise',....
    'w temp & spatial corr noise', 'w temp & spatial corr noise & motion'};

 for type=3:4 % Just running for temporal autocorrelated+ spatial correlated noise  
     % & with motion should be sufficient all AR partial correlation methods
     % correct for temporally correlated noise and effectivness is
     % evaluated in removing the global spatially correlated noise.
     
  mkdir([Curr_Folder filesep 'Parcorr' filesep datalabels{type}])
  cd([Curr_Folder filesep 'Parcorr' filesep datalabels{type}])

   
%% Running ths simulation
    
% Initializing variables to store values obtained from the analysis pipeline with zeros
Rdist_3a = zeros(pathsa,nIters); Rdist_4a = zeros(pathsa,nIters);
Zdist_3a = zeros(pathsa,nIters); Zdist_4a = zeros(pathsa,nIters);
pdist_3a = zeros(pathsa,nIters); pdist_4a = zeros(pathsa,nIters); 
Rdist_31a = zeros(pathsa,nIters); Rdist_41a = zeros(pathsa,nIters);
Zdist_31a = zeros(pathsa,nIters); Zdist_41a = zeros(pathsa,nIters);
pdist_31a = zeros(pathsa,nIters); pdist_41a = zeros(pathsa,nIters);
Rdist_5a = zeros(pathsa,nIters);  Rdist_51a = zeros(pathsa,nIters);
Zdist_5a = zeros(pathsa,nIters); Zdist_51a = zeros(pathsa,nIters);
pdist_5a = zeros(pathsa,nIters); pdist_51a = zeros(pathsa,nIters);
Rdist_32a = zeros(pathsa,nIters); Rdist_33a = zeros(pathsa,nIters);
Zdist_32a = zeros(pathsa,nIters); Zdist_33a = zeros(pathsa,nIters);
pdist_32a = zeros(pathsa,nIters); pdist_33a = zeros(pathsa,nIters); 

Rdist_3b = zeros(pathsb,nIters); Rdist_4b = zeros(pathsb,nIters);
Zdist_3b = zeros(pathsb,nIters); Zdist_4b = zeros(pathsb,nIters);
pdist_3b = zeros(pathsb,nIters); pdist_4b = zeros(pathsb,nIters); 
Rdist_31b = zeros(pathsb,nIters); Rdist_41b = zeros(pathsb,nIters);
Zdist_31b = zeros(pathsb,nIters); Zdist_41b = zeros(pathsb,nIters);
pdist_31b = zeros(pathsb,nIters); pdist_41b = zeros(pathsb,nIters);
Rdist_5b = zeros(pathsb,nIters); Rdist_51b = zeros(pathsb,nIters); 
Zdist_5b = zeros(pathsb,nIters); Zdist_51b = zeros(pathsb,nIters); 
pdist_5b = zeros(pathsb,nIters); pdist_51b = zeros(pathsb,nIters); 
Rdist_32b = zeros(pathsb,nIters); Rdist_33b = zeros(pathsb,nIters);
Zdist_32b = zeros(pathsb,nIters); Zdist_33b = zeros(pathsb,nIters);
pdist_32b = zeros(pathsb,nIters); pdist_33b = zeros(pathsb,nIters);


% Running several iterations
 for iter=1:nIters  

      % Simulate the data based on statisical properites

    switch type-1
        case 1
          [raw1 , ~] = simData_connectivity_shortsep(truth,[1],0);
          % With temprally autocorrelated noise
          data = job.run(raw1);

        case 2
          [raw1 , ~] = simData_connectivity_shortsep(truth,[1],150);
          % With temporall autocorrelated and spatially correlated physiological noise
          data = job.run(raw1);
        
        case 3
          % WIth motion artifacts & spatially correlated physiological noise
          [raw1 , ~] = simData_connectivity_shortsep(truth,[1],150);
          raw1 = nirs.testing.simMotionArtifact(raw1); % adding motion artifacts
          data = job.run(raw1);

        otherwise
         % Random noise
          [raw1 , ~] = simData_connectivity_shortsep(truth,[1],0);
          % Without  spatially or tremporal correlated noise
          data = job.run(raw1);
          data.data = rand(ntp,nchan+nshort);
    end

   % Partial correlation with short channels with robust regression
    [R3, p3, dfe3] = partial_corr(data,'10xFs','short',true, false);
    Z3= atanh(R3); % Fisher's r to z transform
    Rdist_3a(:, iter) = R3(zeroind);
    Zdist_3a(:, iter) = Z3(zeroind);
    pdist_3a(:, iter) = p3(zeroind);
    Rdist_3b(:, iter) = cat(1,R3(nnzind),R3(zeroind));
    Zdist_3b(:, iter) = cat(1,Z3(nnzind),Z3(zeroind));
    pdist_3b(:, iter) = cat(1,p3(nnzind),p3(zeroind));
    
    % Partial correlation with short channels, without robust regression
    [R31, p31, dfe31] = partial_corr(data,'10xFs','short',false, false);
    Z31= atanh(R31); % Fisher's r to z transform
    Rdist_31a(:, iter) = R31(zeroind);
    Zdist_31a(:, iter) = Z31(zeroind);
    pdist_31a(:, iter) = p31(zeroind);
    Rdist_31b(:, iter) = cat(1,R31(nnzind),R31(zeroind));
    Zdist_31b(:, iter) = cat(1,Z31(nnzind),Z31(zeroind));
    pdist_31b(:, iter) = cat(1,p31(nnzind),p31(zeroind));
    
    % Partial correlation with PCA on short channels, and robust regression
    [R32, p32, dfe32] = partial_corr(data,'10xFs','short',true, true);
    Z32 = atanh(R32); % Fisher's r to z transform
    Rdist_32a(:, iter) = R32(zeroind);
    Zdist_32a(:, iter) = Z32(zeroind);
    pdist_32a(:, iter) = p32(zeroind);
    Rdist_32b(:, iter) = cat(1,R32(nnzind),R32(zeroind));
    Zdist_32b(:, iter) = cat(1,Z32(nnzind),Z32(zeroind));
    pdist_32b(:, iter) = cat(1,p32(nnzind),p32(zeroind));
    
    % Partial correlation with PCA on short channels, without robust regression
    [R33, p33, dfe33] = partial_corr(data,'10xFs','short',false, true);
    Z33 = atanh(R33); % Fisher's r to z transform
    Rdist_33a(:, iter) = R33(zeroind);
    Zdist_33a(:, iter) = Z33(zeroind);
    pdist_33a(:, iter) = p33(zeroind);
    Rdist_33b(:, iter) = cat(1,R33(nnzind),R33(zeroind));
    Zdist_33b(:, iter) = cat(1,Z33(nnzind),Z33(zeroind));
    pdist_33b(:, iter) = cat(1,p33(nnzind),p33(zeroind));
    
    % Partial correlation with PCA on long channels and robust regression
    [R4, p4, dfe4] = partial_corr(data,'10xFs','long',true, true);
    Z4= atanh(R4);  % Fisher's r to z transform
    Rdist_4a(:, iter) = R4(zeroind);
    Zdist_4a(:, iter) = Z4(zeroind);
    pdist_4a(:, iter) = p4(zeroind);
    Rdist_4b(:, iter) = cat(1,R4(nnzind),R4(zeroind));
    Zdist_4b(:, iter) = cat(1,Z4(nnzind),Z4(zeroind));
    pdist_4b(:, iter) = cat(1,p4(nnzind),p4(zeroind));
    
    % Partial correlation with PCA on long channels without robust regression
    [R41, p41, dfe41] = partial_corr(data,'10xFs','long',false, true);
    Z41= atanh(R41);  % Fisher's r to z transform
    Rdist_41a(:, iter) = R41(zeroind);
    Zdist_41a(:, iter) = Z41(zeroind);
    pdist_41a(:, iter) = p41(zeroind);
    Rdist_41b(:, iter) = cat(1,R41(nnzind),R41(zeroind));
    Zdist_41b(:, iter) = cat(1,Z41(nnzind),Z41(zeroind));
    pdist_41b(:, iter) = cat(1,p41(nnzind),p41(zeroind));
    

    % Partial correlation with PCA on both long and short channels and robust regression
    [R5, p5, dfe5] = partial_corr(data,'10xFs','all',true, true);
    Z5= atanh(R5);  % Fisher's r to z transform
    Rdist_5a(:, iter) = R5(zeroind);
    Zdist_5a(:, iter) = Z5(zeroind);
    pdist_5a(:, iter) = p5(zeroind);
    Rdist_5b(:, iter) = cat(1,R5(nnzind),R5(zeroind));
    Zdist_5b(:, iter) = cat(1,Z5(nnzind),Z5(zeroind));
    pdist_5b(:, iter) = cat(1,p5(nnzind),p5(zeroind));
    
    % Partial correlation with PCA on long and short channels without robust regression
    [R51, p51, dfe51] = partial_corr(data,'10xFs','all',false, true);
    Z51= atanh(R51);  % Fisher's r to z transform
    Rdist_51a(:, iter) = R51(zeroind);
    Zdist_51a(:, iter) = Z51(zeroind);
    pdist_51a(:, iter) = p51(zeroind);
    Rdist_51b(:, iter) = cat(1,R51(nnzind),R51(zeroind));
    Zdist_51b(:, iter) = cat(1,Z51(nnzind),Z51(zeroind));
    pdist_51b(:, iter) = cat(1,p51(nnzind),p51(zeroind));
    
     fprintf("Iteration completed: %d of %d\n", iter, nIters);
    
 end
 
 %% Calculating the necessary metrics
 
 %  Converting the r values obtrained from anlaysis pipelines into column vector for later plotting of
 % bootstrap -p vs actual p
 Rdist_3a = Rdist_3a(:); Zdist_3a = Zdist_3a(:); pdist_3a = pdist_3a(:); 
 Rdist_31a = Rdist_31a(:); Zdist_31a = Zdist_31a(:); pdist_31a = pdist_31a(:);
 Rdist_32a = Rdist_32a(:); Zdist_32a = Zdist_32a(:); pdist_32a = pdist_32a(:); 
 Rdist_33a = Rdist_33a(:); Zdist_33a = Zdist_33a(:); pdist_33a = pdist_33a(:);
 Rdist_4a = Rdist_4a(:); Zdist_4a = Zdist_4a(:); pdist_4a = pdist_4a(:); 
 Rdist_41a = Rdist_41a(:); Zdist_41a = Zdist_41a(:); pdist_41a = pdist_41a(:);
 Rdist_5a = Rdist_5a(:); Zdist_5a = Zdist_5a(:); pdist_5a = pdist_5a(:); 
 Rdist_51a = Rdist_51a(:); Zdist_51a = Zdist_51a(:); pdist_51a = pdist_51a(:);
     
 % Sorting the R values obtrained in ascending order 
 [R_sorted_3, indxR3]= sort(Rdist_3a, 'ascend'); [R_sorted_31, indxR31]= sort(Rdist_31a, 'ascend'); 
 [R_sorted_32, indxR32]= sort(Rdist_32a, 'ascend'); [R_sorted_33, indxR33]= sort(Rdist_33a, 'ascend'); 
 [R_sorted_4, indxR4]= sort(Rdist_4a, 'ascend'); [R_sorted_41, indxR41]= sort(Rdist_41a, 'ascend'); 
 [R_sorted_5, indxR5]= sort(Rdist_5a, 'ascend'); [R_sorted_51, indxR51]= sort(Rdist_51a, 'ascend'); 
 
 
 % Calculations for generating ROC curves and calculating AUC
 
 indx1 =  randperm(tp);
 indx2 =  randperm(tn);
 
% Since we have less true positives, we randomly sample a subset of true
% negatives to plot ROC curves which require equal number of true positives
% & negatives

 pdist_3b = cat(1,pdist_3b(indx1(1:nvals),:),pdist_3b(indx2(1:nvals),:));
 pdist_31b = cat(1,pdist_31b(indx1(1:nvals),:),pdist_31b(indx2(1:nvals),:));
 pdist_32b = cat(1,pdist_32b(indx1(1:nvals),:),pdist_32b(indx2(1:nvals),:));
 pdist_33b = cat(1,pdist_33b(indx1(1:nvals),:),pdist_33b(indx2(1:nvals),:));
 pdist_4b = cat(1,pdist_4b(indx1(1:nvals),:),pdist_4b(indx2(1:nvals),:));
 pdist_41b = cat(1,pdist_41b(indx1(1:nvals),:),pdist_41b(indx2(1:nvals),:));
 pdist_5b = cat(1,pdist_5b(indx1(1:nvals),:),pdist_5b(indx2(1:nvals),:));
 pdist_51b = cat(1,pdist_51b(indx1(1:nvals),:),pdist_51b(indx2(1:nvals),:));
 
 pdist_3b = pdist_3b(:); pdist_31b = pdist_31b(:); 
 pdist_32b = pdist_32b(:); pdist_33b = pdist_33b(:); 
 pdist_4b = pdist_4b(:); pdist_41b = pdist_41b(:);   
 pdist_5b = pdist_5b(:); pdist_51b = pdist_51b(:); 
 
 % ROC curve analysis
[tp_p3,fp_p3,th_p1] = nirs.testing.roc(truth_all_p, pdist_3b); 
[tp_p31,fp_p31,th_p31] = nirs.testing.roc(truth_all_p, pdist_31b); 
[tp_p32,fp_p32,th_p32] = nirs.testing.roc(truth_all_p, pdist_32b); 
[tp_p33,fp_p33,th_p33] = nirs.testing.roc(truth_all_p, pdist_33b);
[tp_p4,fp_p4,th_p4] = nirs.testing.roc(truth_all_p, pdist_4b);
[tp_p41,fp_p41,th_p41] = nirs.testing.roc(truth_all_p, pdist_41b);
[tp_p5,fp_p5,th_p5] = nirs.testing.roc(truth_all_p, pdist_5b);
[tp_p51,fp_p51,th_p51] = nirs.testing.roc(truth_all_p, pdist_51b); 

% AUC Metrics 
AUC_p3 = trapz(fp_p3, tp_p3); AUC_p31 = trapz(fp_p31, tp_p31); 
AUC_p32 = trapz(fp_p32, tp_p32); AUC_p33 = trapz(fp_p33, tp_p33); 
AUC_p4 = trapz(fp_p4, tp_p4); AUC_p41 = trapz(fp_p41, tp_p41); 
AUC_p5 = trapz(fp_p5, tp_p5); AUC_p51 = trapz(fp_p51, tp_p51); 

 
%% Plotting Figures

% Type-1 Error control plot (Actual p vs bootstraped-p)
titlea= sprintf("FP plot for %s",datalabels{type});
figure('Name',titlea,'NumberTitle','off','Position',[100 100 800 600]);
plot(quantile(pdist_3a(indxR3),int), pvals_int, 'r','LineWidth',2); hold on; 
plot(quantile(pdist_31a(indxR31),int), pvals_int, 'r--','LineWidth',2); hold on; 
plot(quantile(pdist_32a(indxR32),int), pvals_int, 'm','LineWidth',2); hold on; 
plot(quantile(pdist_33a(indxR33),int), pvals_int, 'm--','LineWidth',2); hold on; 
plot(quantile(pdist_4a(indxR4),int), pvals_int, 'g','LineWidth',2); hold on; 
plot(quantile(pdist_41a(indxR41),int), pvals_int, 'g--','LineWidth',2); hold on; 
plot(quantile(pdist_5a(indxR5),int), pvals_int, 'b','LineWidth',2); hold on;
plot(quantile(pdist_51a(indxR51),int), pvals_int, 'b--','LineWidth',2);

xlabel("p-hat"); ylabel(" False positive rate"); pbaspect([1 1 1]);
l= legend(["Short AR Par-corr w/ robust" ,"Short AR Par-corr",...
    "Short AR Par-corr w/ PCA & w/ robust","Short AR Par-corr w/ PCA",...
    "Long AR Par-corr w/ PCA & w/ robust","Long AR Par-corr w/ PCA", ...
    "Short+Long AR Par-corr w/ PCA & w/ robust","Short+Long AR Par-corr w/ PCA"], 'Location','eastoutside');
l.FontSize = 12;
% Plotting a x=y ref line to compare the false positive(FP) rate vs expected FP rate 
hline = refline([1,0]);
hline.LineStyle = ':';
hline.LineWidth = 2;
hline.Color = 'k';
hline.HandleVisibility = 'off';
ax =gca;
ax.FontSize = 15;
ax.FontWeight = 'bold';
saveas(gcf,titlea,'tiffn'); 
saveas(gcf,titlea,'mfig');

% Plotting the ROC curves to assess performance
titleb= sprintf("ROC curves for %s",datalabels{type});
figure('Name',titleb,'NumberTitle','off','Position',[100 100 800 600]);
plot(fp_p3,tp_p3,'r','LineWidth',2); hold on; 
plot(fp_p31,tp_p31,'r--','LineWidth',2); hold on; 
plot(fp_p32,tp_p32,'m','LineWidth',2); hold on; 
plot(fp_p33,tp_p33,'m--','LineWidth',2); hold on; 
plot(fp_p4,tp_p4, 'g','LineWidth',2); hold on; 
plot(fp_p41,tp_p41, 'g--','LineWidth',2); hold on; 
plot(fp_p5,tp_p5, 'b','LineWidth',2); hold on;
plot(fp_p51,tp_p51, 'b--','LineWidth',2);

xlabel('False positive rate'); ylabel('True positive rate'); pbaspect([1 1 1]);
l= legend(["Short AR Par-corr w/ robust" ,"Short AR Par-corr",...
    "Short AR Par-corr w/ PCA & w/ robust","Short AR Par-corr w/ PCA",...
    "Long AR Par-corr w/ PCA & w/ robust","Long AR Par-corr w/ PCA", ...
    "Short+Long AR Par-corr w/ PCA & w/ robust","Short+Long AR Par-corr w/ PCA"], 'Location','eastoutside');
l.FontSize = 12;
% Plotting a x=y ref line to show performance obtained by random chance
hline = refline([1,0]);
hline.LineStyle = ':';
hline.Color = 'k';
hline.LineWidth = 2;
hline.HandleVisibility = 'off';
ax =gca;
ax.FontSize = 15;
ax.FontWeight = 'bold';
saveas(gcf,titleb,'tiffn'); 
saveas(gcf,titleb,'mfig');
cd(Curr_Folder)
end

%close all;


% This code creates the figure 10 in the manuscript

% Creating folder to save results

Curr_Folder = pwd; % Get path of working directory
% This code requires NIRS Brain AnalyzIR toolbox to run. The toolbox can be
% downloaded from  https://github.com/huppertt/nirs-toolbox. 
% Please download the toolbox and add to path

addpath(genpath('nirs-toolbox')); % Add toolbox to path

addpath(genpath('functions')); % Add to path folder containing important functions

%Creating folder to save results
fig_dir  = [Curr_Folder filesep 'Grangers'];
mkdir(fig_dir)


% Preprocessing pipeline to preprocess resting-state fNIRS data

job = nirs.modules.OpticalDensity;
job = nirs.modules.BeerLambertLaw(job);
job = nirs.modules.KeepTypes(job);
job.types = 'hbo';  % Keep only HBO as similar results should apply to HbR
job = nirs.modules.LabelShortSeperation(job);

% Generate the simulated ground truth. This ground truth will be used to
% evaluate the performance of partial correlation pipelines
[raw , truth] = simData_connectivity_shortsep([],[0,1],0); % Data simulated with lag relationship
% truth contains the generated covariance matrix and raw the corresponding data
hboRest = job.run(raw);


% Setting simulation Parameters
nIters =100; % No of iterations, no of datasets to be generated
Fs = raw.Fs;  %Sampling rate
Pmax = round(10*Fs);% Max model order for autoregressive (AR) models
criterion = 'BIC'; %Informtation criteria to compare model fits to find optimal model orders in AR and MVAR models 

%Extract info from data
nchan = sum(~raw.probe.link.ShortSeperation)/2;  %No of long channels
nshort= sum(raw.probe.link.ShortSeperation)/2; %No of short channels
ntp = size(raw.data, 1);
truth_lag = squeeze(truth(:,:,2));
zeroind = triu(true(size(truth_lag)),1)&(~truth_lag);%Indices of elements with true negatives in the data covariance  matrix
dataind = triu(true(nchan),1);% getting the upper right hand corner of the data covariance matrix.

pathsa = nnz(zeroind);%True negative paths
pathsb = 0.5*nchan*(nchan-1);% All paths


% Truth distribution
% This will later be used to compare the performance of simulation results

t = sqrt(ntp-Pmax).*((truth_lag)./(sqrt(1-truth_lag.^2)));
pvals = 2*nirs.math.tpvalue(-abs(t),ntp-Pmax);
truth_binp  = pvals<0.05;%Binarizing the "ground truth".
% ones in the truth_binp denote a connection between channels and 0 denotes no connection
nnzind = triu(true(size(truth_binp)),1)&(truth_binp);%Indices of elements with true positives
%in the data covariance  matrix
tp = nnz(nnzind); % No of true positives
tn = nnz(zeroind); % No of true negatives
nvals = min(tp,tn); %Often are data has less true positives than true negatives
truth_all_p = repmat(vertcat(ones(nvals,1),zeros(nvals,1)),nIters,1);% Will later be used to generate ROC
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
    'w temp & spatial corr noise', 'w temp & spatial corr noise & motion',};

for type=3:4 % Just running for temporal autocorrelated+ spatial correlated noise
     % & with motion should be sufficient all MVGC methods
     % correct for temporally correlated noise and effectivness is
     % evaluated in removing the global spatially correlated noise.
     
  mkdir([Curr_Folder filesep 'Grangers' filesep datalabels{type}])
  cd([Curr_Folder filesep 'Grangers' filesep datalabels{type}])
   
%% Running ths simulation
    
% Initializing variables to store values obtained from the analysis pipeline with zeros
Fdist_21a = zeros(pathsa,nIters); Fdist_2a = zeros(pathsa,nIters);
pdist_21a = zeros(pathsa,nIters); pdist_2a = zeros(pathsa,nIters); 
Gdist_21a = zeros(pathsa,nIters); Gdist_2a = zeros(pathsa,nIters);

Fdist_3a = zeros(pathsa,nIters); Fdist_4a = zeros(pathsa,nIters);
pdist_3a = zeros(pathsa,nIters); pdist_4a = zeros(pathsa,nIters);
Gdist_3a = zeros(pathsa,nIters); Gdist_4a = zeros(pathsa,nIters);

Fdist_31a = zeros(pathsa,nIters); Fdist_41a = zeros(pathsa,nIters);
pdist_31a = zeros(pathsa,nIters); pdist_41a = zeros(pathsa,nIters);
Gdist_31a = zeros(pathsa,nIters); Gdist_41a = zeros(pathsa,nIters);

Fdist_21b = zeros(pathsb,nIters); Fdist_2b = zeros(pathsb,nIters);
pdist_21b = zeros(pathsb,nIters); pdist_2b = zeros(pathsb,nIters); 
Gdist_21b = zeros(pathsb,nIters); Gdist_2b = zeros(pathsb,nIters);

Fdist_3b = zeros(pathsb,nIters); Fdist_4b = zeros(pathsb,nIters);
pdist_3b = zeros(pathsb,nIters); pdist_4b = zeros(pathsb,nIters);
Gdist_3b = zeros(pathsb,nIters); Gdist_4b = zeros(pathsb,nIters);

Fdist_31b = zeros(pathsb,nIters); Fdist_41b = zeros(pathsb,nIters);
pdist_31b = zeros(pathsb,nIters); pdist_41b = zeros(pathsb,nIters);
Gdist_31b = zeros(pathsb,nIters); Gdist_41b = zeros(pathsb,nIters);


% Running several iterations

 for iter=1:nIters  
   
   % Simulate the data
   
% Simulate the data based on statisical properites
    switch type-1
        case 1
          [raw1 , ~] = simData_connectivity_shortsep(truth,[0,1],0);
          % With temporally autocorrelated noise
          data = job.run(raw1);

        case 2
          [raw1 , ~] = simData_connectivity_shortsep(truth,[0,1],150);
          % With temporally autocorrelated and spatially correlated physiological noise
          data = job.run(raw1);

        case 3
          [raw1 , ~] = simData_connectivity_shortsep(truth,[0,1],150);
          raw1 = nirs.testing.simMotionArtifact(raw1);  
          % WIth motion artifacts & spatially correlated physiological noise
          data = job.run(raw1);

        otherwise
          % Random noise
          [raw1 , ~] = simData_connectivity_shortsep(truth,[0,1],0);
          % Without  spatially or tremporal correlated noise   
          data = job.run(raw1);
          data.data = rand(ntp,nchan+nshort);
    end


 % Run multivariate gc code controlling for shortchannels (After PCA) & with robust regression
    [G2, F2, df12, df22, p2] = grangercausality(data,Pmax,'multivariate','short',false,true,true,criterion);
    Gdist_2a(:, iter) = G2(zeroind);
    Fdist_2a(:, iter) = F2(zeroind);
    pdist_2a(:, iter) = p2(zeroind);
    Gdist_2b(:, iter) = cat(1,G2(nnzind),G2(zeroind));
    Fdist_2b(:, iter) = cat(1,F2(nnzind),F2(zeroind));
    pdist_2b(:, iter) = cat(1,p2(nnzind),p2(zeroind));
    
 % Run multivariate gc code controlling for shortchannels (After PCA) & without robust regression
    [G21, F21, df121, df221, p21] = grangercausality(data,Pmax,'multivariate','short',false,false,true,criterion);
    Gdist_21a(:, iter) = G21(zeroind);
    Fdist_21a(:, iter) = F21(zeroind);
    pdist_21a(:, iter) = p21(zeroind);
    Gdist_21b(:, iter) = cat(1,G21(nnzind),G21(zeroind));
    Fdist_21b(:, iter) = cat(1,F21(nnzind),F21(zeroind));
    pdist_21b(:, iter) = cat(1,p21(nnzind),p21(zeroind));

 % Run multivariate gc code controlling for long channels (After PCA) & with robust regression
    [G3, F3, df13, df23, p3] = grangercausality(data,Pmax,'multivariate','long',false,true,true,criterion);
    Gdist_3a(:, iter) = G3(zeroind);
    Fdist_3a(:, iter) = F3(zeroind);
    pdist_3a(:, iter) = p3(zeroind);
    Gdist_3b(:, iter) = cat(1,G3(nnzind),G3(zeroind));
    Fdist_3b(:, iter) = cat(1,F3(nnzind),F3(zeroind));
    pdist_3b(:, iter) = cat(1,p3(nnzind),p3(zeroind));
    
  % Run multivariate gc code  controlling for long channels (After PCA) & without robust regression
    [G31, F31, df131, df231, p31] = grangercausality(data,Pmax,'multivariate','long',false,false,true,criterion);
    Gdist_31a(:, iter) = G31(zeroind);
    Fdist_31a(:, iter) = F31(zeroind);
    pdist_31a(:, iter) = p31(zeroind);
    Gdist_31b(:, iter) = cat(1,G31(nnzind),G31(zeroind));
    Fdist_31b(:, iter) = cat(1,F31(nnzind),F31(zeroind));
    pdist_31b(:, iter) = cat(1,p31(nnzind),p31(zeroind));
    
 % Run multivariate gc code controlling for short+long channels (After PCA) & with robust regression
    [G4, F4, df14, df24, p4] = grangercausality(data,Pmax,'multivariate','all',false,true,true,criterion);
    Gdist_4a(:, iter) = G4(zeroind);
    Fdist_4a(:, iter) = F4(zeroind);
    pdist_4a(:, iter) = p4(zeroind);
    Gdist_4b(:, iter) = cat(1,G4(nnzind),G4(zeroind));
    Fdist_4b(:, iter) = cat(1,F4(nnzind),F4(zeroind));
    pdist_4b(:, iter) = cat(1,p4(nnzind),p4(zeroind));
    
  % Run multivariate gc code controlling for short+long channels (After PCA) & without robust regression
    [G41, F41, df141, df241, p41] = grangercausality(data,Pmax,'multivariate','all',false,false,true,criterion);
    Gdist_41a(:, iter) = G41(zeroind);
    Fdist_41a(:, iter) = F41(zeroind);
    pdist_41a(:, iter) = p41(zeroind);
    Gdist_41b(:, iter) = cat(1,G41(nnzind),G41(zeroind));
    Fdist_41b(:, iter) = cat(1,F41(nnzind),F41(zeroind));
    pdist_41b(:, iter) = cat(1,p41(nnzind),p41(zeroind));
    
fprintf("Iteration completed: %d of %d\n", iter, nIters);
 end
 
 %% Calculating the necessary metrics
 
 %  Converting the r values obtrained from anlaysis pipelines into column vector for later plotting of
 % bootstrap -p vs actual p

 Gdist_2a = Gdist_2a(:); pdist_2a = pdist_2a(:);
 Gdist_3a = Gdist_3a(:); pdist_3a = pdist_3a(:); 
 Gdist_4a = Gdist_4a(:); pdist_4a = pdist_4a(:);
 Gdist_21a = Gdist_21a(:); pdist_21a = pdist_21a(:); 
 Gdist_31a = Gdist_31a(:); pdist_31a = pdist_31a(:); 
 Gdist_41a = Gdist_41a(:); pdist_41a = pdist_41a(:);
 
 % Sorting the R values obtrained in ascending order 
 [G_sorted_2, indxG2]= sort(Gdist_2a, 'ascend');
 [G_sorted_21, indxG21]= sort(Gdist_21a, 'ascend'); 
 [G_sorted_3, indxG3]= sort(Gdist_3a, 'ascend'); 
 [G_sorted_4, indxG4]= sort(Gdist_4a, 'ascend'); 
 [G_sorted_31, indxG31]= sort(Gdist_31a, 'ascend'); 
 [G_sorted_41, indxG41]= sort(Gdist_41a, 'ascend'); 
 
 
 % Calculations for generating ROC curves and calculating AUC
 indx1 =  randperm(tp);
 indx2 =  randperm(tn);
 
 % Since we have less true positives, we randomly sample a subset of true
 % negatives to plot ROC curves which require equal number of true positives
 % & negatives
 
 pdist_21b = cat(1,pdist_21b(indx1(1:nvals),:),pdist_21b(indx2(1:nvals),:));
 pdist_2b = cat(1,pdist_2b(indx1(1:nvals),:),pdist_2b(indx2(1:nvals),:));
 pdist_3b = cat(1,pdist_3b(indx1(1:nvals),:),pdist_3b(indx2(1:nvals),:));
 pdist_4b = cat(1,pdist_4b(indx1(1:nvals),:),pdist_4b(indx2(1:nvals),:));
 pdist_31b = cat(1,pdist_31b(indx1(1:nvals),:),pdist_31b(indx2(1:nvals),:));
 pdist_41b = cat(1,pdist_41b(indx1(1:nvals),:),pdist_41b(indx2(1:nvals),:));
 
 
 pdist_2b = pdist_2b(:); pdist_21b = pdist_21b(:); 
 pdist_3b = pdist_3b(:); pdist_31b = pdist_31b(:);
 pdist_4b = pdist_4b(:); pdist_41b = pdist_41b(:); 
 
 % ROC curve analysis  
[tp_p21,fp_p21,th_p21] = nirs.testing.roc(truth_all_p, pdist_21b); 
[tp_p2,fp_p2,th_p2] = nirs.testing.roc(truth_all_p, pdist_2b);
[tp_p3,fp_p3,th_p3] = nirs.testing.roc(truth_all_p, pdist_3b);
[tp_p4,fp_p4,th_p4] = nirs.testing.roc(truth_all_p, pdist_4b);
[tp_p31,fp_p31,th_p31] = nirs.testing.roc(truth_all_p, pdist_31b);
[tp_p41,fp_p41,th_p41] = nirs.testing.roc(truth_all_p, pdist_41b);

% AUC Metrics
AUC_p2 = trapz(fp_p2, tp_p2); AUC_p21 = trapz(fp_p21, tp_p21); 
AUC_p3 = trapz(fp_p3, tp_p3); AUC_p31 = trapz(fp_p31, tp_p31);
AUC_p4 = trapz(fp_p4, tp_p4); AUC_p41 = trapz(fp_p41, tp_p41);


 
%% Plotting Figures

% Type-1 Error control plot (Actual p vs bootstraped-p)
titlea= sprintf("FP plot for %s",datalabels{type});
figure('Name',titlea,'NumberTitle','off','Position',[100 100 800 600]);
plot(quantile(pdist_2a(indxG2),int), pvals_int, 'c','LineWidth',2); hold on; 
plot(quantile(pdist_21a(indxG21),int), pvals_int, 'c--','LineWidth',2); hold on; 
plot(quantile(pdist_3a(indxG3),int), pvals_int, 'm','LineWidth',2); hold on; 
plot(quantile(pdist_31a(indxG31),int), pvals_int, 'm--','LineWidth',2); hold on; 
plot(quantile(pdist_4a(indxG4),int), pvals_int, 'Color',[1, 0.5, 0],'LineWidth',2); hold on;
plot(quantile(pdist_41a(indxG41),int), pvals_int, 'Color',[1, 0.5, 0], 'LineStyle','--','LineWidth',2); hold on;

xlabel("p-hat"); ylabel(" False positive rate"); pbaspect([1 1 1]);
l = legend(["Short MVGC w/ robust", "Short MVGC" ,"Long MVGC w/ robust",...
    "Long MVGC", "Short+Long MVCG w/ robust","Short+Long MVCG"], 'Location','eastoutside');
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
saveas(gcf,titlea,'mfig')

% Plotting the ROC curves to assess performance
titleb= sprintf("ROC curves for %s",datalabels{type});
figure('Name',titleb,'NumberTitle','off','Position',[100 100 800 600]);
plot(fp_p2,tp_p2,'c','LineWidth',2); hold on;  
plot(fp_p21,tp_p21,'c--','LineWidth',2); hold on;  
plot(fp_p3,tp_p3, 'm','LineWidth',2); hold on; 
plot(fp_p31,tp_p31, 'm--','LineWidth',2); hold on; 
plot(fp_p4,tp_p4, 'Color',[1, 0.5, 0],'LineWidth',2); hold on;
plot(fp_p41,tp_p41, 'Color',[1, 0.5, 0], 'LineStyle','--','LineWidth',2); hold on;

xlabel('False positive rate'); ylabel('True positive rate'); pbaspect([1 1 1]);
l = legend(["Short MVGC w/ robust", "Short MVGC" ,"Long MVGC w/ robust",...
    "Long MVGC", "Short+Long MVCG w/ robust","Short+Long MVCG"], 'Location','eastoutside');
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
saveas(gcf,titleb,'mfig')
cd(Curr_Folder)
end
%close all;

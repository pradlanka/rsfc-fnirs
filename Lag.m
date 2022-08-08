% This code creates the figures 11 & 12 in the manuscript

Curr_Folder = pwd; % Get path of working directory
% Creating folder to save results

% This code requires NIRS Brain AnalyzIR toolbox to run. The toolbox can be
% downloaded from  https://github.com/huppertt/nirs-toolbox. 


% Please download the toolbox and add to path

addpath(genpath('nirs-toolbox')); % Add toolbox to path

addpath(genpath('Other Codes')); % Add to path folder containing important functions


%Creating folder to save results
fig_dir  = [Curr_Folder filesep 'Lag'];
mkdir(fig_dir)


% Preprocessing pipeline to preprocess resting-state fNIRS data

job = nirs.modules.OpticalDensity;
job = nirs.modules.BeerLambertLaw(job);
job = nirs.modules.KeepTypes(job);
job.types = 'hbo'; % Keep only HBO as similar results should apply to HbR
job = nirs.modules.LabelShortSeperation(job);

% Data is simulated with relationship between channels at either zeroth lag
% or first lag
laglabels  = {'Zeroth lag','First lag'};
lags={[1,0],[0,1]}; % [1 0] indicates zeroth lag and [0 1] indicates first lag

for lag =1:length(lags)
     mkdir([Curr_Folder filesep 'Lag' filesep laglabels{lag}]);
     cd([Curr_Folder filesep 'Lag' filesep laglabels{lag}]);
      
% Generate the simulated ground truth at either the zeroth or first lag. This ground truth will be used to
% evaluate the performance of several analysis pipelines
[raw , truth] = simData_connectivity_shortsep([],lags{lag},0); 
hboRest = job.run(raw);

% Setting simulation Parameters
nIters =100; % No of iterations, no of datasets to be generated
Fs = raw.Fs; %Sampling rate
Pmax = round(10*Fs); % Max model order for autoregressive (AR) models
nchan = sum(~raw.probe.link.ShortSeperation)/2; %No of long channels
nshort= sum(raw.probe.link.ShortSeperation)/2; %No of short channels
ntp = size(raw.data, 1);  %No of timepoints in the data
criterion = 'BIC'; %Informtation criteria to compare model fits to find optimal model orders in AR and MVAR models 
truth_lag = squeeze(truth(:,:,lag));
zeroind = triu(true(size(truth_lag)),1)&(~truth_lag); %Indices of elements with true negatives in the data covariance  matrix
dataind = triu(true(nchan),1); % getting the upper right hand corner of the data covariance matrix.


pathsa = nnz(zeroind); %True negative paths
pathsb = 0.5*nchan*(nchan-1);  % All paths

% Truth distribution
% This will later be used to compare the performance of simulation results
t = sqrt(ntp-Pmax).*((truth_lag)./(sqrt(1-truth_lag.^2)));
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
datalabels = {'random data','w temp corr noise', ....
    'w temp & spatial corr noise', 'w temp & spatial corr noise & motion',};


for type=3:4 % Just running for temporal autocorrelated+ spatial correlated noise  
     % & with motion should be sufficient as both AR partial correlation methods
     % & MVGC methods correct for temporally correlated noise
     
  mkdir([Curr_Folder filesep 'Lag' filesep laglabels{lag} filesep datalabels{type}])
  cd([Curr_Folder filesep 'Lag' filesep laglabels{lag} filesep datalabels{type}])
   
%% Running ths simulation
    
% Initializing variables to store values obtained from the analysis pipeline with zeros

Rdist_3a = zeros(pathsa,nIters); Rdist_31a = zeros(pathsa,nIters);
Zdist_3a = zeros(pathsa,nIters); Zdist_31a = zeros(pathsa,nIters);
pdist_3a = zeros(pathsa,nIters); pdist_31a = zeros(pathsa,nIters);

Fdist_6a = zeros(pathsa,nIters); Fdist_61a = zeros(pathsa,nIters);
Zdist_6a = zeros(pathsa,nIters); Zdist_61a = zeros(pathsa,nIters);
pdist_6a = zeros(pathsa,nIters); pdist_61a = zeros(pathsa,nIters);
Gdist_6a = zeros(pathsa,nIters); Gdist_61a = zeros(pathsa,nIters);

Fdist_7a = zeros(pathsa,nIters); Fdist_71a = zeros(pathsa,nIters);
Zdist_7a = zeros(pathsa,nIters); Zdist_71a = zeros(pathsa,nIters);
pdist_7a = zeros(pathsa,nIters); pdist_71a = zeros(pathsa,nIters);
Gdist_7a = zeros(pathsa,nIters); Gdist_71a = zeros(pathsa,nIters);


Rdist_3b = zeros(pathsb,nIters); Rdist_31b = zeros(pathsb,nIters); 
Zdist_3b = zeros(pathsb,nIters); Zdist_31b = zeros(pathsb,nIters); 
pdist_3b = zeros(pathsb,nIters); pdist_31b = zeros(pathsb,nIters); 


Fdist_6b = zeros(pathsb,nIters); Fdist_61b = zeros(pathsb,nIters);
Zdist_6b = zeros(pathsb,nIters); Zdist_61b = zeros(pathsb,nIters);
pdist_6b = zeros(pathsb,nIters); pdist_61b = zeros(pathsb,nIters);
Gdist_6b = zeros(pathsb,nIters); Gdist_61b = zeros(pathsb,nIters);

Fdist_7b = zeros(pathsb,nIters); Fdist_71b = zeros(pathsb,nIters);
Zdist_7b = zeros(pathsb,nIters); Zdist_71b = zeros(pathsb,nIters);
pdist_7b = zeros(pathsb,nIters); pdist_71b = zeros(pathsb,nIters);
Gdist_7b = zeros(pathsb,nIters); Gdist_71b = zeros(pathsb,nIters);

% Running several iterations
 for iter=1:nIters  
   
   % Simulate the data based on statisical properites

    switch type-1
        case 1
          % With temporal autocorrelated noise 
          [raw1 , truthx] = simData_connectivity_shortsep(truth,lags{lag},0);
          data = job.run(raw1);
 

        case 2
           % With termporal autocorrelation & spatially correlated physiological noise
          [raw1 , truthx] = simData_connectivity_shortsep(truth,lags{lag},150);
          data = job.run(raw1);

        case 3
          % WIth motion artifacts & spatially correlated physiological noise
          [raw1 , truthx] = simData_connectivity_shortsep(truth,lags{lag},150);
          raw1 = nirs.testing.simMotionArtifact(raw1); % adding motion artifacts
          data = job.run(raw1);

        otherwise
          % Random noise
          [raw1 , truthx] = simData_connectivity_shortsep(truth,lags{lag},0);     
          data = job.run(raw1);
          data.data = rand(ntp,nchan+nshort);
    end

  

    % AR Partial correlation with PCA filtering of short channels with robust
    % regression
    [R3, p3, dfe3] = partial_corr(data,'10xFs','short',true, true);
    Z3= atanh(R3);
    Rdist_3a(:, iter) = R3(zeroind);
    Zdist_3a(:, iter) = Z3(zeroind);
    pdist_3a(:, iter) = p3(zeroind);
    Rdist_3b(:, iter) = cat(1,R3(nnzind),R3(zeroind));
    Zdist_3b(:, iter) = cat(1,Z3(nnzind),Z3(zeroind));
    pdist_3b(:, iter) = cat(1,p3(nnzind),p3(zeroind));
    
    % AR Partial correlation with PCA filtering of short channels without robust
    % regression
    [R31, p31, dfe31] = partial_corr(data,'10xFs','short',false, true);
    Z31= atanh(R31);
    Rdist_31a(:, iter) = R31(zeroind);
    Zdist_31a(:, iter) = Z31(zeroind);
    pdist_31a(:, iter) = p31(zeroind);
    Rdist_31b(:, iter) = cat(1,R31(nnzind),R31(zeroind));
    Zdist_31b(:, iter) = cat(1,Z31(nnzind),Z31(zeroind));
    pdist_31b(:, iter) = cat(1,p31(nnzind),p31(zeroind));
    
    % Modified MVGC with zero-lag after PCA filtering of short channels
    % with robust regression
     [G6, F6, df16, df26, p6] = grangercausality(data,Pmax,'multivariate','short',true,true,true,criterion);
     Z6 = log( F6 ) / 2;
     mu6 = 0.5*((1/df16)-(1/df26));
     sd6= sqrt(0.5*((1/df16)+(1/df26)));
     Z6 = ((Z6-mu6)/sd6);
     Gdist_6a(:, iter) = G6(zeroind);
     Fdist_6a(:, iter) = F6(zeroind);
     Fdist_6a(:, iter) = Z6(zeroind);
     pdist_6a(:, iter) = p6(zeroind);
     Gdist_6b(:, iter) = cat(1,G6(nnzind),G6(zeroind));
     Fdist_6b(:, iter) = cat(1,F6(nnzind),F6(zeroind));
     Zdist_6b(:, iter) = cat(1,Z6(nnzind),Z6(zeroind));
     pdist_6b(:, iter) = cat(1,p6(nnzind),p6(zeroind));
     
    % Modified MVGC with zero-lag after PCA filtering of short channels
    % without robust regression
     [G61, F61, df161, df261, p61] = grangercausality(data,Pmax,'multivariate','short',true,false,true,criterion);
     Z61 = log( F61 ) / 2;
     mu61 = 0.5*((1/df161)-(1/df261));
     sd61= sqrt(0.5*((1/df161)+(1/df261)));
     Z61 = ((Z61-mu61)/sd61);
     Gdist_61a(:, iter) = G61(zeroind);
     Fdist_61a(:, iter) = F61(zeroind);
     Fdist_61a(:, iter) = Z61(zeroind);
     pdist_61a(:, iter) = p61(zeroind);
     Gdist_61b(:, iter) = cat(1,G61(nnzind),G61(zeroind));
     Fdist_61b(:, iter) = cat(1,F61(nnzind),F61(zeroind));
     Zdist_61b(:, iter) = cat(1,Z61(nnzind),Z61(zeroind));
     pdist_61b(:, iter) = cat(1,p61(nnzind),p61(zeroind));

    
    % MVGC after PCA filtering of short channels with robust regression
     [G7, F7, df17, df27, p7] = grangercausality(data,Pmax,'multivariate','short',false,true,true,criterion);
     Z7 = log( F7 ) / 2;
     mu7 = 0.5*((1/df17)-(1/df27));
     sd7= sqrt(0.5*((1/df17)+(1/df27)));
     Z7 = ((Z7-mu7)/sd7);
     Gdist_7a(:, iter) = G7(zeroind);
     Fdist_7a(:, iter) = F7(zeroind);
     Fdist_7a(:, iter) = Z7(zeroind);
     pdist_7a(:, iter) = p7(zeroind);
     Gdist_7b(:, iter) = cat(1,G7(nnzind),G7(zeroind));
     Fdist_7b(:, iter) = cat(1,F7(nnzind),F7(zeroind));
     Zdist_7b(:, iter) = cat(1,Z7(nnzind),Z7(zeroind));
     pdist_7b(:, iter) = cat(1,p7(nnzind),p7(zeroind));
     
     % MVGC after PCA filtering of short channels without robust regression
     [G71, F71, df171, df271, p71] = grangercausality(data,Pmax,'multivariate','short',false,false,true,criterion);
     Z71 = log( F71 ) / 2;
     mu71 = 0.5*((1/df171)-(1/df271));
     sd71= sqrt(0.5*((1/df171)+(1/df271)));
     Z71 = ((Z71-mu71)/sd71);
     Gdist_71a(:, iter) = G71(zeroind);
     Fdist_71a(:, iter) = F71(zeroind);
     Fdist_71a(:, iter) = Z71(zeroind);
     pdist_71a(:, iter) = p71(zeroind);
     Gdist_71b(:, iter) = cat(1,G71(nnzind),G71(zeroind));
     Fdist_71b(:, iter) = cat(1,F71(nnzind),F71(zeroind));
     Zdist_71b(:, iter) = cat(1,Z71(nnzind),Z71(zeroind));
     pdist_71b(:, iter) = cat(1,p71(nnzind),p71(zeroind));
   fprintf("Iteration completed: %d of %d\n", iter, nIters);
 end
 
 %% Calculating the necessary metrics
 
 % Converting the r values obtrained from anlaysis pipelines into column vector for later plotting of
 % bootstrap -p vs actual p

 Rdist_3a = Rdist_3a(:); Zdist_3a = Zdist_3a(:); pdist_3a = pdist_3a(:); 
 Rdist_31a = Rdist_31a(:); Zdist_31a = Zdist_31a(:); pdist_31a = pdist_31a(:);
 Gdist_6a = Gdist_6a(:); Fdist_6a = Fdist_6a(:); pdist_6a = pdist_6a(:); 
 Gdist_61a = Gdist_61a(:); Fdist_61a = Fdist_61a(:); pdist_61a = pdist_61a(:);  
 Gdist_7a = Gdist_7a(:); Fdist_7a = Fdist_7a(:); pdist_7a = pdist_7a(:); 
 Gdist_71a = Gdist_71a(:); Fdist_71a = Fdist_71a(:); pdist_71a = pdist_71a(:);  

 [R_sorted_3, indxR3]= sort(Rdist_3a, 'ascend'); [R_sorted_31, indxR31]= sort(Rdist_31a, 'ascend'); 
 [G_sorted_6, indxG6]= sort(Gdist_6a, 'ascend'); [G_sorted_61, indxG61]= sort(Gdist_61a, 'ascend'); 
 [G_sorted_7, indxG7]= sort(Gdist_7a, 'ascend'); [G_sorted_71, indxG71]= sort(Gdist_71a, 'ascend'); 
 
 % Calculations for generating ROC curves and calculating AUC
  indx1 =  randperm(tp);
  indx2 =  randperm(tn);
  
% Since we have less true positives, we randomly sample a subset of true
% negatives to plot ROC curves which require equal number of true positives
% & negatives
 pdist_3b = cat(1,pdist_3b(indx1(1:nvals),:),pdist_3b(indx2(1:nvals),:));
 pdist_31b = cat(1,pdist_31b(indx1(1:nvals),:),pdist_31b(indx2(1:nvals),:)); 
 pdist_6b = cat(1,pdist_6b(indx1(1:nvals),:),pdist_6b(indx2(1:nvals),:));
 pdist_61b = cat(1,pdist_61b(indx1(1:nvals),:),pdist_61b(indx2(1:nvals),:));  
 pdist_7b = cat(1,pdist_7b(indx1(1:nvals),:),pdist_7b(indx2(1:nvals),:));
 pdist_71b = cat(1,pdist_71b(indx1(1:nvals),:),pdist_71b(indx2(1:nvals),:)); 
 
 pdist_3b = pdist_3b(:); pdist_31b = pdist_31b(:); 
 pdist_6b = pdist_6b(:); pdist_61b = pdist_61b(:); 
 pdist_7b = pdist_7b(:); pdist_71b = pdist_71b(:); 
 
 % ROC curve analysis
[tp_p3,fp_p3,th_p3] = nirs.testing.roc(truth_all_p, pdist_3b);
[tp_p31,fp_p31,th_p31] = nirs.testing.roc(truth_all_p, pdist_31b); 
[tp_p6,fp_p6,th_p6] = nirs.testing.roc(truth_all_p, pdist_6b);
[tp_p61,fp_p61,th_p61] = nirs.testing.roc(truth_all_p, pdist_61b); 
[tp_p7,fp_p7,th_p7] = nirs.testing.roc(truth_all_p, pdist_7b);
[tp_p71,fp_p71,th_p71] = nirs.testing.roc(truth_all_p, pdist_71b);

% AUC metrics 
AUC_p3 = trapz(fp_p3, tp_p3); AUC_p31 = trapz(fp_p31, tp_p31); 
AUC_p6 = trapz(fp_p6, tp_p6); AUC_p61 = trapz(fp_p61, tp_p61); 
AUC_p7 = trapz(fp_p7, tp_p7); AUC_p71 = trapz(fp_p71, tp_p71); 

 
%% Plotting Figures

% Type-1 Error control plot (Actual p vs bootstraped-p)
titlea= sprintf("FP plot for %s",datalabels{type});
figure('Name',titlea,'NumberTitle','off','Position',[100 100 800 600]);

plot(quantile(pdist_3a(indxR3),int), pvals_int, 'b','LineWidth',2); hold on;
plot(quantile(pdist_31a(indxR31),int), pvals_int, 'b--','LineWidth',2);
plot(quantile(pdist_6a(indxG6),int), pvals_int, 'm','LineWidth',2); hold on;
plot(quantile(pdist_61a(indxG61),int), pvals_int, 'm--','LineWidth',2);
plot(quantile(pdist_7a(indxG7),int), pvals_int, 'c','LineWidth',2); hold on;
plot(quantile(pdist_71a(indxG71),int), pvals_int, 'c--','LineWidth',2); hold on;

xlabel("p-hat"); ylabel(" False positive rate"); pbaspect([1 1 1]);
l =legend(["AR Par-corr w/ robust", "AR Par-corr", ...
    "Modified MVGC w/ robust", "Modified MVGC",...
    "MVGC w/ robust", "MVGC"], 'Location','eastoutside');
l.FontSize = 15;
% Plotting a x=y ref line to compare the false positive(FP) rate vs expected FP rate 
hline = refline([1,0]);
hline.LineStyle = ':';
hline.LineWidth = 2;
hline.Color = 'k';
hline.HandleVisibility = 'off';
ax =gca;
ax.FontSize = 18;
ax.FontWeight = 'bold';
saveas(gcf,titlea,'tiffn'); 
saveas(gcf,titlea,'mfig');


% Plotting the ROC curves to assess performance
titleb= sprintf("ROC curves for %s",datalabels{type});
figure('Name',titleb,'NumberTitle','off','Position',[100 100 800 600]);
 
plot(fp_p3,tp_p3, 'b','LineWidth',2); hold on;
plot(fp_p31,tp_p31, 'b--','LineWidth',2);
plot(fp_p6,tp_p6, 'm','LineWidth',2); hold on;
plot(fp_p61,tp_p61, 'm--','LineWidth',2);
plot(fp_p7,tp_p7, 'c','LineWidth',2); hold on;
plot(fp_p71,tp_p71, 'c--','LineWidth',2); hold on;

xlabel('False positive rate'); ylabel('True positive rate'); pbaspect([1 1 1]);
l = legend(["AR Par-corr w/ robust", "AR Par-corr", ...
    "Modified MVGC w/ robust", "Modified MVGC",...
    "MVGC w/ robust", "MVGC"], 'Location','eastoutside');
l.FontSize = 15;
% Plotting a x=y ref line to show performance obtained by random chance
hline = refline([1,0]);
hline.LineStyle = ':';
hline.LineWidth = 2;
hline.Color = 'k';
hline.HandleVisibility = 'off';
ax =gca;
ax.FontSize = 18;
ax.FontWeight = 'bold';
saveas(gcf,titleb,'tiffn'); 
saveas(gcf,titleb,'mfig');
cd(Curr_Folder)
end
end
%close all;


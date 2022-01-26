%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% CODE TO BUILD THE MULTI-MODAL AND MULTI-SUBJECT NETWORK %%%%%%%%%%%%%
%%%%% AND RUN THE EXTENDED MULTI-LAYER MODULARITY OPTIMIZATION %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars
clc

%% set paths and directories

% add to the path the genLouvain package
addpath(genpath('...\GenLouvain-2.1'));

% folder with preprocessed sc and fc matrices
dataDir = '...'; 

% create a folder where we'll save the output of the extended multi-layer modularity optimization
outputDir = '...'; 
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end

%% load data

load(fullfile(dataDir, 'SCFC_schaefer100.mat'));

sc = data.sc;  % NB: sc is already in the form of a correlation matrix
fc = data.fc;
cellage = data.age;
clear data

% select subjects between 20 and 40 yo
age = zeros(length(cellage),1);
for i=1:length(cellage)
    age(i) = str2num(cellage{i});
end
clear data

adults = find(age>20 & age<40);
sc = sc(:,:,adults); 
fc = fc(:,:,adults);


minW = min(sc(sc>0)); % minimum weight of the matrices

N = size(fc, 1);    % nodes
T = size(fc, 3);    % subjects

%% 

% set parameters values for genLouvain

etaVals = linspace(minW, 1, 50);
omegaVals = linspace(minW, 1, 50);
gammaVals = linspace(minW, 1, 50);
limit = 10000;
verbose = false;

% build the multi-modal and multi-subject modularity matrix
B = spalloc(N*T*2, N*T*2, (N+T*2)*N*T*2);
gammaMask = B;
gammaMaskWeighted = B;
etaMask = B;
omegaMask = B;
NN = N*2;
for t = 1:T
    idx = (1:N) + (t - 1)*N*2;
    jdx = idx + N;
    B(idx,idx) = sc(:,:,t);
    B(jdx,jdx) = fc(:,:,t);
    
    
    % gammaMask: uniform null model
    gammaMask(idx,idx) = ones(N).*~eye(N);
    gammaMask(jdx,jdx) = ones(N).*~eye(N);
   
    
    %etaMask
    etaMask(idx,jdx) = eye(N);
    etaMask(jdx,idx) = eye(N);
    
    % omegaMask: all-to-all
    for s = 1:T
        idx = (1:NN) + (s - 1)*NN;
        jdx = (1:NN) + (t - 1)*NN;
        if s ~= t
            omegaMask(idx,jdx) = eye(NN);
        end
    end
    
    clear idx jdx kfc ksc divfc divsc
end

% loop over gamma, omega, eta values
for igamma = 1:length(gammaVals)
    gamma = gammaVals(igamma);
    
    for iomega = 1:length(omegaVals)
        omega = omegaVals(iomega);
        
        for ieta = 1:length(etaVals)
            sprintf('g = %g, w = %g, e = %g',...
                gammaVals(igamma),omegaVals(iomega),etaVals(ieta))
            eta = etaVals(ieta);
            
            r = randperm(N*T*2);
            [~, idx] = sort(r);
            
            Bfull = (B - (gammaMask*gamma)) + (etaMask*eta) + (omegaMask*omega);
            [ci, Q(iomega,ieta)] = genlouvain(Bfull(r,r), limit, verbose);
            ci = ci(idx);
            
            % reshape communities
            ci = reshape(ci, [N,T*2]);
            cir(:,:,iomega,ieta) = ci;
            cifc(:,:,iomega,ieta) = ci(:,2:2:end);
            cisc(:,:,iomega,ieta) = ci(:,1:2:end);
            
            clear ci eta r idx Bfull
        end
        clear omega ieta
    end
    clear gamma iomega
    
    % save partitions
    save(fullfile(outputDir, sprintf('CommFC_adults_g%g', igamma)), 'cifc', '-v7.3')
    save(fullfile(outputDir, sprintf('CommSC_adults_g%g', igamma)), 'cisc', '-v7.3')
    save(fullfile(outputDir, sprintf('Q_adults_g%g', igamma)), 'Q', '-v7.3')
end
clear igamma


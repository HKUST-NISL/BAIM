%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                   BAIM code for paper's graph                       %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load GASSOM basis functions and example binocular images
clear all;
load('example-data/GASSOMbases.mat');
load('example-data/IcubImage4.mat');
[MM,NN] = size(imL);
histwin_half = 20;
binslocal = 55;
binsglobal = round(sqrt(231*311));

subspaceind = 287;
tempbasis1 = GASSOMbases{1}(:,subspaceind);
tempbasis2 = GASSOMbases{2}(:,subspaceind);
% tempbasis1 = GASSOMbases{1};
% tempbasis2 = GASSOMbases{2};


%% Show the coresponding coefficient map,histogram and saliency map for each subspace

[infomap_LBAIM,response,histo] = LBAIM(imL,imR,size(imL)/2,tempbasis1,tempbasis2);
infomap_LBAIM = (infomap_LBAIM-min(min(infomap_LBAIM)))/(max(max(infomap_LBAIM))-min(min(infomap_LBAIM)));
figure();imagesc(squeeze(response));hold on;plot(histwin_half*[-1,-1,1,1,-1]+160,histwin_half*[1,-1,-1,1,1]+120,'r');axis equal tight off;
figure();bar(linspace(0,1,length(histo)),histo,0.6);axis([-0.01 find(histo<0.005*max(histo),1)/length(histo) 0 max(histo)]);
figure();imagesc(infomap_LBAIM,[0 1]);axis equal tight off;

[infomap_GBAIM,response,histo] = GBAIM(imL,imR,tempbasis1,tempbasis2);
infomap_GBAIM = (infomap_GBAIM-min(min(infomap_GBAIM)))/(max(max(infomap_GBAIM))-min(min(infomap_GBAIM)));
figure();imagesc(squeeze(response));axis equal tight off;
figure();bar(linspace(0,1,length(histo)),histo,0.8);axis([-0.01 find(histo<0.005*max(histo),1)/length(histo) 0 max(histo)]);
figure();imagesc(infomap_GBAIM,[0 1]);axis equal tight off;


%% Show the 2 basis vectors in the corresponding sub-space
basisL1 = reshape(tempbasis1(1:100),[10,10]);
basisR1 = reshape(tempbasis1(101:200),[10,10]);
tempshowbasis1 = [basisL1;-1*ones(1,10);basisR1];

basisL2 = reshape(tempbasis2(1:100),[10,10]);
basisR2 = reshape(tempbasis2(101:200),[10,10]);
tempshowbasis2 = [basisL2;-1*ones(1,10);basisR2];

basismax = max([tempbasis1;tempbasis2]);
basismin = min([tempbasis1;tempbasis2]);

totalshowbasis = ones(23,23);
totalshowbasis(2:end-1,2:11) = tempshowbasis1;
totalshowbasis(2:end-1,13:end-1) = tempshowbasis2;

figure();imagesc(totalshowbasis,[basismin,basismax]);colormap(gray);
set(gcf,'position',[500,400,200,200]);axis equal tight off;
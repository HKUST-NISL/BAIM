%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                                                                     %%%
%%% GBAIM (global binocular attention based on information maximization)%%%
%%%                                                                     %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [infomap,response,histo] = GBAIM(imL,imR,bases1,bases2,convolve)

% Set to defaults or input values
nargs = nargin;
if nargs < 6, convolve = 1; end

% For convolve = 1 set these
sigval = 5; % How many pixels correspond to 1 degree visual angle
sigwin = [30 30]; % What size of window is needed to contain the above

bins = 100;

% Output related
contrastval = 2; % How much contrast in the "transparent" representation

% Read image
imL = im2double(imL);
imR = im2double(imR);

p = 10;  pm = p-1;  ph = p/2;   
                 

col_imL = im2col(imL,[p p],'sliding');
col_imR = im2col(imR,[p p],'sliding');
col_bino = [col_imL;col_imR];
response = reshape((bases1'*col_bino).^2+(bases2'*col_bino).^2,[size(bases1',1),size(imL,1)-p+1,size(imL,2)-p+1]);

minscale = min(min(min(response)));
maxscale = max(max(max(response)));

for z=1:size(response,1)
   % Translate image from 1xMxN to MxN for ease of processing
   tempim = zeros(size(imL,1)-pm,size(imL,2)-pm);
   tempim(1:size(imL,1)-pm,1:size(imL,2)-pm) = response(z,1:size(imL,1)-pm,1:size(imL,2)-pm); 

   
   % Scale temporary copy of feature plane
   stempim = tempim-minscale;
   stempim = stempim./(maxscale-minscale);
   

   % Compute histogram
   histo = imhist(stempim(),bins);
   % Rescale values based on histogram to reflect likelihood
   ts(z,:,:) = histo(round(stempim*(bins-1)+1))./sum(histo);
end


% Overall information content is product of individual feature maps with
% Shannon definition of Self-Information applied - Can do a log of products
% or a sum of logs...

% Width and height of information map
wid = size(ts,2);
hei = size(ts,3);

% Initialize to information gained from 1st feature domain
infomapt = -log(reshape(ts(1,:,:),wid,hei)+0.000001);

% Add information gained from remaining features
for z=2:size(response,1)
infomapt = infomapt - log(reshape(ts(z,:,:),wid,hei)+0.000001);
end


if (convolve == 1)
    infomapt = filter2(fspecial('gaussian',sigwin,sigval),infomapt);
end

infomapt = infomapt.^contrastval;

% Pad the final map so that its size matches the input image
infomap = zeros(size(imL,1),size(imL,2))+min(min(infomapt));
infomap(ph+1:size(imL,1)-ph+1,ph+1:size(imL,2)-ph+1)=infomapt;

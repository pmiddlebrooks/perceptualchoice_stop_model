

addpath(genpath('/Users/bramzandbelt/Dropbox/Code/'))
hFig = set_figure({982,663,'pixels'},{'USLetter','landscape'});
modelVar = cell2mat(arrayfun(@(in1) find(model==in1),[1,5,13,2,7,6,33,10,21,20,23,35,34,12,11],'Uni',0));


p = panel;
p.margin = [5 5 2 2];

p.pack(6,15);


CT = othercolor('Greens3',64);
CT = flipud(CT);
colormap(CT);
for iSubject = 1:6
  
  % Scale data by subject
  
  y = Y(:,:,iSubject);
  
  % Compute ranking
  [srtY,iY] = sort(y(:));
  rnkY = nan(size(y));
  rnkY(iY(1:9)) = 1:9;
  
  y = y - min(y(:));
  y(y < 1) = 1;
  y = log10(y);
  
  
  
  mn = min(y(:));
  
  % Set cut off at 3 (i.e. dBic = 1000)
  y(y > 3) = 3;
  
  
  rng = max(y(:))-mn;
  y = 1+63*(y-mn)/rng; % Self scale data
  

  for iModel = 1:15
    
    thisY = y(modelVar(iModel),:);
    thisY = reshape(thisY,3,3);
    
    p(iSubject,iModel).select();
    sc(thisY,CT,[1,64],'w')
    
%     image(thisY);
    
%     thisImage = pcolor([thisY,nan(3,1);nan(1,3+1)]);
%     
%     heatmaptext(rand(3,3),'fontcolor','k','precision',2);
%     
    
    thisRnk = rnkY(modelVar(iModel),:);
    thisRnk = reshape(thisRnk,3,3);

    iRnk = find(~isnan(thisRnk));
    rnk = thisRnk(iRnk);
    [iRow,iCol] = ind2sub([3,3],iRnk);

%     iRow = [1 2 3 1 2 3 1 2 3];iCol = [1 1 1 2 2 2 3 3 3];
    arrayfun(@(in1,in2,in3) text(in1,in2,num2str(in3),'HorizontalAlignment','Center'),iCol,iRow,rnk);
%     arrayfun(@(in1,in2,in3) text(in1,in2,num2str(in3),'HorizontalAlignment','Center'),iCol(:),iRow(:),thisRnk(:))
    
    axis ij
    axis tight
    axis square
    set(gca,'XTick',[],'YTick',[])
    
  end
end
% 
% p(2).select();
% 
% hC = colorbar;
% 
% L = [1 10 100 1000];
% 
% % Choose appropriate
% 
% % or somehow auto generate colorbar labels
% 
% l = 1+63*(log10(L)-mn)/rng; % Tick mark positions
% 
% set(hC,'Ytick',l,'YTicklabel',L);
% 
% 
% 
% % showPixelValues; 
% heatmapDemo
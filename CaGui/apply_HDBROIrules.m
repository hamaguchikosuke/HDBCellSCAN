function [stat, res, ops] = apply_HDBROIrules(ops, stat0, res0,clustrules)

% because HDBCellScan only detects spatially connected pixels, 
% no need for further dividing. 

% divide cluster into subregions.
% stat0    = get_regions_kh(stat0, res0);

%% I can add more filter function here
% if ops.splitROIs
%     [stat, res] = get_validregions(stat0,res0, clustrules);
% else
   stat = stat0;
   res = res0;
   
   
   for ii = 1:length(stat)
      if (stat(ii).npix > clustrules.MinNpix) & (stat(ii).npix <= clustrules.MaxNpix)
          stat(ii).igood = 1;
      else
          stat(ii).igood= 0;
      end
   end
% end




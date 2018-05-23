function [ops1,U,Sv] = get_svd_20171102(ops1)

for ii = 1:numel(ops1)   

    if numel(ops1{ii}.yrange)<10 || numel(ops1{ii}.xrange)<10
        warning('valid range after registration very small, continuing to next plane')
%         continue;
    end
    
    if getOr(ops1{ii}, {'getSVDcomps'}, 0)
        ops1{ii}    = get_svdcomps(ops1{ii});
    end
    
    if ops1{ii}.getROIs || getOr(ops1{ii}, {'writeSVDroi'}, 0)
        [ops1{ii}, U, Sv]    = get_svdForHDBSCAN(ops1{ii}); % -> SVDroi_<mouse_name>_<date>_plane<#>_ch<#>
    end
   
end
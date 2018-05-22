function ops1 = get_roi_20171102(ops1,U,Sv,clustrules)

clustModel=ops1.clustModel;
clustrules = get_clustrules(clustrules);

if ops1.getROIs
    switch clustModel
        case 'HDBCellScan' 
            % control the dimension of nSVD here
            
            [Ly,Lx,nSVD]=size(U);
            
            U = U(:, :, 1:min(nSVD,ops1.nSVDforROI));
            Sv              = Sv(1:min(nSVD,ops1.nSVDforROI));
              
            % smooth svd 
            U=single(imgaussfilt(double(U),1.5));
            [ops1, stat, res] = get_HDBCellScan(ops1, U, Sv,clustrules);
                  
        case 'standard'
            [ops1, stat, res]  = fast_clustering(ops1,U, Sv);
        case 'neuropil'
            
            ops1.diameter = getOr(ops1,'diameter',clustrules.diameter);
            [ops1, stat, res]  = fastClustNeuropilCoef_kh(ops1,U, Sv); % save data in '%s/F_%s_%s_plane%d_Nk%d.mat'
    end
    
    non_empty_cluster = find([stat.npix]>0);
    stat = stat(non_empty_cluster);
    
    switch clustModel
        case 'HDBCellScan'
            % in case I want to add more filter in here
            stat0 = stat;
            res0  = res;
            [stat, res,ops1] = apply_HDBROIrules(ops1, stat, res,clustrules);
           
                      
          %% Save the ROI and related information.
            if ~exist(ops1.ResultsSavePath, 'dir'),  mkdir(ops1.ResultsSavePath),    end
            
            MinClustSize = clustrules.MinClust;
            ops1.ROISaveFile = sprintf('%s/ROI_%s_%s_plane%d_Ch%d_MinClust%d.mat', ...
            ops1.ResultsSavePath, ...
                ops1.mouse_name, ops1.date, ops1.PlaneID,ops1.ChannelID, MinClustSize);
            ops1.ROISummaryName = sprintf('%s/ROI_%s_%s_%s_plane%d_Ch%d_MinClust%d.mat', ops1.ResultsSavePath, ...
                ops1.mouse_name, ops1.date, ops1.SubDirs{1},ops1.PlaneID,ops1.ChannelID, MinClustSize);
            
            ops = ops1;
            fprintf('Saving ROI in %s\n',ops1.ROISaveFile);
            save(ops.ROISaveFile,  'ops', 'res', 'stat', 'stat0', 'res0', 'clustrules');

           %% show beautiful ROIs 
           
            myfigure('ROI');cla;
            opt.hue_range = [0 1];
            IMG=ROI_gem_img(res.iclust,res.probabilities,res.M,opt);
            
            imshow(IMG);
            drawnow;
            [savepath,savename,~]=fileparts(ops1.ROISummaryName);
            savepath = fullfile(ops1.RootStorage,ops1.mouse_name,'Summary');
            if ~exist(savepath,'dir'), mkdir(savepath); end
            title(savename,'FontSize',9,'Interpreter','none');
            
            print(gcf,'-dpng',fullfile(savepath,[savename,'.png']));

        case {'standard', 'neuropil'}

            [~, ~,ops1] = apply_ROIrules_kh(ops1, stat, res, clustrules);
            %    inside this code, ROI is saved in -> 
            % ops.ROISaveFile = sprintf('%s/F_%s_%s_plane%d_Ch%d_Nk%d.mat', ops.ResultsSavePath, ...
            %     ops.mouse_name, ops.date, ops.PlaneID,ops.ChannelID, Nk);
            % save(ops.ROISaveFile,  'ops', 'res', 'stat', 'stat0', 'res0', 'clustrules')
    end

end


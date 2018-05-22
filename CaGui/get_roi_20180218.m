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
            
            U=reshape(U,[],ops1.nSVDforROI);
           
            do_subtract_neuropil = 1; % 
            % well, I need to do something to subtract neuropil for deeper imaging, or
            % higher expression brain.
            
            if (do_subtract_neuropil)
               ops1.diameter = getOr(ops1,'diameter',clustrules.diameter);
               ops1.ratioNeuropil = getOr(ops1,'ratioNeuropil',10); % 10 times 
               
               % S is normalized npix x Tile^s matrix
                S = getNeuropilBasis(ops1, Ly, Lx, 'Fourier'); % 'raisedcosyne', 'Fourier'
                [npix,nTilesSq]=   size(S);     
                nBasis = sqrt(nTilesSq) ;
                
                %% method 1: U = U + neu; neu = S'*lambda 
                StS = S'*S;
                StU = S'*U;
                lambda = (StS+1e-4*eye(nTilesSq))\StU;
                neu = S*lambda;
                U = U-neu; 
                
%                 %% 
%                 PixL = ones(1,npix)';
%                 
%                 Uneu = U';
%                 
%                 Sm = bsxfun(@times, S, PixL);
%                 StS = Sm' * Sm;
%                 StU = Sm' * Uneu';
%                 Lam = (StS + 1e-4 * eye(nBasis)) \ StU;
%                 
%                 % recompute neuropil pixel contribution
%                 neuropil = Lam' * S';
%                 PixL = mean(bsxfun(@times, neuropil, Uneu), 1);
%                 PixL = bsxfun(@rdivide, PixL, mean(neuropil.^2,1)); % normalize variance to 1
%                 PixL = max(0, PixL);
%                 neuropil = bsxfun(@times, neuropil, PixL);
%                 U = U - neuropil'; %what's left over for cell model
%                 
                %     U = reshape(U,[],nSVD);
                fig_title = 'U-neuropil';
                Sneuropil = S;
            else
                Sneuropil = [];
                fig_title = 'U-only';
            end
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

           %% print a beautiful picture of detected cells
%             L=res.iclust;
%             [Ly,Lx]=size(L);
%             L(L==0)=length(unique(L))+1;
%             
%             rng('default');
%             rand_hue  = rand(1,length(unique(L)));
%             rng('shuffle');
%             rand_hue= [rand_hue,0];
%             Hue = rand_hue(L);
%             Sat = reshape(res.probabilities,Ly,Lx);
%             Val = res.M;
%             IMG=hsv2rgb(cat(3,Hue,Sat,Val));
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


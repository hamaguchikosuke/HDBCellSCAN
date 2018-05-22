% Density-Based Spatial Clustering of Applications with Noise (DBSCAN) 
% Distance matrix version 
% 
% Usage:  [IDX, isnoise]=dbscan_D(D,epsilon,MinPts);
% === inputs ===
% D: distance matrix. It needs to be a symmetric square matrix.
% 
% epsilon: cut-off distance to be a cluster
% MinPts:  minimum points. Core points need to have at least MinPts of
% neighbors within epsilon distance.
% 
% modified from DBSCAN by Yarpiz so that it can analyze larger data set
% 
% % Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% % All rights reserved. Please read the "license.txt" for license terms.
%
% % Project Code: YPML110
% % Project Title: Implementation of DBSCAN Clustering in MATLAB
% % Publisher: Yarpiz (www.yarpiz.com)
% % 
% % Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% % 
% % Contact Info: sm.kalami@gmail.com, info@yarpiz.com
% 
% 
% modified by Kosuke Hamaguchi, 2017/09/09

function [IDX, isnoise]=dbscan_D(D,epsilon,MinPts)

    C=0;
    
    n=size(D,1);
    IDX=zeros(n,1);
    
       
    visited=false(n,1);
    isnoise=false(n,1);
    
    percent_done =0;
    for ii=1:n
        if ~visited(ii)
            visited(ii)=true;
            
            Neighbors_vec=RegionQuery(ii);
            if nnz(Neighbors_vec)<MinPts
                % X(i,:) is NOISE
                isnoise(ii)=true;
            else
                C=ceil(C+1);
                ExpandCluster(ii,Neighbors_vec,C);
            end
            
        end
      if ii/n > percent_done+0.1
          fprintf('%d/%d,',ii,n);
          percent_done= percent_done+0.1;
      end
    end
    fprintf('Done\n');
    
    function ExpandCluster(ii,Neighbors_vec,C)
        IDX(ii)=C;
        
        kk = 1;
        
        search_index =find(Neighbors_vec);  
        
        while true
            jj = search_index(kk);
            
            if ~visited(jj)
                visited(jj)=true;
                NewNeighbors_vec=RegionQuery(jj);
                if nnz(NewNeighbors_vec)>=MinPts
                    Diff_NewNeighbors= (NewNeighbors_vec-Neighbors_vec)>0;
                    search_index = cat(1,search_index,find(Diff_NewNeighbors));
                    Neighbors_vec = Neighbors_vec + Diff_NewNeighbors;
                end
            end
            if IDX(jj)==0
                IDX(jj)=C;
            end
            
            kk = kk + 1;
            if kk > length(search_index)
                break;
            end
        end
    end

    function Neighbors_vec=RegionQuery(ii)     
        Neighbors_vec=(D(:,ii)>0 & D(:,ii)<=epsilon );
    end

end




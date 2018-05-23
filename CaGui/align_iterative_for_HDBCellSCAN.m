function ops = align_iterative_for_HDBCellSCAN(data, ops)

% uu = squeeze(sum(sum(data(:,:,:).^2,1),2));
% [~, isort] = sort(uu, 'descend');
% ops.mimg        = data(:,:,isort(50));

fracImgPreAlign = getOr(ops, 'fracImgPreAlign', 1/2);
maxImgPreAlign = round(size(data,3) * fracImgPreAlign);

ops.mimg = pick_reg_init(data);

dsold = zeros(size(data,3), 2);
err = zeros(ops.NiterPrealign, 1);
%%
tempSubPixel = ops.SubPixel;
ops.SubPixel = Inf;
 fprintf('Reg D(Error)=:')
for i = 1:ops.NiterPrealign    
    
    [dsnew, Corr]  = registration_offsets(data, ops, 1);
    dreg  = register_movie(data, ops, dsnew);
    [~, igood] = sort(Corr, 'descend');
    if i<floor(ops.NiterPrealign/2)        
        igood = igood(1:min(maxImgPreAlign,100));  
    else
        igood = igood(1:maxImgPreAlign);  
    end
    ops.mimg = mean(dreg(:,:,igood),3);
    
    err(i) = mean(sum((dsold - dsnew).^2,2)).^.5;
    fprintf('%3.3f, ',err(i)); 
    dsold = dsnew;
   
     if ops.showTargetRegistration
        myfigure('RefImage');cla;
        
        imagesc(ops.mimg)
        colormap('gray')
        title(sprintf('RefImage: Itr%d, mouse %s, date %s', ...
            i, ops.mouse_name, ops.date),'Interpreter','none')
        drawnow
     end
    if err(i)<1
        fprintf('Converged.\n')
        break;
    end
end
ops.SubPixel = tempSubPixel;

ops.AlignNanThresh = median(Corr) - 4*std(Corr);
ops.ErrorInitialAlign = err;
ops.dsprealign = dsnew;
ops.dsprealign_Corr = Corr;
switch ops.RegShape
    case 'Square'
        ops.Ly = size(data,1);
        ops.Lx = size(data,2);
    case 'GrinLens'
        % do not change Ly,Lx
    otherwise
        
end
end 

function F = ResonanceImagingPhaseShifter(F,shift)
[H,~,~]=size(F);

yrange = 2:2:H;


if shift>0
    F(yrange, (1+shift):end,:,:) = F(yrange, 1:(end-shift),:,:);
else
    F(yrange, 1:(end+shift),:,:) = F(yrange, (1-shift):end,:,:);
end

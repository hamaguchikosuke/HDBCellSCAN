function TiffWriter2(image,fname,bitspersamp,writeopt)

if nargin<=3
    writeopt = 'w'; %  open Tiff file for writing; discard existing contents.
%     writeopt = 'a'; % open or create Tiff file for writing; created files will 
                    % be in 32-bit Tiff format; any existing file format will
                    % be preserved; append image data to end of file 
end

t = Tiff(fname,writeopt);
tagstruct.ImageLength = size(image,1);
tagstruct.ImageWidth = size(image,2);
tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
if bitspersamp==8
    tagstruct.BitsPerSample = 8;
end
if bitspersamp==16
    tagstruct.BitsPerSample = 16;
end
if bitspersamp==32
    tagstruct.BitsPerSample = 32;
end
tagstruct.SamplesPerPixel = 1;
tagstruct.RowsPerStrip = size(image,1);
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software = 'TiffWriter2';
t.setTag(tagstruct);
t.write(image(:,:,1));
for i=2:size(image,3)
    t.writeDirectory();
    t.setTag(tagstruct);
    t.write(image(:,:,i));
end
t.close();
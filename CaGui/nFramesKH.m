function n = nFramesKH(tiff)
%nFrames find the number of frames in the Tiff

%keep guessing until we seek too far
guess = 1000;
overSeeked = false;

f = dir(tiff);
if isempty(f)
    error('File not found!');
end
 info=imfinfo(tiff);
if f.bytes>4e9 && length(info)==1
   
    
    % This is a big Tiff, use different strategy to read this file.
    expression = ['\w*images=(?<n>\d+)\n\w*'];
    tokenNames = regexp(info.ImageDescription,expression,'names');
    if ~isempty(tokenNames.n)
        n = str2double(tokenNames.n);
    end
    return;
else
  
    n = length(info);
    return;
end

% %%
% if ischar(tiff)
%   tiff = Tiff(tiff, 'r');
%   closeTiff = onCleanup(@() close(tiff));
% end
% 
% while ~overSeeked
%   try
%     tiff.setDirectory(guess);
%     guess = 2*guess; %double the guess
%   catch ex
%     overSeeked = true; %we tried to seek past the last directory
%   end
% end
% %when overseeking occurs, the current directory/frame will be the last one
% n = tiff.currentDirectory;

end

classdef imcircbuff < handle
% Image Circular buffer class.
% define cb=imcircbuff(height,width,Nch,nframe,type);
% height, width : size of image. 
% Nch           : number of channels.
% nframe        : number of frames to record
% type          : 'uint8','uint16'
% 
% In the circular buffer, height x width x nframe 3D matrix 
% is prepared. The latest data is stored from the end of the data,
% and shifted to earlier index. 
% 
% To construct image circular buffer class, 
% >> imcb = imcircbuff(512,512,1,100,'uint16'); 
% 
% to put data,
%   imcb.putdata(randi(2^16,512));
% 
% to get data 
%   data=imcb.getdata;  
% which returns the whole data in the order where the last added chunk
% comes to the end of data.
%   data=cb.getdata(index); 
% index=  1:nFrame is equal to imcb.getdata;
% index:  nFrame or 0 returns the most recently added data.
% index: nFrame-1 or -1                                                                                                                                                                                                          returns the second latest data.
% index: 1  returns the oldest data added to the circular buffer.
%                           
% data=imcb.getdata(index,roi); returns the image of a specific roi.
% roi is defined as a single vector. 
% 1 is the upper-left corner, 
% height is the bottom-left corner,
% width x height is the bottom-right corner. 
% You can generate this index by using ind2sub function.
% 
% ex) 
% imcb=imcircbuff(512,512,1,100,'uint16');
% for ii=1:100
% imcb.putdata(ii*ones(512));
% end
% 
% ex) 
% out1=imcb.getdata([1:10]);
% out2=imcb.getdata([-10:10],[1:2]);

    properties (SetAccess = private)
        width
        height
        Nch
        nFrame 
        type
        data
        ptr     = 1;
    end
    
    properties (Dependent)
        TotalSize
    end
    
    
    methods
        function obj = imcircbuff(height,width,Nch,nFrame,type)
         obj.height  = height;
         obj.width   = width;
         obj.Nch     = Nch;
         obj.nFrame  = nFrame;
         obj.type    = type;
%          obj.data    = zeros(width,height,nFrame)
         switch type
             case 'uint8'
                 obj.data = uint8(zeros(height,width,Nch,nFrame));
             case 'uint16'
                 obj.data = uint16(zeros(height,width,Nch,nFrame));
             otherwise
                 error('Unknown data type %s!',type)
         end
      end
    end
    
    methods
        function TotalSize = get.TotalSize(obj)
            TotalSize = [obj.height,obj.width,obj.Nch,obj.nFrame];
        end %      
        function obj = set.TotalSize(obj,~)
            fprintf('Total size: [%d,%d,%d,%d]\n',obj.TotalSize)
            error('You cannot set TotalSize explicitly');
        end 
    end
    
    methods 
        function putdata(obj,arg1)
%             size(arg1)
%             size(obj.data)
            obj.data(:,:,:,obj.ptr)=arg1;
            
            obj.ptr=obj.ptr+1;
            if (obj.ptr>obj.nFrame)
                obj.ptr=1;
            end
        end
        
        function out=getdata(obj,varargin)
            pos=(obj.ptr-1);
            out=circshift(obj.data,[0,0,0,-pos]); %  output is circshifted to place the oldest on the first, latest on the end. 
            
            if nargin>=2
                arg1 = varargin{1};
                nega_index=(arg1<=0);
                arg1(nega_index)=arg1(nega_index)+obj.nFrame;
                arg1 = mod(arg1,obj.nFrame);
                arg1(arg1==0)=obj.nFrame;
            end
            
            if nargin==1
%                 out = out; % do nothing;
            elseif nargin==2
                out = out(:,:,:,arg1);
            elseif nargin==3
                out = out(:,:,varargin{2},arg1);
            elseif nargin==4
                out = out(:,varargin{3},varargin{2},arg1);
            elseif nargin==5
                out = out(varargin{4},varargin{3},varargin{2},arg1);
            end
                
            
%             if nargin==1
%                 out=circshift(obj.data,[0,0,0,-pos]);
%             elseif nargin==2
%                 arg1 = arg1+pos;
%                 nega_index=(arg1<=0);
%                 % wrap around the negative index 
%                 arg1 = mod(arg1,obj.nFrame);
%                 arg1(arg1==0)=obj.nFrame;
% %                 arg1(nega_index)=arg1(nega_index)+obj.nFrame;
%                 out=obj.data(:,:,:,arg1);
%             elseif nargin==3
%                 arg1 = arg1+pos;
%                 nega_index=(arg1<=0);
%                  % wrap around the negative index 
%                   arg1 = mod(arg1,obj.nFrame);
%                 arg1(arg1==0)=obj.nFrame;
%                [indx,indy]=ind2sub([obj.height,obj.width],arg2);
%                 out=obj.data(indx,indy,:,arg1);
%             end
        end
        
        function dispdata(obj)
            obj.data
        end
    end
    
end
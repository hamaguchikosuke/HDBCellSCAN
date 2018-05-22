classdef circbuff < handle
% Circular buffer class.
% define cb=circbuff(ChunkSize,N_Chunk,n_channel);
% ChunkSize: length of the data chunk saved each time,
% N_Chunk:   the number of chunk in the circular buffer
% N_Channel: number of channel to be recorded.
% 
% In the circular buffer, N_Channel x (ChunkSize x N_Chunk) 2D matrix 
% is prepared, and N_Channel x ChunkSize 2D matrix is inserted 
% in each putdata(). 
% 
% to put data,
%   cb.putdata(randn(chunksise,n_channel))
% to get data 
%   data=cb.getdata; 
% which returns the whole data in the order where the last added chunk
% comes to the end of data.
%   data=cb.getdata(index); 
% where index = (-chunksize+1):0 returns the most recently added data.
% index = 1:chunksize is the oldest data added to the circular buffer.
% 
% data=cb.getdata(index,channel); returns the specific channels index. 
% 
% ex) 
% cb=circbuff(10,20,2);
% for ii=1:21
% cb.putdata(repmat([[1:10]+(ii-1)*10]',1,2));
% end
% 
% ex) 
% out1=cb.getdata([1:10]);
% out2=cb.getdata([-10:10],[1:2]);

    properties (SetAccess = private)
        ChunkSize
        N_Chunk
        N_Channel 
        data
        ptr     = 1;
    end
    
    properties (Dependent)
        TotalSize
    end
    
    
    methods
        function obj = circbuff(chunksize,n_chunk,n_channel)
         obj.ChunkSize   = chunksize;
         obj.N_Chunk     = n_chunk;
         obj.N_Channel  = n_channel;
         obj.data    = nan(chunksize*n_chunk,n_channel);
      end
    end
    
    methods
        function TotalSize = get.TotalSize(obj)
            TotalSize = [obj.ChunkSize*obj.N_Chunk,obj.N_Channel];
        end %      
        function obj = set.TotalSize(obj,~)
            fprintf('%s[%d,%d]\n','TotalSize is: ',obj.TotalSize)
            error('You cannot set TotalSize explicitly');
        end 
    end
    
    methods 
        function putdata(obj,arg1)
            index=[1:obj.ChunkSize]+(obj.ptr-1)*obj.ChunkSize;
            obj.data(index,:)=arg1;
            obj.ptr=obj.ptr+1;
            if (obj.ptr>obj.N_Chunk)
                obj.ptr=1;
            end
        end
        
        function out=getdata(obj,arg1,arg2)
            pos=(obj.ptr-1)*obj.ChunkSize;
            if nargin==1
                out=circshift(obj.data,[-pos,0]);
            elseif nargin==2
                arg1 = arg1+pos;
                nega_index=(arg1<=0);
                arg1(nega_index)=arg1(nega_index)+obj.TotalSize(1);
                arg1 = mod(arg1,obj.TotalSize(1));
                arg1(arg1==0)=obj.TotalSize(1);
                out=obj.data(arg1,:);
            elseif nargin==3
                arg1 = arg1+pos;
                nega_index=(arg1<=0);
                arg1(nega_index)=arg1(nega_index)+obj.TotalSize(1);
                out=obj.data(arg1,arg2);
            end
        end
        
        function dispdata(obj)
            obj.data
        end
    end
    
end
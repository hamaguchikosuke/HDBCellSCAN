classdef mlttiff < handle
% multiple tiff file class, to treat multiple tiff files as a single series of multipage tiff.
% 
% 
% To construct multiple tiff file class (mlttiff), 
% >> mt=mlttiff(filenames);
% filenames : a cell array of filenames. 
% note that, the width and height of the image size needs to be the same.  
% 
% to get the whole file information,
% mt.height: height of a slice of image (size of 1st dim)
% mt.width: width of a slice of image    (size of 2nd dim)
% mt.depth : depth of a slice of image (size of 3rd dim)
% mt.nSeries =number of  time series (size of 4th dim)
% 
% if the image is a gray scale, 3rd dimension is squeezed and series becomes 3rd dim. 
% 
% to get an image,
%  mt.img(:,:,1);
% by Kosuke Hamaguchi 2015

    properties 
           MaxMemory = 2e9; % By default, limit the internal buffer size as 2GB.
           SIlentMode = 0;      % 1: show progress bar. 0: no show. 
    end
    
    properties (SetAccess = private)
        width
        height
        nSeries
        type
        buffer
        filenames
        CumNSeries
        Loaded
        indices_covered = [];
    end
    
     properties (Dependent)
        TotalNSeries
        MaxSlice
        BytePerPixel
    end
    
    methods
        function obj = mlttiff(f)
        
          obj.Loaded = 0;
      
          if iscell(f)
                 obj.filenames = f;
                 obj.nSeries = zeros(1,length(f));
                 fprintf('Checking tiff pages (mlttiff)...')
                 for ii=1:length(f)
                     obj.nSeries(ii)=nFramesKH(f{ii});
                     fprintf('%d,',obj.nSeries(ii));
                 end
                 fprintf('done\n')
                 tiff = Tiff(f{1}, 'r');
          else
              obj.filenames =cell(1);
              obj.filenames{1}=f;
              obj.nSeries=nFramesKH(f);
               tiff = Tiff(f, 'r');
          end
%           obj.nSeries=cellfun(@nFrames, f);
           
          obj.CumNSeries = [0 cumsum(obj.nSeries)]; % add zero for later purpose.
   
          
          
          obj.width = tiff.getTag('ImageWidth');
          obj.height = tiff.getTag('ImageLength');
          obj.type = class(read(tiff));
          
         switch obj.type
             case { 'uint8',  'uint16'}
                  % do nothing
             otherwise
                 error('Unknown data type %s!',type)
         end
         
        end
    end
    
    methods
        function TotalNSeries = get.TotalNSeries(obj)
            TotalNSeries = obj.CumNSeries(end);
        end %
        
         function BytePerPixel = get.BytePerPixel(obj)
                  switch obj.type
                      case 'uint16'
                          BytePerPixel =2;
                      case 'uint8'
                          BytePerPixel =1;
                      otherwise
                          error('Unknown file type %s',obj.type);
            end
            
        end %
        
        function MaxSlice = get.MaxSlice(obj)
            SingleSliceSize = obj.width * obj.height * obj.BytePerPixel;            
            MaxSlice = floor(obj.MaxMemory/SingleSliceSize);
            MaxSlice = min(obj.TotalNSeries,MaxSlice); 
        end %
    end
    
    methods   
              
        function out=get_single_slice(obj,index)
            
            
            already_read_ind = find(obj.indices_covered==index);
            
            if ~isscalar(index)
                error('get_single_slice(index) must have a scalar index');
            end
            
            if isempty(already_read_ind)
                index = max(index,ceil(obj.MaxSlice/2)); % in case hitting lower bound
                index = min(index,obj.TotalNSeries-ceil(obj.MaxSlice/2)); % in case hitting upper bound
                load_index = [-ceil(obj.MaxSlice/2):floor(obj.MaxSlice/2)]+index;
                load_index = load_index( load_index>0 & load_index <= obj.TotalNSeries);
                load_index = load_index(1:obj.MaxSlice);
                obj.set_buffer(load_index);
                already_read_ind = find(obj.indices_covered==index);
            end
            
            out=obj.buffer(:,:,already_read_ind);
        end
        
        function  out = get_stacks(obj,indices)
            % get stack of data directly from disk.
            Expected_DataSize = obj.width * obj.height * length(indices) * obj.BytePerPixel;
            if Expected_DataSize > obj.MaxMemory
                error('Exceeding max memory size %f',obj.MaxMemory);
            end
        
            out = ones(obj.height,    obj.width, length(indices), obj.type);
             
            StartInd=obj.CumNSeries(1:end-1)+1;
            EndInd = obj.CumNSeries(2:end);
            FilesToRead = [];
            
            for ii=1:length(StartInd)
                if any(indices>=StartInd(ii) & indices <=EndInd(ii))
                    FilesToRead =cat(1, FilesToRead,ii);
                end
            end
            
            ObjIndex =0;
            
            for FileInd=1:length(FilesToRead)
                
                LocalInd = indices;
                % in case local index started more earlier
                LocalInd = LocalInd(LocalInd>=StartInd(FilesToRead(FileInd)));
                % in case local index extend to the next file
                LocalInd = LocalInd(LocalInd<=EndInd(FilesToRead(FileInd)));
                
                % to make it local index
                LocalInd = LocalInd-StartInd(FilesToRead(FileInd))+1;
                ObjIndex = ObjIndex(end)+[1:length(LocalInd)];
                out(:,:,ObjIndex)=BigTiffReader(obj.filenames{FilesToRead(FileInd)} ,LocalInd);
                
            end
        end
      
        function  out = init_buffer(obj)
            % for a given set of indices, check the necessary memory size,
            Expected_DataSize = obj.width * obj.height * obj.MaxSlice* obj.BytePerPixel;
            obj.buffer = ones(obj.height,    obj.width,  obj.MaxSlice, obj.type);

            out = 1;
        end
                
        
        function  out = set_buffer(obj,indices)
%             to set buffer 

            % index inside the buffer.
            % ex) indices = [10:20], then index within the buffer, buf_ind=1:11.
            if nargin==1
                indices = 1:obj.MaxSlice;
            elseif length(indices)>obj.MaxSlice
                error('length of indices (%d) is larger than allocated buffer size (%d)',...
                    length(indices),obj.MaxSlice);
            end
            out  = 1;
            buff_ind = 1:length(indices);
            
            if isempty(obj.indices_covered)
                obj.init_buffer;
            elseif length(obj.indices_covered) > length(indices)
                % obj.buffer will be trimmed at the end of the code.
            elseif length(obj.indices_covered) < length(indices)
                n =   length(indices)- length(obj.indices_covered) ;
                obj.buffer = cat(3,obj.buffer,zeros(obj.height,obj.width,n,obj.type)); % expand the buffer.
            end
            
            % recycle the data that we already read.
            [indices_common,IA,IB]=intersect(obj.indices_covered,indices,'stable');
            obj.buffer(:,:,buff_ind(IB))=obj.buffer(:,:,IA);
            indices_to_read = setdiff(indices,indices_common,'stable');
            
            % now read the data that we yet to read
            StartInd=obj.CumNSeries(1:end-1)+1;
            EndInd = obj.CumNSeries(2:end);
            FilesToRead = [];
            
            for ii=1:length(StartInd)
                if any(indices_to_read>=StartInd(ii) & indices_to_read <=EndInd(ii))
                    FilesToRead =cat(1, FilesToRead,ii);
                end
            end
            
            ObjIndex = length(indices_common);
            
            for FileInd=1:length(FilesToRead)
                
                % limit the index within this file
                WIthinLogical = (indices_to_read>=StartInd(FilesToRead(FileInd))) & ...
                    (indices_to_read<=EndInd(FilesToRead(FileInd)));
                
                % fill in
                LocalInd = indices_to_read(WIthinLogical) - StartInd(FilesToRead(FileInd))+1;
                ObjIndex = ObjIndex(end)+[1:length(LocalInd)];
                obj.buffer(:,:,ObjIndex)=BigTiffReader(obj.filenames{FilesToRead(FileInd)} ,LocalInd);
                
            end
            
            obj.buffer = obj.buffer(:,:,1:ObjIndex(end)); 
            
            obj.indices_covered = indices;
        end
        
    end
end
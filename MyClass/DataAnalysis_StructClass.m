classdef DataAnalysis_StructClass < handle
    % DataAnalysis Struct Class
    %
    % The purpose of this class is to provide an easy way to store
    % data analysis results. Typical flow is
    %
    % 1) initialize structure,
    %   >> obj=DataAnalysis_StructClass();
    %
    % 2) put data. 1st index must be always identifier.
    %   >> obj.put('kh001','data1name',[1:5],'data2name',rand(1),...)
    %
    %
    %   When you put first data into obj, it initialize the data structure.
    %   Importantly, when the first identifier is a string (numeric),
    %     the following identifier must be always a string (numeric).
    %
    %   it puts data in a hidden structure
    %   S(*).identifier= 'kh001';
    %   S(*).data1name    =[1:5];
    %   S(*).data2name    =0.3321;
    %
    % 3) delete data
    %   >> obj.delete('identifier','kh001');
    %   >> obj.delete('kh001'); You can omit 'identifier'.
    %
    % 4) save data
    %   >> obj.save(filename)
    %
    % 5) load data
    %   >>  obj=DataAnalysis_StructClass('load',filename);
    %
    % 6) retrieve data in structure format.
    %   >> data=obj.get('kh001');
    %   >> data=obj.get('identiifer',10011); OR
    %   >> data=obj.getall; % to return all the data.
    %
    % An example of usage of this class.
    % hint) put data in column wise.
    %
    %     data= DataAnalysis_StructClass;
    %     data.put('kh001','meanX',1,'X',[1 1 1]');
    %     data.put('kh002','meanX',2,'X',[1 2 3]');
    %     data.put('kh003','meanX',3,'X',[2 3 4]');
    %     data.put('kh004','meanX',4,'X',[4 5 3]');
    %     mydata=data.getall;
    % %     now you can retrieve meanX of all the file by
    %     [mydata(:).meanX]
    % %     and also, you can retrieve the all X data if they have the same
    % % size.
    %     [mydata(:).X]
    % by KH 20121117
    %
    
    properties
        S
    end
    
    properties (SetAccess = private)
        identifier_class
        identifiers
    end
    
    properties (Dependent)
        
    end
    
    
    methods
        function obj = DataAnalysis_StructClass(opt,val1)
            
            if nargin==0
                opt='init';
            end
            
            switch opt
                case 'init'
                    obj.S=[];
                case 'load'
                    filename=val1;
                    reload(obj,filename);
                otherwise
                    error('Unknown option!');
            end
        end
        
        function put(obj,new_identifier,varargin)
            K=length(obj.S);
            %             fprintf('identifier:%s\n',new_identifier);
            
            
            if K==0 % First time to put data.
                obj.sanity_check_identifier(new_identifier);
                
                obj.identifier_class=class(new_identifier);
                obj.identifiers={new_identifier};
                I=1;
            elseif ~isa(new_identifier,obj.identifier_class)
                error('This identifier has different class from the original (%s)',obj.identifier_class);
            else
                
                I=return_index(obj,new_identifier);
                
                if isempty(I)
                    I=K+1;
                else
                    fprintf('%s found. Overwritten.\n',obj.S(I).identifier);
                    obj.S(I)=[];
                end
            end
            
            
            obj.S(I).identifier = new_identifier;
            for ii=1:(nargin-2)/2
                new_var=varargin{2*ii-1};
                new_val=varargin{2*ii};
                obj.S(I).(new_var)=new_val;
            end
            obj.identifiers={obj.S(:).identifier};
        end
        
        
        
        function out = get(obj,key_identifier)
            
            I=return_index(obj,key_identifier);
            
            if isempty(I)
                out=[];
            else
                out=obj.S(I);
            end
        end
        
        function SUCCESS = delete(obj,key_identifier,arg1)
            
            if nargin==3
                if strmatch(key_identifier,'identifier','exact')
                    key_identifier=arg1;
                else
                    error('cannot use other data as a key');
                end
            end
            I=return_index(obj,key_identifier);
            
            if isempty(I)
                fprintf('not found.\n')
                SUCCESS=0;
            else
                SUCCESS=1;
                fprintf('deleted.\n')
                newind=setdiff([1:length(obj.S)],I);
                obj.S=obj.S(newind);
                obj.identifiers=obj.identifiers(newind);
                
            end
        end
        
        function out = getall(obj)
            out=obj.S;
        end
        
        function I=return_index(obj,key)
            switch obj.identifier_class
                case 'char'
                    I=strmatch(key, obj.identifiers,'exact');
                otherwise
                    I=find(obj.identifiers==key);
            end
        end
        
        function save(obj,filename)
            S=obj.S;
            identifier_class=obj.identifier_class;
            identifiers     =obj.identifiers;
            fprintf('saving DataAnalysis_StructClass component into %s...\n',filename);
            save(filename,'S','identifier_class','identifiers');
        end
        
        function reload(obj,filename)
            load(filename,'S','identifier_class','identifiers');
            obj.S=S;
            obj.identifier_class=identifier_class;
            obj.identifiers     =identifiers;
        end
        
        
        function sanity_check_identifier(obj,new_identifier)
%             class(new_identifier)
            if ~ischar(new_identifier) && ~isnumeric(new_identifier)
                error('Identifier cannot be a cell. It either char or numeric');
            end
        end
    end
    
    
    
    
end
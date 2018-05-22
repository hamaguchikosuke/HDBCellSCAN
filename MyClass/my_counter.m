classdef my_counter < handle
% my_counter class.
% define ctr=my_counter();
% It has three states, stop, start, pause. 
% initial state is stop, and then once it started, it starts to count until
% paused or stopped.
% 
% >> ctr.stop; % stops counter and reset to 0.
% >> ctr.pause; % pause counter, keep counter value.
% >> ctr.start; % start counter. 
% >> ctr.update; % increment counter when ctr.state is 'start'.
% otherwise, do nothing. you can increment n by 
% >> ctr.update(n);
% 
% >> ctr.set(n); %you can set counter at any value, at any state.
% 
% >> ctr.parameter = 10; % you can set a parameter for the counter.
%    % this parameter is not updated and fixed. 
% ex)
% ctr=my_counter();
% ctr.start
% for ii=1:10
%   ctr.update;
% end
% ctr
% ctr.set(1000)
% ctr

    properties (SetAccess = private)
        state
        counter
        running
        stopping
        pausing
    end
    
    properties (Dependent)
        
    end
    
    properties
          parameter
    end
    
    methods
        function obj = my_counter()
         obj.state          = 'stop';
         obj.counter        = 0;
         obj.parameter      = [];
         obj.running        = 0;
         obj.stopping       = 1;
         obj.pausing        = 0;
      end
    end
    
   
    
    methods 
        function stop(obj)
%             fprintf('Counter stopped and reset to 0.\n',obj.counter)
            obj.state = 'stop';
            obj.running        = 0;
            obj.stopping       = 1;
            obj.pausing        = 0;
            obj.counter = 0;
        end
        
        function start(obj,arg1)
%             fprintf('Counter start from %d.\n',obj.counter)
            obj.state = 'start';
            obj.running        = 1;
            obj.stopping       = 0;
            obj.pausing        = 0;
%             obj.counter = 0;
        end
        
        function pause(obj)
            %              fprintf('Counter paused at %d.\n',obj.counter)
            obj.state = 'pause';
            obj.running        = 0;
            obj.stopping       = 0;
            obj.pausing        = 1;
            %             obj.counter = 0;
        end
        
        function out=update(obj,arg1)
            
            switch obj.state
                case 'start'
                    if nargin==1
                        obj.counter = obj.counter+1;
                        
                    elseif nargin==2
                        obj.counter = obj.counter+arg1;
                    end
                case {'stop','pause'}
                    % do nothng;
                otherwise
                    error('No such a state %s',obj.state);
            end
            
        end
        
         function out=set(obj,arg1)
                       
             if nargin==2
                 obj.counter = arg1;
             else
                 error('Please assign set value obj.set(n).')
             end
  
        end
        
        function disp(obj)
            fprintf('Counter:%d\n',obj.counter)
        end
    end
    
end
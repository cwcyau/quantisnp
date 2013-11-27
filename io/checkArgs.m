function [trueOrFalse, arg] = checkArgs(arg, argType, argBounds, varargin)

nStrings = nargin - 3;

switch argType

    case 'numeric'

        argFail = 0;
        
        if strmatch( class(arg), 'char', 'exact' ) % if the input is a string
            
            loc = regexp(arg, '\d');	% find the number of digits
            
            if length(loc) == length(arg)	% if the number of digits equals the length of the string
                arg = str2num(arg); % convert to a numeral
            else
                argFail = 1;
            end
            
        end

        % if string was numeric, check that the value is greater than 0
        if argFail == 0
            
            if arg >= argBounds(1) & arg <= argBounds(2)
                
            else
                argFail = 1;
                
            end
            
        end

        % string contains letters, or is less than 0 exit with an error
        if argFail == 1
            
            for i = 1 : nStrings
                
                disp(varargin{i});
                
            end
            
            trueOrFalse = 1;
            
            return;
            
        end

end

trueOrFalse = 0;
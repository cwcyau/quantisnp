function dispstr(varargin)

nArgs = nargin;

commandStr = ['sprintf(''' varargin{1}];

for j = 2 : nArgs
	if strmatch('char', class(varargin{j}))
		commandStr = [ commandStr ''', ''' varargin{j} ];
	else
		commandStr = [ commandStr ''', ''' num2str(varargin{j}) ];
	end
end
commandStr = [ commandStr ''' );'];

str = eval(commandStr);

disp(str);
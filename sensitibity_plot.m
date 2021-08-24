function [mag] = sensitibity_plot( G, freq, vars, varvals, dt)
% SENSITIVITY_PLOT  Display a plot of parameter sensitivity
% [mag] = param_sensitibity( G, freq, vars, varvals, dt)

arguments
    G sym
    freq (1,:) {mustBeNumeric}
    vars (1,:) cell
    varvals (1,:)  cell {mustBeEqualSize(vars,varvals)}
    dt (1,1) {mustBePositive} = 1
end

syms omega s z

temp=subs(G,vars,varvals);
varlist=symvar(temp);

if numel(varlist)>1
    errmsg=sprintf('(%s)',varlist);
    error(strcat('incomplete variable list:',errmsg));
end

if strcmp(char(varlist),'s')
    q='s';
    qval=1i*omega;
    dt=eps;
else
    q='z';
    qval=exp(1i*omega*dt);
end

numvar=numel(vars);
mag=ones(numel(freq),numvar)*nan;
for k=1:numvar
    var=vars{k};
    dGdth=diff(G,var);
    ps=subs(dGdth*var,vars,varvals);
    psbode=symfun(subs(20*log10(abs(ps)),q,qval),omega);
    for l=1:find(freq<(1/dt)/2,true,'last')
        try
            mag(l,k)=double(psbode(2*pi*freq(l)));
        catch
            % division by zero error is possible
            mag(l,k)=nan;
        end
    end
end
lname=cellfun(@char,vars,'UniformOutput',false);
zx=semilogx(freq,mag);
for k=1:numvar
    zx(k).DisplayName = lname{k};
end


Gs=subs(G,vars,varvals);
Gsbode=symfun(subs(20*log10(abs(Gs)),q,qval),omega);
mag=double(Gsbode(2*pi*freq));
lname{end+1} = 'nominal';
hold on
semilogx(freq,mag,'k--');

legend(lname)
xlim(freq([1,end]))
xlabel('Frequency [Hz]')
ylabel('Parameter Sensitivity [dB]')

hold off
end

function mustBeEqualSize(a,b)
% Test for equal size
if ~isequal(size(a),size(b))
    eid = 'Size:notEqual';
    msg = 'Size of fourth input must equal size of third input.';
    throwAsCaller(MException(eid,msg))
end
end


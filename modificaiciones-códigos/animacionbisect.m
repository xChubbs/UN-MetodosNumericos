function c = animacionbisect(f,a,b,delta)
% Animacion del comportamiento del metodo de biseccion
% Entradas:
% f = funcion de iteracion de punto fijo
% a = extremo inferior (izquierdo) del intervalo de busqueda
% b = extremo superior (derecho) del intervalo de busqueda
% delta = error admisible, tolerancia
% Salida:
% c = aprocimacion del cero de la funcion
%
% Desarrollado en la Universidad Nacional de Colombia Sede Medellin 2007
% Basado en:
%  METODOS NUMERICOS: Programas en Matlab
% (c) 2004 por John H. Mathews y Kurtis D. Fink
%  Software complementario acompa�ando al texto:
%  METODOS NUMERICOS con Matlab, Cuarta Edicion
%  ISBN: 0-13-065248-2
%  Prentice-Hall Pub. Inc.
%  One Lake Street
%  Upper Saddle River, NJ 07458

ya = feval(f, a);
yb = feval(f, b);

if  ya*yb > 0

    hh = errordlg(['La funcion tiene el mismo...' ...
        'signo en a y en b f(a)*f(b)>0'],'Error');
    waitfor(hh)

    return
end

max1 = 1 + round((log(b-a) - log(delta)) / log(2));

close all

fplot(f,[0.1*round(10*a-1) 0.1*round(10*b+1)],'Color','g');
grid on
hold on

fplot(@(x)0*ones(size(x)),[0.1*round(10*a-1) 0.1*round(10*b+1)],'Color','b');

hh = warndlg({'Oprima OK para empezar la animacion',...
    'No cierre la ventana hasta que finalice la animacion'},'Instruccion');

waitfor(hh)
k=0;
c = (a + b) / 2;

yc = feval(f, c);

ap=plot([a a],[0 feval(f,a)],'--ks');
bp=plot([b b],[0 feval(f,b)],'--kd');
cp=plot([c c],[0 feval(f,c)],'--ro');

legend({'funcion g(x)','recta y=0','a_k','b_k','c_k'},'Location','Best')

for  k = 1:max1

    c = (a + b) / 2;
    yc = feval(f, c);

    pause(0.3)
    set(ap,'XData',[a a],'YData',[0 feval(f,a)])
    set(bp,'XData',[b b],'YData',[0 feval(f,b)])
    pause(0.5)
    set(cp,'XData',[c c],'YData',[0 feval(f,c)])
    pause(0.3)

    if  yc == 0

        a = c;
        b = c;
    elseif  yb*yc > 0

        b = c;
        yb = yc;
    else

        a = c;
        ya = yc;
    end
    xlabel(['aproximacion al cero de f(x), c(' ...
                        mat2str(k) ') = ' mat2str(c)])

    if  b-a < delta, break, end
end

c = (a + b) / 2;

pause(0.3)

set(ap,'XData',[a a],'YData',[0 feval(f,a)])
set(bp,'XData',[b b],'YData',[0 feval(f,b)])

pause(0.3)

set(cp,'XData',[c c],'YData',[0 feval(f,c)])

err = abs(b - a)
yc = feval(f, c)
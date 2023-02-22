function  [c, err, yc] = bisectn (f, a, b, delta, Iter)

% Entrada - f es la funcion introducida con @
%	      - a y b son los extremos izquierdo y derecho
%	      - delta es la tolerancia
%         - Iter Iteración
% Salida  - c es el cero
%	      - yc = f(c)
% 	      - err es el error estimado para  c

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
    return
end

max1 = 1 + round((log(b-a) - log(delta)) / log(2));

for  k = 1:max1
    
    k
	c = (a + b) / 2

	yc = feval(f, c);

    error = b - a
    
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

    if k == Iter, break, end

	if  b-a < delta, break, end

end

c = (a + b) / 2;
err = abs(b - a);
yc = feval(f, c);

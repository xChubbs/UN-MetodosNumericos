"""
 ******************************************************************************
 * @file           : rutinas_MN.py
 * @author         : alujan
 * @brief          : Header - Rutinas Métodos Numéricos
 *
 ******************************************************************************
 * @authors
 *
 * Universidad Nacional de Colombia
 * Sede Medellín
 * 2023-01S
 *
 ******************************************************************************
"""
# Importación de librarias externas
import numpy.linalg as al
import numpy as np


# ------------------------ Rutinas Univariable ------------------------------ #

def bisect(f, a, b, delta):
  """
  Entrada - f es la funcion introducida con @
          - a y b son los extremos izquierdo y derecho
          - delta es la tolerancia
  Salida  - c es el cero
          - yc = f(c)
          - err es el error estimado para  c
  """
  
  ya = f(a)
  yb = f(b)

  if ya*yb >0: 
    print("La función en el intervalo" +
              "no tiene cambio de signo")

  maxIter = 1 + round((np.log(b-a) -
              np.log(delta)) / np.log(2))

  for k in range(1, maxIter +1):
    
    c = (a + b) / 2
    yc = f(c)

    if yc == 0:
      a = c
      b = c

    elif yb * yc > 0:
      b = c
      yb = yc

    else:
      a = c
      ya = yc

    if b - a < delta:
      break

  c = (a + b)/2
  yc = f(c)
  error = abs(b - a)

  return c, yc, error, k

def fixpt(f, p0 : float, delta : float, maxIter : int):

  for iter in range(2, maxIter +1):
    p = f(p0)
    error = abs(p0 - p)
    errorRel = error / abs(p)
    
    p0 = p

    if (error < delta) or (errorRel < delta): break
    if iter == maxIter: print("Máximo de iteraciones alcanzado")

  return p, error, errorRel


def newtonMod(f, df, p0, delta, epsilon, max1):

    """
    Entrada - f funcion creada con @
            - df funcion derivada creada con @
            - d2f función segunda derivada creada con @
            - p0 es la aproximacion inicial a cero de  f
            - delta es la tolerancia para  p0
            - epsilon es la tolerancia para los valores de la funcion  y
            - max1 es el numero maximo de iteraciones
    Salida  - p0 es la aproximacion de Newton-Raphson hacia cero
            - err es el error estimado para  p0
            - k es el numero de iteraciones
            - y es el valor de la funcion  f(p0)
    """

    for k in range(1, max1 +1):

        # --- Tú código aquí! --- #
        """
        Recuerda que podemos realizar un cambio de variable el cual nos
        permite reducir la multiplicidad del cero sin necesidad de
        conocer la multiplicidad de la raiz.
        ¡¡Debes reescribir la ecuación de avance para este método!!
        """

        err = abs(p1 - p0)
        relerr = 2 * err /(abs(p1) + delta)

        p0 = p1
        y = mu(p0)

        if (err < delta) or (relerr < delta) or (abs(y) < epsilon):
            break

    return p0, err, k, y


# ------------------------ Rutinas Multivariable --------------------------- #

def jacobi(A, b, X0, delta, maxI):

    # Entrada  - A es una matriz no singular  N x N
    #          - B es una matriz  N x 1
    #          - P es una matriz  N x 1; los supuestos iniciales
    #	       - delta es la tolerancia para  P
    #	       - max1 es el numero maximo de iteraciones
    # Salida   - X es una matriz  N x 1: la aproximacion de jacobi a
    #	         la solucion de  AX = B

    #  METODOS NUMERICOS: Programas en Matlab
    # (c) 2004 por John H. Mathews y Kurtis D. Fink
    #  Software complementario acompañando al texto:
    #  METODOS NUMERICOS con Matlab, Cuarta Edicion
    #  ISBN: 0-13-065248-2
    #  Prentice-Hall Pub. Inc.
    #  One Lake Street
    #  Upper Saddle River, NJ 07458

    D = np.diag(np.diag(A))
    L = -np.tril(A, -1)
    U = -np.triu(A, 1)

    TJ = al.inv(D) @ (L + U)
    CJ = al.inv(D) @ b

    for k in range(1, maxI +1):

        X = TJ @ X0 + CJ

        err = abs(al.norm(X - X0))
        rerr = err / (al.norm(X) + delta)

        X0 = X

        if err < delta or rerr < delta: break


    return X, err, rerr, k

def sor(A, b, X0, w, delta, maxI):

    # Entrada  - A es una matriz no singular  N x N
    #          - B es una matriz  N x 1
    #          - P es una matriz  N x 1; el supuesto inicial
    #          - w parametro de sobrerelajacion (0<w<2)
    #	       - delta es la tolerancia para  P
    #	       - max1 es el numero maximo de iteraciones
    # Salida   - Y es una matriz  N x 1: la aproximacion de SOR a
    #	         la solucion de  AX = B

    #  METODOS NUMERICOS: Programas en Matlab
    # (c) 2004 por John H. Mathews y Kurtis D. Fink
    #  Software complementario acompañando al texto:
    #  METODOS NUMERICOS con Matlab, Cuarta Edicion
    #  ISBN: 0-13-065248-2
    #  Prentice-Hall Pub. Inc.
    #  One Lake Street
    #  Upper Saddle River, NJ 07458

    D = np.diag(np.diag(A))
    L = -np.tril(A, -1)
    U = -np.triu(A, 1)

    # --- Tú código aquí! --- #

    for k in range(1, maxI +1):

        # --- Tú código aquí! --- #

        err = abs(al.norm(X - X0))
        rerr = err / (la.norm(X) + delta)

        X0 = X

        if err < delta or rerr < delta: break


    return X, err, rerr, k

def newdim(F, JF, P, delta, epsilon, max1):

    """
    Entrada  - F funcion del sistema creada con @
               JF matriz jacobiana, funcion creada con @
             - P es la aproximacion inicial a la solucion
             - delta es la tolerancia para  P
             - epsilon es la tolerancia para  F(P)
             - max1 es el numero maximo de iteraciones
    Salida   - P es la aproximacion a la solucion
             - iter es el numero de iteraciones realizadas
             - err es el error estimado para  P

    METODOS NUMERICOS: Programas en Matlab
    (c) 2004 por John H. Mathews y Kurtis D. Fink
    Software complementario acompa�ando al texto:
    METODOS NUMERICOS con Matlab, Cuarta Edicion
    ISBN: 0-13-065248-2
    Prentice-Hall Pub. Inc.
    One Lake Street
    Upper Saddle River, NJ 07458
    """

    Y = F(P)

    for iter in range(1, max1 +1):

        J = JF(P)

        Q = P - (al.inv(J) @ Y.T).T

        Z = F(Q)

        err = al.norm(Q - P)
        relerr = err / (al.norm(Q) + delta)

        P = Q
        Y = Z

        if (err < delta)  or (relerr < delta) or (abs(Y).any() < epsilon):
            break
    
    return P, iter, err


# --------------------- Rutinas Aprox. Polinomial --------------------------- #

def lagran (X, Y):

# Entrada  - X es un vector que contiene una lista de las abscisas
#          - Y es un vector que contiene una lista de las ordenadas
# Salida   - C vector que contiene los coeficientes del polinomio
#          - L es una matriz que contiene los coeficientes de los
#            polinomios de Lagrange

#   METODOS NUMERICOS: Programas en Matlab
#   (c) 2004 por John H. Mathews y Kurtis D. Fink
#   Software complementario acompañando al texto:
#   METODOS NUMERICOS con Matlab, Cuarta Edicion
#   ISBN: 0-13-065248-2
#   Prentice-Hall Pub. Inc.
#   One Lake Street
#   Upper Saddle River, NJ 07458

    w = len(X)
    n = w - 1
    L = np.zeros((w, w))

    for k in range(n +1):

        V = 1

        for j in range(n +1):

            if k != j:

                temp = np.poly([X[j]])
                V = np.convolve(V, temp) / (X[k] - X[j])

        L[k, :] = V

    C = Y @ L
   
    return C, L


# --------------------- Rutinas Aprox. EDO --------------------------- #

def euler (f, a, b, ya, M):

# Entrada  - f funcion creada con @
#          - a y b son los extremos izquierdo y derecho
#          - ya es la condicion inicial  y(a)
#          - M es el numero de pasos
# Salida   - E = [T', Y'] donde T es el vector de abscisas y
#            Y es el vector de ordenadas

#  METODOS NUMERICOS: Programas en Matlab
# (c) 2004 por John H. Mathews y Kurtis D. Fink
#  Software complementario acompa�ando al texto:
#  METODOS NUMERICOS con Matlab, Cuarta Edicion
#  ISBN: 0-13-065248-2
#  Prentice-Hall Pub. Inc.
#  One Lake Street
#  Upper Saddle River, NJ 07458

    M = int(M)
    h = (b - a) / M
    T = np.arange(a, b +1, h)
    T = T.tolist()

    Y = []
    Y += [ya]

    for  j in range(0, M):
        Y  += [Y[j] + h * f(T[j], Y[j])]

    return T, Y

def eulers (f, a, b, Za, M):

# Entrada  - f funcion creada con @
#          - a y b son los extremos izquierdo y derecho
#          - ya es la condicion inicial  y(a)
#          - M es el numero de pasos
# Salida   - E = [T', Y'] donde T es el vector de abscisas y
#            Y es el vector de ordenadas

#  METODOS NUMERICOS: Programas en Matlab
# (c) 2004 por John H. Mathews y Kurtis D. Fink
#  Software complementario acompa�ando al texto:
#  METODOS NUMERICOS con Matlab, Cuarta Edicion
#  ISBN: 0-13-065248-2
#  Prentice-Hall Pub. Inc.
#  One Lake Street
#  Upper Saddle River, NJ 07458

    M = int(M)
    h = (b - a) / M
    T = np.arange(a, b +1, h)

    Z = np.zeros( (M+1, len(Za)))
    Z[0, :] = Za

    for j in range(0, M):
        Z[j +1, :] = Z[j, :] + h * f(T[j], Z[j,:])

    return T, Z
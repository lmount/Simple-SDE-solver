# Description
Kloeden - Numerical Solution of stochastic differential
            equations (Springer 1992)  page XXX.
            Strong order 1.0 Runge Kutta scheme.
            http://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_method_%28SDE%29
            
            dX = a(X,t)*dt + b(X, t)*dW
            
##    Syntax:
    Solve_SDE(alfa=None, beta=None, X0=None, dt=None, N=100, t0=0, DW=None) 
        
        
##    Parameters:
        alfa  : a lambda function with two arguments, the X state and the time
                defines the differential equation.
        beta  : a lambda function with two arguments, the X state and the time
                defines the stochastic part of the SDE.
        X0    : The initial state of the SDE. Mandatory for SDE with variables > 1
                (default: gaussian random)
        dt    : The timestep of the solution
                (default: 1)
        N     : The number of timesteps (defines the length of the timeseries)
                (default: 100)
        t0    : The initial time of the solution 
                (default: 0)
        DW    : The Wiener function in lambda notation
                (default: gaussian random number generator, [lambda Y, dt: randn(len(X0)) * sqrt(dt)] )
        
        
##    Examples:
    
###  Simple Wiener Process:
    dX = 0 + 1*dW
    
    
        alfa = lambda X,t: 0
        beta = lambda X,t: 1
        Y = Solve_SDE(alfa=alfa, beta=beta, dt=1, N=1000)
    
        
        
### Stochastic Lorenz Equation:
    dX = s (Y - X) + Y * dW1
    dY = (r X - Y - X*Z) + dW2
    dZ = (X*Y - b Z)  + dW3
    
    
    xL = lambda X, t: 10.0 * (X[1] - X[0])  ;
    yL = lambda X, t: 28.0 * X[0] - X[1] - X[0] * X[2] ;
    zL = lambda X, t: X[0] * X[1] - 8.0/3.0 * X[2] ;

    alfa = lambda X, t: array( [xL(X,t), yL(X,t), zL(X,t)] ); 
    beta = lambda X, t: array( [     X[1],      1,      1] ); 
    X0 = [3.4, -1.3, 28.3];
    Y = Solve_SDE(alfa=alfa, beta=beta, X0=X0, dt=0.01, N=10000)


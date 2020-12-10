# HeatDistributionCalculation
My implementation of iterative heat distribution calculation algorithm on _C_ with the support of multiprocessing by using MPI

## Running
```
mpicc heat_equation.c -lm -o  heat_equation
mpirun -np 4 ./heat_equation
```


## Proof of correctness
Programm result for `T = 0.1, k = 1, h = 0.02, dt = 0.0002`

```
Calculation time: 0.001698s
At x =   0:     u(x, y) = 0.030889, u_accurate = 0.000000
At x = 0.1:     u(x, y) = 0.181060, u_accurate = 0.146691
At x = 0.2:     u(x, y) = 0.315403, u_accurate = 0.278987
At x = 0.3:     u(x, y) = 0.420015, u_accurate = 0.383934
At x = 0.4:     u(x, y) = 0.488039, u_accurate = 0.451286
At x = 0.5:     u(x, y) = 0.510215, u_accurate = 0.474487
At x = 0.6:     u(x, y) = 0.488039, u_accurate = 0.451286
At x = 0.7:     u(x, y) = 0.420015, u_accurate = 0.383934
At x = 0.8:     u(x, y) = 0.315403, u_accurate = 0.278987
At x = 0.9:     u(x, y) = 0.181060, u_accurate = 0.146691
At x =   1:     u(x, y) = 0.030889, u_accurate = 0.000000

Absolute deviation = 0.384742
Mean absolute deviation = 0.034977
```

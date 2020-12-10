#include "math.h"
#include "stdio.h"
#include "stdlib.h"

#include "mpi.h"

double total_length = 1;
double start_temp = 1;
double outside_temp = 0;

int N_POINTS;
double h;
int total_points;

double  T = 0.1;
double dt = 0.0002;
double k = 1;
double start_time;

double accurate_temperature(double x, double t) {
    double total_sum = 0;
    for (int i = 0; i < 10000000; ++i) {
        total_sum += exp(-k * pow(M_PI, 2) * pow(2 * i + 1, 2) * t / pow(total_length, 2))
            * sin(M_PI * (2 * i + 1) * x / total_length) / (2 * i + 1);
    }
    return total_sum * 4 * start_temp / M_PI;
}

void pointer_swap(double** a, double **b) {
    double* temp = *a;
    *a = *b;
    *b = temp;
}

int main(int argc, char* argv[]) {

    // Default initialization for simple test
    N_POINTS = 51;
    h = total_length / (N_POINTS - 1);
    total_points = N_POINTS + 2;


    // This block is for heavy testing
    // Uncomment if you want test it yourself
    //
    // {
    //     N_POINTS = 2000 + 1;
    //     h = total_length / (N_POINTS - 1);
    //     total_points = N_POINTS + 2;

    //     T = 0.0003;
    //     dt = h * h / (2 * k);
    // }

    int rank, size;
    int next, prev, tag_l = 1, tag_r = 2;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) 
        start_time = MPI_Wtime();

    int l = (total_points / size) * rank;
    int r = (total_points / size) * (rank + 1) - 1;
    if (rank + 1 == size)
        r = total_points - 1;
    
    double send_l, send_r;
    MPI_Status stats;

    // Initialization of temperature array and its buffer
    double* temperatures = (double*) calloc(r - l + 1, sizeof(double));
    double* new_temperatures = (double*) calloc(r - l + 1, sizeof(double));
    if (rank == 0) {
        temperatures[0] = outside_temp;
        for (int i = 1; i <= r; ++i)
            temperatures[i] = start_temp;
    } else if (rank + 1 == size) {
        temperatures[r - l] = outside_temp;
        for (int i = 0; i < r - l; ++i) 
            temperatures[i] = start_temp;
    } else {
        for (int i = 0; i <= r - l; ++i)
            temperatures[i] = start_temp;
    }

    double l_add = start_temp;
    double r_add = start_temp;

    // Main loop
    for (double current_temp = 0; current_temp < T; current_temp += dt) {
        if (rank == 0) {
            new_temperatures[0] = outside_temp;
            for (int i = 1; i < r; ++i)
                new_temperatures[i] = temperatures[i] + k * dt / (h * h) * 
                    (temperatures[i - 1] - 2 * temperatures[i] + temperatures[i + 1]);
            new_temperatures[r] = temperatures[r] + k * dt / (h * h) * 
                (temperatures[r - 1] - 2 * temperatures[r] + r_add);
            
            prev = size - 1;
            next = rank + 1;
        } else if (rank + 1 == size) {
            new_temperatures[r - l] = outside_temp;
            for (int i = 1; i < r - l; ++i) 
                new_temperatures[i] = temperatures[i] + k * dt / (h * h) * 
                    (temperatures[i - 1] - 2 * temperatures[i] + temperatures[i + 1]);
            new_temperatures[0] = temperatures[0] + k * dt / (h * h) * 
                    (l_add - 2 * temperatures[0] + temperatures[1]);
            prev = rank - 1;
            next = 0;
        } else {
            for (int i = 1; i < r - l; ++i)
                new_temperatures[i] = temperatures[i] + k * dt / (h * h) * 
                    (temperatures[i - 1] - 2 * temperatures[i] + temperatures[i + 1]);
            
            new_temperatures[0] = temperatures[0] + k * dt / (h * h) * 
                    (l_add -2 * temperatures[0] + temperatures[1]);

            new_temperatures[r - l] = temperatures[r - l] + k * dt / (h * h) * 
                (temperatures[r - l - 1] - 2 * temperatures[r - l] + r_add);

            prev = rank - 1;
            next = rank + 1;
        }
        pointer_swap(&new_temperatures, &temperatures);

        send_l = temperatures[0];
        send_r = temperatures[r - l];

        MPI_Send(&send_l, 1, MPI_DOUBLE, prev, tag_l, MPI_COMM_WORLD);
        MPI_Send(&send_r, 1, MPI_DOUBLE, next, tag_r, MPI_COMM_WORLD);

        MPI_Recv(&l_add, 1, MPI_DOUBLE, prev, tag_r, MPI_COMM_WORLD, &stats);
        MPI_Recv(&r_add, 1, MPI_DOUBLE, next, tag_l, MPI_COMM_WORLD, &stats);
    }

    if (rank == 0) {
        printf("Calculation time: %lfs\n", MPI_Wtime() - start_time);
    }

    // Printing info about some points and calculating errors
    double total_error = 0;
    for (int i = l; i <= r; ++i) {
        double x = (i - 1) * h;
        if ((int)(10 * x) == 10 * x) {
            double accurate_u = accurate_temperature(x, T);
            printf("At x = %*.3g:\tu(x, y) = %.6f, u_accurate = %.6f\n",
                3, x,
                temperatures[i - l], accurate_u);
            total_error += fabs(temperatures[i - l] - accurate_u);
        }
    }

    // Aggregating errors and printing total statistics
    if (rank == 0) {
        for (int i = 1; i < size; ++i) {
            double s;
            MPI_Recv(&s, 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &stats);
            total_error += s;
        }

        printf("Absolute deviation = %lf\n", total_error);
        printf("Mean absolute deviation = %lf\n", total_error / 11);
    } else {
        MPI_Send(&total_error, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
    }

    MPI_Finalize();
}
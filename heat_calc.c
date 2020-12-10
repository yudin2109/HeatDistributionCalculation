#include "stdio.h"
#include "math.h"
#include "mpi.h"
#include "stdlib.h"

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

    N_POINTS = 51;
    h = total_length / (N_POINTS - 1);
    total_points = N_POINTS + 2;

    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) 
        start_time = MPI_Wtime();

    int l = (total_points / size) * rank;
    int r = (total_points / size) * (rank + 1) - 1;
    if (rank + 1 == size)
        r = total_points - 1;
    

    int next, prev, tag_l = 1, tag_r = 2;
    double send_l, send_r;
    MPI_Request reqs[4];
    MPI_Status stats[4];


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
    for (double current_temp = 0; current_temp < T; current_temp += dt) {
        // printf("Process %d got: %lf %lf\n", rank, l_add, r_add);
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

        // MPI_Irecv(&r_add, 1, MPI_DOUBLE, prev, tag_r, MPI_COMM_WORLD, &reqs[0]);
        // MPI_Irecv(&l_add, 1, MPI_DOUBLE, next, tag_l, MPI_COMM_WORLD, &reqs[1]);

        // MPI_Isend(&send_l, 1, MPI_DOUBLE, prev, tag_l, MPI_COMM_WORLD, &reqs[2]);
        // MPI_Isend(&send_r, 1, MPI_DOUBLE, next, tag_r, MPI_COMM_WORLD, &reqs[3]);

        // MPI_Waitall(4, reqs, stats);

        MPI_Send(&send_l, 1, MPI_DOUBLE, prev, tag_l, MPI_COMM_WORLD);
        MPI_Send(&send_r, 1, MPI_DOUBLE, next, tag_r, MPI_COMM_WORLD);

        // printf("Process %d send %lf %lf\n", rank, send_l, send_r);

        MPI_Recv(&l_add, 1, MPI_DOUBLE, prev, tag_r, MPI_COMM_WORLD, &stats[0]);
        MPI_Recv(&r_add, 1, MPI_DOUBLE, next, tag_l, MPI_COMM_WORLD, &stats[1]);
    }

    if (rank == 0) {
        printf("Calculation time: %lf\n", MPI_Wtime() - start_time);
    }

    for (int i = l; i <= r; ++i) {
        double x = (i - 1) * h;
        if ((int)(10 * x) == 10 * x) {
            printf("At x = %*.3g: \nu(x, y) = %.6f, u_accurate = %.6f\n",
                3, x,
                temperatures[i - l], accurate_temperature(x, T));
            // total_error += fabs(temperatures[i] - accurate_temperature(x, T));
        }
    }

    // auto* temperatures = new double[n_points + 2];
    // temperatures[0] = temperatures[n_points + 1] = outside_temp;
    // for (size_t i = 1; i != n_points + 1; ++i)
    //     temperatures[i] = start_temp;
    
    // auto* new_temperatures = new double[n_points + 2];
    // new_temperatures[0] = new_temperatures[n_points + 1] = outside_temp;
    
    // for (double current_time = 0; current_time < T; current_time += dt) {
    //     for (size_t i = 1; i != n_points + 1; ++i) {
    //         new_temperatures[i] = temperatures[i] + k * dt / (h * h) * 
    //             (temperatures[i - 1] - 2 * temperatures[i] + temperatures[i + 1]);
    //         if (isnanf(-new_temperatures[i])) {
    //             std::cout << current_time << " " << i << std::endl;
    //             std::printf("%.5f %.5f %.5f\n", temperatures[i-1], temperatures[i], temperatures[i+1]);
    //             exit(0);
    //         }
    //     }
    //     std::swap(temperatures, new_temperatures);
    // }
    
    // double x = -h;
    // double total_error = 0;
    // int share = (n_points) / 10;
    // for (int i = 0; i < n_points + 2; ++i) {
    //     if (i % share == 1)  {
    //         std::printf("At x = %*.3g: ", 3, x);
    //         std::printf("u(x, y) = %.6f, u_accurate = %.6f\n",
    //             temperatures[i], accurate_temperature(x, T));
    //         total_error += fabs(temperatures[i] - accurate_temperature(x, T));
    //     }
    //     x += h;
    // }
    // std::cout << std::endl;
    // std::printf("TOTAL ERROR = %.4g\n", total_error / 11);

    MPI_Finalize();
}
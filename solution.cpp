#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <chrono>
#include <immintrin.h>

#define NREPET 1

void printUsage(int argc, const char **argv)
{
    printf("Usage: %s N\n", argv[0]);
    printf("Example: %s 1024\n", argv[0]);
}

void verify(const float *A, const float *B, const float *C, int N)
{
    int correct = 1;
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if (C[i * N + j] != N)
            {
                printf(
                    "C(%d, %d) = %f is incorrect; C(%d, %d) should be %d\n",
                    i, j, C[i * N + j], i, j, N);
                correct = 0;
                break;
            }
        }
    }
    if (correct)
    {
        printf("The result is correct!\n\n");
    }
    else
    {
        printf("The result is not correct!\n\n");
    }
}

int main(int argc, char const *argv[])
{
    if (argc != 2)
    {
        printUsage(argc, argv);
        return 0;
    }
    int N = std::atoi(argv[1]);
    const int B1 = 64;
    const int B2 = 256;

    // Allocate and initialize the matrix A and vectors x, b
    // Allouer et initialiser la matrice A et matrices x, b
    float *A = (float *)malloc(N * N * sizeof(float));
    float *B = (float *)malloc(N * N * sizeof(float));
    float *C = (float *)malloc(N * N * sizeof(float));
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            A[i * N + j] = 1.0f;
            B[i * N + j] = 1.0f;
            C[i * N + j] = 0.0f;
        }
    }

    // Sequential and scalar matrix-matrix multiplication code with loop order i->j->k
    // Code sequentiel et scalaire produit matrice-matrice avec l'ordre de boucles i->j->k
    {
        auto start = std::chrono::high_resolution_clock::now();
        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j < N; j++)
            {
                for (int k = 0; k < N; k++)
                {
                    C[i * N + j] += A[i * N + k] * B[k * N + j]; // C(i, j) = C(i, j) + A(i, k) * B(k, j)
                }
            }
        }
        std::chrono::duration<double> timeDiff = std::chrono::high_resolution_clock::now() - start;
        std::cout << std::scientific << "Sequential scalar matmat i->j->k took " << timeDiff.count() << "s." << std::endl;
        std::cout << std::fixed << std::setprecision(2) << "Performance: " << 2.0 * N * N * (N - 1) / ((1e6) * timeDiff.count()) << "Mflops/s" << std::endl;
        verify(A, B, C, N);
    }

    // Sequential and scalar matrix-matrix multiplication code with loop order i->k->j
    // Code sequentiel et scalaire produit matrice-matrice avec l'ordre de boucles i->k->j
    {
        auto start = std::chrono::high_resolution_clock::now();
        for (int repet = 0; repet < NREPET; repet++)
        {
            memset(&C[0], 0, N * N * sizeof(float));
            for (int i = 0; i < N; i++)
            {
                for (int k = 0; k < N; k++)
                {
                    for (int j = 0; j < N; j++)
                    {
                        C[i * N + j] += A[i * N + k] * B[k * N + j]; // C(i, j) = C(i, j) + A(i, k) * B(k, j)
                    }
                }
            }
        }
        std::chrono::duration<double> timeDiff = (std::chrono::high_resolution_clock::now() - start) / NREPET;
        std::cout << std::scientific << "Sequential scalar matmat i->k->j took " << timeDiff.count() << "s." << std::endl;
        std::cout << std::fixed << std::setprecision(2) << "Performance: " << 2.0 * N * N * (N - 1) / ((1e6) * timeDiff.count()) << "Mflops/s" << std::endl;
        verify(A, B, C, N);
    }

    // Sequential and scalar matrix-matrix multiplication code with loop order i->k->j and single level tiling
    // Code sequentiel et scalaire produit matrice-matrice avec l'ordre de boucles i->k->j et tuilage d'un niveau
    {
        auto start = std::chrono::high_resolution_clock::now();
        for (int repet = 0; repet < NREPET; repet++)
        {
            memset(&C[0], 0, N * N * sizeof(float));
            for (size_t i = 0; i < N; i += B2)
            {
                for (size_t k = 0; k < N; k += B2)
                {
                    for (size_t j = 0; j < N; j += B2)
                    {
                        float *At = A + i * N + k;
                        float *Bt = B + k * N + j;
                        float *Ct = C + i * N + j;

                        for (size_t i2 = 0; i2 < B2; i2++)
                        {
                            for (size_t k2 = 0; k2 < B2; k2++)
                            {
                                for (size_t j2 = 0; j2 < B2; j2++)
                                {
                                    Ct[i2 * N + j2] += At[i2 * N + k2] * Bt[k2 * N + j2];
                                }
                            }
                        }
                    }
                }
            }
        }

        std::chrono::duration<double> timeDiff = (std::chrono::high_resolution_clock::now() - start) / NREPET;
        std::cout << std::scientific << "Single tile scalar matmat i->k->j took " << timeDiff.count() << "s." << std::endl;
        std::cout << std::fixed << std::setprecision(2) << "Performance: " << 2.0 * N * N * (N - 1) / ((1e6) * timeDiff.count()) << "Mflops/s" << std::endl;
        verify(A, B, C, N);
    }

    // Sequential and scalar matrix-matrix multiplication code with loop order i->k->j and two level tiling
    // Code sequentiel et scalaire produit matrice-matrice avec l'ordre de boucles i->k->j et tuilage de deux niveaux
    {
        auto start = std::chrono::high_resolution_clock::now();
        for (int repet = 0; repet < NREPET; repet++)
        {
            memset(&C[0], 0, N * N * sizeof(float));

            __m256 

            for (size_t i = 0; i < N; i += B2)
            {
                for (size_t k = 0; k < N; k += B2)
                {
                    for (size_t j = 0; j < N; j += B2)
                    {
                        // iterate on nested tile
                        for (size_t i2 = 0; i2 < B2; i2 += B1)
                        {
                            for (size_t k2 = 0; k2 < B2; k2 += B1)
                            {
                                for (size_t j2 = 0; j2 < B2; j2 += B1)
                                {
                                    float *tA = A + (i + i2) * N + k + k2;
                                    float *tB = B + (k + k2) * N + j + j2;
                                    float *tC = C + (i + i2) * N + j + j2;

                                    for (size_t i3 = 0; i3 < B1; i3++)
                                    {
                                        for (size_t k3 = 0; k3 < B1; k3++)
                                        {
                                            for (size_t j3 = 0; j3 < B1; j3++)
                                            {
                                                tC[i3 * N + j3] += tA[i3 * N + k3] * tB[k3 * N + j3];
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        std::chrono::duration<double> timeDiff = (std::chrono::high_resolution_clock::now() - start) / NREPET;
        std::cout << std::scientific << "Double tile scalar matmat i->k->j took " << timeDiff.count() << "s." << std::endl;
        std::cout << std::fixed << std::setprecision(2) << "Performance: " << 2.0 * N * N * (N - 1) / ((1e6) * timeDiff.count()) << "Mflops/s" << std::endl;
        verify(A, B, C, N);
    }

    // Vectorized matrix-matrix multiplication code with loop order i->k->j and two level tiling + AVX
    // Produit matrice-matrice vectorise avec l'ordre de boucles i->k->j et tuilage de deux niveaux + AVX
    {
        auto start = std::chrono::high_resolution_clock::now();
        for (int repet = 0; repet < NREPET; repet++)
        {
            memset(&C[0], 0, N * N * sizeof(float));
            for (size_t i = 0; i < N; i += B2)
            {
                for (size_t k = 0; k < N; k += B2)
                {
                    for (size_t j = 0; j < N; j += B2)
                    {
                        // iterate on nested tile
                        for (size_t i2 = 0; i2 < B2; i2 += B1)
                        {
                            for (size_t k2 = 0; k2 < B2; k2 += B1)
                            {
                                for (size_t j2 = 0; j2 < B2; j2 += B1)
                                {
                                    float *tA = A + (i + i2) * N + k + k2;
                                    float *tB = B + (k + k2) * N + j + j2;
                                    float *tC = C + (i + i2) * N + j + j2;

                                    for (size_t i3 = 0; i3 < B1; i3++)
                                    {
                                        for (size_t k3 = 0; k3 < B1; k3++)
                                        {
                                            for (size_t j3 = 0; j3 < B1; j3++)
                                            {
                                                tC[i3 * N + j3] += tA[i3 * N + k3] * tB[k3 * N + j3];
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        std::chrono::duration<double> timeDiff = (std::chrono::high_resolution_clock::now() - start) / NREPET;
        std::cout << std::scientific << "Double tile AVX matmat i->k->j took " << timeDiff.count() << "s." << std::endl;
        std::cout << std::fixed << std::setprecision(2) << "Performance: " << 2.0 * N * N * (N - 1) / ((1e6) * timeDiff.count()) << "Mflops/s" << std::endl;
        verify(A, B, C, N);
    }
}

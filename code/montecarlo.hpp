#ifndef _MONTECARLO_prog_
#define _MONTECARLO_prog_

const double pi = 3.141592653589793;

#include "randnumgen.hpp"
#include <chrono>
#include <initializer_list>
#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <fstream>
#include <thread>
#include <iomanip>
#include <algorithm>
#include <mutex>
#include <complex>

using std::cout;
using std::endl;
using ntype = double;
struct complex_n
{
    double a;
    double b;
};
struct local_H
{
    double H_x;
    double H_y;
};
struct index_FNN
{
    int i;
    int j;
    int k;
};
struct FNN_index
{
    index_FNN left;
    index_FNN right;
    index_FNN up;
    index_FNN down;
    index_FNN front;
    index_FNN behind;
};
struct FNN_index_trasformed
{
    int left;
    int right;
    int up;
    int down;
    int front;
    int behind;
};
struct metropolis_output
{
    int b;
    double v;
};
struct config
{
    int L;
    int mc_moves;
    int N_T;
    double fraction_of_pi;
    double T_min;
    double T_max;
    int N_sample;
};
struct metropolis_output2
{
    int acc;
    std::vector<double> E_i_story;
};
struct J_sparse_element
{
    int n;
    double value;
};
struct s_chi
{
    double _0;
    double min;
};
void printConfig(const config &configurazione)
{
    std::cout << "Configurazione letta:\n";
    std::cout << "Spin per lato (L): " << configurazione.L << "\n";
    std::cout << "Mosse Monte Carlo: " << configurazione.mc_moves << "\n";
    std::cout << "Numero di sistemi (N_T): " << configurazione.N_T << "\n";
    std::cout << "Frazione di pi greco: " << configurazione.fraction_of_pi << "\n";
    std::cout << "Temperatura minima (T_min): " << configurazione.T_min << "\n";
    std::cout << "Temperatura massima (T_max): " << configurazione.T_max << "\n"; 
    std::cout << "Numero di campioni: " << configurazione.N_sample << "\n";

}
void showProgressBar(int progress, int total, double elapsed_seconds, double estimated_total_seconds)
{
    int barWidth = 70;

    float progressRatio = static_cast<float>(progress) / total;

    std::cout << "\r["; // Usa '\r' per tornare all'inizio della stessa riga
    int pos = static_cast<int>(barWidth * progressRatio);
    for (int i = 0; i < barWidth; ++i)
    {
        if (i < pos)
            std::cout << "=";
        else if (i == pos)
            std::cout << ">";
        else
            std::cout << " ";
    }
    std::cout << "] " << int(progressRatio * 100.0) << "%";

    // Calcolo del tempo rimanente
    double remaining_seconds = estimated_total_seconds - elapsed_seconds;
    if (remaining_seconds < 0)
        remaining_seconds = 0;

    std::cout << " ETA: " << std::fixed << std::setprecision(1) << remaining_seconds << "s   ";

    // Forza l'output immediato
    std::cout.flush();
}
double R(int i, int j, int k)
{
    double distance = 0;
    distance = pow(i * i + j * j + k * k, 0.5);
    return distance;
}
double dot_product(double v_x, double v_y, double v_z, double k_x, double k_y, double k_z)
{
    return v_x * k_x + v_y * k_y + v_z * k_z;
}
config uploadConfig()
{
    config configurazione;
    std::ifstream configFile("config.txt");

    std::unordered_map<std::string, std::string> configMap; // Crea una mappa non ordinata dove ogni chiave (std::string) è associata a un valore (std::string). Questa mappa servirà per memorizzare le coppie chiave=valore dal file.
    std::string line;                                       // Crea una stringa per memorizzare ogni riga del file
    while (std::getline(configFile, line))                  // Legge ogni riga del file e la memorizza nella stringa line
    {
        std::istringstream iss(line);                                // Crea un oggetto di tipo std::istringstream per poter leggere la stringa line
        std::string key, value;                                      // Crea due stringhe per memorizzare la chiave e il valore
        if (std::getline(iss, key, '=') && std::getline(iss, value)) // Legge la stringa iss fino al carattere '=' e memorizza il risultato nella stringa key, poi legge il resto della stringa iss e memorizza il risultato nella stringa value
        {
            configMap[key] = value;
        }
    }

    configurazione.L = std::stoi(configMap["spin_per_lato"]);
    configurazione.mc_moves = std::stoi(configMap["mosse_montecarlo"]);
    configurazione.N_T = std::stoi(configMap["numero_sistemi"]);
    configurazione.fraction_of_pi = std::stod(configMap["fraction_of_pi"]);
    configurazione.T_min = std::stod(configMap["T_min"]);
    configurazione.T_max = std::stod(configMap["T_max"]);
    configurazione.N_sample = std::stoi(configMap["N_sample"]);

    configFile.close();

    return configurazione;
}

class sys
{
public:
    std::vector<std::vector<std::vector<double>>> spin;  // spins' value
    int L;                                               // lattice size
    std::vector<std::vector<double>> J;                  // whole J matrix
    std::vector<std::vector<J_sparse_element>> J_sparse; // sparse J matrix
    sys(int L_)
        : L(L_), spin(L_, std::vector<std::vector<double>>(L_, std::vector<double>(L_, 0.0))), J(L_ * L_ * L_, std::vector<double>(L_ * L_ * L_, 0.0)), J_sparse(L_ * L_ * L_, std::vector<J_sparse_element>(6))
    {
    }
    ~sys()
    {
    }
    // Inizializza la matrice spin con valori casuali tra 0 e 2pi
    void initial_condition()
    {
        for (int k = 0; k < L; ++k)
        {
            for (int j = 0; j < L; ++j)
            {
                for (int i = 0; i < L; ++i)
                {
                    spin[i][j][k] = rng.ranf_double(0, 2 * pi); // Usa il membro `spin` per salvare i valori
                }
            }
        }
        return;
    }
    // Stampa i valori della matrice spin
    void visualize_spin()
    {
        for (int k = 0; k < L; ++k)
        {
            for (int j = 0; j < L; ++j)
            {
                for (int i = 0; i < L; ++i)
                {
                    std::cout << "spin[" << i << "][" << j << "][" << k << "] = " << spin[i][j][k] << std::endl;
                }
            }
        }
        return;
    }
    // Stampa i valori della matrice spin su file
    void write_spin()
    {
        std::ofstream outfile("spin_data.txt"); // Crea un oggetto di tipo std::ofstream per scrivere su file

        for (int k = 0; k < L; ++k)
        {
            for (int j = 0; j < L; ++j)
            {
                for (int i = 0; i < L; ++i)
                {
                    if (outfile.is_open())
                    {
                        outfile << i << " " << j << " " << k << " " << spin[i][j][k] << "\n";
                    }
                }
            }
        }

        outfile.close();
    }
    // Stampa i valori della matrice J su file
    void write_J()
    {
        std::ofstream outfile("J_data.txt");

        for (int j = 0; j < L * L * L; ++j)
        {
            for (int i = 0; i < L * L * L; ++i)
            {
                if (outfile.is_open())
                {
                    outfile << J[i][j] << " ";
                }
            }
            outfile << "\n";
        }

        outfile.close();
    }
    // Genera un numero casuale distribuito secondo una gaussiana
    ntype gauss_pdf(void)
    {
        // gaussian of variance 1 and 0 mean
        double x1, x2;
        do
        {
            x1 = rng.ranf();
            x2 = rng.ranf();
        } while (x1 == 0 || x2 == 0);
        return cos(2 * pi * x2) * sqrt(-2.0 * log(x1));
    }
    // Mi trova i primi vicini di un certo spin
    FNN_index find_FNN_indexes(int i, int j, int k)
    {
        index_FNN temp;
        FNN_index result;

        int i_prime = i;
        int j_prime = j;
        int k_prime = k;
        // left
        if (i == 0)
            i_prime = L;

        result.left.i = i_prime - 1;
        result.left.j = j;
        result.left.k = k;

        // right
        i_prime = i;
        if (i == (L - 1))
            i_prime = -1;

        result.right.i = i_prime + 1;
        result.right.j = j;
        result.right.k = k;
        i_prime = i;

        // down
        if (k == 0)
            k_prime = L;

        result.down.i = i;
        result.down.j = j;
        result.down.k = k_prime - 1;

        k_prime = k;
        // up
        if (k == (L - 1))
            k_prime = -1;

        result.up.i = i;
        result.up.j = j;
        result.up.k = k_prime + 1;
        k_prime = k;

        // behind
        if (j == 0)
            j_prime = (L);

        result.behind.i = i;
        result.behind.j = j_prime - 1;
        result.behind.k = k;

        j_prime = j;

        // front
        if (j == (L - 1))
            j_prime = -1;

        result.front.i = i;
        result.front.j = j_prime + 1;
        result.front.k = k;
        j_prime = j;

        // std::cout<<"left FNN: "<<result.left.i<<result.left.j<<result.left.k<<std::endl;
        // std::cout<<"right FNN: "<<result.right.i<<result.right.j<<result.right.k<<std::endl;
        // std::cout<<"up FNN: "<<result.up.i<<result.up.j<<result.up.k<<std::endl;
        // std::cout<<"down FNN: "<<result.down.i<<result.down.j<<result.down.k<<std::endl;
        // std::cout<<"front FNN: "<<result.front.i<<result.front.j<<result.front.k<<std::endl;
        // std::cout<<"behind FNN: "<<result.behind.i<<result.behind.j<<result.behind.k<<std::endl;
        return result;
    }
    // Applico la trasformazione degli indici per poter usare una matrice 1D
    int index_transform(int i, int j, int k)
    {
        return L * L * i + L * j + k;
    }
    FNN_index_trasformed find_FNN_indexes_transformed(FNN_index FNN)
    {
        FNN_index_trasformed result;
        // IMPLEMENTARE CON FUNZIONE INDEX_TRANSFORM
        result.behind = FNN.behind.i * L * L + FNN.behind.j * L + FNN.behind.k;
        result.front = FNN.front.i * L * L + FNN.front.j * L + FNN.front.k;
        result.up = FNN.up.i * L * L + FNN.up.j * L + FNN.up.k;
        result.down = FNN.down.i * L * L + FNN.down.j * L + FNN.down.k;
        result.left = FNN.left.i * L * L + FNN.left.j * L + FNN.left.k;
        result.right = FNN.right.i * L * L + FNN.right.j * L + FNN.right.k;

        return result;
    }

    void set_J()
    {
        // I had to create a matrix like this because, in general, I might have had 1/2/../N neighbour interaction.
        // So what I did is unwrap the matrix of spin into a 1D array of dimensione L^3. I have to write down the interaction for
        // each spin. The selected spin moves on one row, which has L^3 elements. So any spin can interacts with all spins. That's the reason due to I need a L^3 x L^3 matrix ( in general ).
        // In the next method I will create a sparse matrix to save memory and computational time.
        int index;
        FNN_index_trasformed result;
        for (int i = 0; i < L; ++i)
        {
            for (int j = 0; j < L; ++j)
            {
                for (int k = 0; k < L; ++k)
                {
                    result = find_FNN_indexes_transformed(find_FNN_indexes(i, j, k));
                    index = index_transform(i, j, k);
                    if (index <= result.left)
                    {
                        J[index][result.left] = gauss_pdf();
                        // J[index][result.left] = 1;
                        J[result.left][index] = J[index][result.left];
                    }
                    if (index <= result.right)
                    {
                        J[index_transform(i, j, k)][result.right] = gauss_pdf();
                        // J[index_transform(i, j, k)][result.right] = 1;
                        J[result.right][index] = J[index][result.right];
                    }
                    if (index <= result.up)
                    {
                        J[index_transform(i, j, k)][result.up] = gauss_pdf();
                        // J[index_transform(i, j, k)][result.up] = 1;
                        J[result.up][index] = J[index][result.up];
                    }
                    if (index <= result.down)
                    {
                        J[index_transform(i, j, k)][result.down] = gauss_pdf();
                        // J[index_transform(i, j, k)][result.down] = 1;
                        J[result.down][index] = J[index][result.down];
                    }
                    if (index <= result.behind)
                    {
                        J[index_transform(i, j, k)][result.behind] = gauss_pdf();
                        // J[index_transform(i, j, k)][result.behind] = 1;
                        J[result.behind][index] = J[index][result.behind];
                    }
                    if (index <= result.front)
                    {
                        J[index_transform(i, j, k)][result.front] = gauss_pdf();
                        // J[index_transform(i, j, k)][result.front] = 1;
                        J[result.front][index] = J[index][result.front];
                    }
                }
            }
        }
        return;
    }
    void set_sparse_J()
    {
        int counter = 0;
        for (int i = 0; i < L * L * L; i++)
        {
            for (int j = 0; j < L * L * L; j++)
            {
                if (J[i][j] != 0)
                {
                    J_sparse[i][counter].n = j;
                    J_sparse[i][counter].value = J[i][j];
                    counter++;
                }
            }
            counter = 0;
        }
    }
    // Visualizza la matrice J
    void visualize_J()
    {
        int dim = L * L * L;
        for (int i = 0; i < dim; ++i)
        {
            for (int j = 0; j < dim; ++j)
            {
                std::cout << std::setw(10) << J[i][j] << " "; // setw() sets the width of the next element to be printed
            }
            std::cout << std::endl;
            std::cout << std::endl;
        }
        return;
    }
    // Visualizza la matrice J_sparse
    void visualize_sparse_J()
    {
        int dim = L * L * L;
        for (int i = 0; i < dim; ++i)
        {
            for (int j = 0; j < 6; ++j)
            {
                std::cout << std::setw(10) << J_sparse[i][j].value << " ";
            }
            std::cout << std::endl;
            std::cout << std::endl;
        }
        return;
    }
    // original H(i,j,k) with whole J
    /*
    local_H H(int i, int j, int k)
    {
        local_H H_ijk;
        int index;
        int dim = L * L * L;
        int x, y, z;

        index = index_transform(k, j, i);
        // std::cout<<"Spin selezionato "<<index<< " spin"<<std::endl;

        H_ijk.H_x = 0.0;
        H_ijk.H_y = 0.0;

        for (int m = 0; m < dim; ++m)
        {
            z = std::floor(m / (L * L));
            y = std::floor((m % (L * L)) / L);
            x = m % L;
            // std::cout<<"Dovrebbero uscirmi le coordinate del "<<m<< " spin"<<std::endl;
            // std::cout<<"i: "<<x<<std::endl;
            // std::cout<<"j: "<<y<<std::endl;
            // std::cout<<"k: "<<z<<std::endl;

            H_ijk.H_x += J[index][m] * cos(spin[x][y][z]);
            H_ijk.H_y += J[index][m] * sin(spin[x][y][z]);
        }

        return H_ijk;
    }
    */

    // new H(i,j,k) with J_sparse
    local_H H(int i, int j, int k)
    {
        local_H H_ijk;
        int index;
        int dim = L * L * L;
        int x, y, z;

        index = index_transform(k, j, i);

        H_ijk.H_x = 0.0;
        H_ijk.H_y = 0.0;

        for (int m = 0; m < 6; ++m)
        {
            z = std::floor(J_sparse[index][m].n / (L * L));
            y = std::floor((J_sparse[index][m].n % (L * L)) / L);
            x = J_sparse[index][m].n % L;

            H_ijk.H_x += J_sparse[index][m].value * cos(spin[x][y][z]);
            H_ijk.H_y += J_sparse[index][m].value * sin(spin[x][y][z]);
        }

        return H_ijk;
    }
    // old method to calculate energy for each spin
    /*
        double E_i(int i, int j, int k)
        {
            double E=0;
            int index;
            int dim = L * L * L;
            int x, y, z;

            index = index_transform(k, j, i);
            for (int m = 0; m < dim; ++m)
            {
                z = std::floor(m / (L * L));
                y = std::floor((m % (L * L)) / L);
                x = m % L;
                // std::cout<<"Dovrebbero uscirmi le coordinate del "<<m<< " spin"<<std::endl;
                // std::cout<<"i: "<<x<<std::endl;
                // std::cout<<"j: "<<y<<std::endl;
                // std::cout<<"k: "<<z<<std::endl;

                E += J[index][m] * cos(spin[i][j][k] - spin[x][y][z]);
            }

            return -E;
        }
        double total_e()
        {
            double E = 0;
            for (int i = 0; i < L; ++i)
            {
                for (int j = 0; j < L; ++j)
                {
                    for (int k = 0; k < L; ++k)
                    {
                        E+= E_i(i,j,k);
                    }
                }
            }
            return E/2;
        }*/

    // Starting  method total_e() for the most general way on n- neighours interactions
    /*
        double total_e()
        {
            double E = 0;
            int x, y, z;
            int x_prime, y_prime, z_prime;
            for (int m = 0; m < L * L * L; m++) // I sum only over upper triangular matrix such that I don't count same terms twice.
            {
                z_prime = std::floor(m / (L * L));
                y_prime = std::floor((m % (L * L)) / L);
                x_prime = m % L;
                for (int n = 0; n < m; n++)
                {
                    z = std::floor(n / (L * L));
                    y = std::floor((n % (L * L)) / L);
                    x = n % L;

                    E += J[m][n] * cos(spin[x_prime][y_prime][z_prime] - spin[x][y][z]);
                }
            }
            return -E;
        }

    */

    // total_e() optimized for J_sparse
    double total_e()
    {
        double E = 0;
        int x, y, z;
        int x_prime, y_prime, z_prime;
        for (int i = 0; i < L * L * L; i++) // I sum only over upper triangular matrix such that I don't count same terms twice.
        {
            z_prime = std::floor(i / (L * L));
            y_prime = std::floor((i % (L * L)) / L);
            x_prime = i % L;
            for (int j = 0; j < 6; j++)
            {
                z = std::floor(J_sparse[i][j].n / (L * L));
                y = std::floor((J_sparse[i][j].n % (L * L)) / L);
                x = J_sparse[i][j].n % L;

                E += J_sparse[i][j].value * cos(spin[x_prime][y_prime][z_prime] - spin[x][y][z]);
            }
        }
        return -E / 2;
    }

    metropolis_output metropolis_singlestep(double theta, double T)
    {
        double V_o = 0, V_n = 0;
        int choose_x, choose_y, choose_z;
        double temp_o;
        int x, y, z;
        double delta_theta;
        double delta_V = 0;
        double p = 0;
        double rnd = 0;
        int b = 0;
        metropolis_output risultati;

        V_o = total_e();

        // std::cout << "v_o: " << V_o << std::endl;

        choose_x = rng.ranf_int(0, L);
        choose_y = rng.ranf_int(0, L);
        choose_z = rng.ranf_int(0, L);

        temp_o = spin[choose_x][choose_y][choose_z];
        // std::cout << "E' stato scelto lo spin: " << choose_x << " " << choose_y << " " << choose_z << " con valorE: " << spin[choose_x][choose_y][choose_z] << std::endl;

        delta_theta = rng.ranf_double(-1.0, 1.0) * theta;
        spin[choose_x][choose_y][choose_z] = spin[choose_x][choose_y][choose_z] + delta_theta;

        V_n = total_e();
        // std::cout << "v_n: " << V_n << std::endl;

        delta_V = V_n - V_o;

        if (delta_V >= 0)
        {
            p = exp(-(1 / T) * delta_V);
            // std::cout << "exp: " << p << std::endl;
            rnd = rng.ranf();
            // std::cout << "rnd: " << rnd << std::endl;

            if (rnd <= p)
            {
                b = 1;
            }
            else
                b = 0;
        }
        else
        {
            b = 1;
        }

        if (b == 0)
        {
            spin[choose_x][choose_y][choose_z] = temp_o;
        }
        // std::cout << "Esito della mossa: " << b << std::endl;
        if (b == 0)
        {
            risultati.v = V_o;
            risultati.b = b;
        }
        else
        {
            risultati.v = V_n;
            risultati.b = b;
        }
        return risultati; // Mi ritorna il valore dell'energia e se la mossa è stata accettata o meno (con 1 o 0)
    }
    template <typename T>
    // Scrive un vettore su file
    void write_vector_to_file(const std::vector<T> vec, const std::string filename)
    {
        std::ofstream file(filename);

        if (!file)
        {
            std::cerr << "Errore nell'apertura del file!" << std::endl; // Stampa un messaggio di errore se non è possibile aprire il file
            return;
        }

        for (const T &valore : vec) // Per ogni valore nel vettore
        {
            file << valore << std::endl;
        }

        file.close();
    }
    // Legge una configurazione di spin da file
    void read_spin_from_file(const std::string &filename)
    {
        std::ifstream infile(filename);
        if (!infile.is_open())
        {
            std::cerr << "Errore: impossibile aprire il file " << filename << std::endl;
            return;
        }

        int x, y, z;  // Coordinate
        double value; // Valore dello spin
        while (infile >> x >> y >> z >> value)
        {
            spin[x][y][z] = value; // Aggiorna il valore nella matrice
        }

        infile.close();
    }
    // Legge la matrice J da file
    void read_J_from_file(const std::string &filename)
    {
        std::ifstream infile(filename);
        if (!infile.is_open())
        {
            std::cerr << "Errore: impossibile aprire il file " << filename << std::endl;
            return;
        }
        for (int i = 0; i < L; ++i)
        {
            for (int j = 0; j < L; ++j)
            {
                infile >> J[i][j];
            }
        }

        infile.close();
    }
    void metropolis(double theta, double T, int mc_step, std::string filename)
    {
        int acc = 0;
        double eff = 0;
        metropolis_output results;
        auto start_time = std::chrono::steady_clock::now();
        std::vector<double> v_o_story; // energy story

        for (int i = 0; i < mc_step; i++)
        {
            results = metropolis_singlestep(theta, T);
            acc += results.b;
            v_o_story.push_back(results.v);

            auto current_time = std::chrono::steady_clock::now();
            std::chrono::duration<double> elapsed = current_time - start_time;

            // Stima del tempo totale
            double elapsed_seconds = elapsed.count();
            double estimated_total_seconds = (i > 0) ? (elapsed_seconds / i) * mc_step : 0; // Stima del tempo totale in base al tempo impiegato per le iterazioni precedenti
            showProgressBar(i, mc_step, elapsed_seconds, estimated_total_seconds);
        }

        // std::cout << "mosse accettate: " << acc << std::endl;
        eff = static_cast<double>(acc) / mc_step;
        write_vector_to_file(v_o_story, filename);

        std::cout << "Frazioni di mosse accettate: " << std::fixed << std::setprecision(3) << eff << std::endl;
        auto end_time = std::chrono::steady_clock::now();

        std::chrono::duration<double> total_elapsed = end_time - start_time;
        std::cout << "Esecuzione completata in " << total_elapsed.count() << " secondi." << std::endl;
        std::cout << "\n"
                  << endl;
        return;
    }

    void over_relaxation_singlestep(int i, int j, int k)
    {
        double s_prime_x = 0;
        double s_prime_y = 0;
        local_H H_i;
        H_i = H(i, j, k);
        s_prime_x = -cos(spin[i][j][k]) + 2 * H_i.H_x / (H_i.H_x * H_i.H_x + H_i.H_y * H_i.H_y) * (cos(spin[i][j][k]) * H_i.H_x + sin(spin[i][j][k]) * H_i.H_y);
        s_prime_y = -sin(spin[i][j][k]) + 2 * H_i.H_y / (H_i.H_x * H_i.H_x + H_i.H_y * H_i.H_y) * (cos(spin[i][j][k]) * H_i.H_x + sin(spin[i][j][k]) * H_i.H_y);
        // std::cout<<" s e' "<<spin[i][j][k]<<std::endl;
        // std::cout<<"componente x di s_prime e' "<<s_prime_x<<std::endl;
        // std::cout<<"componente y di s_prime e' "<<s_prime_y<<std::endl;
        // std::cout<<"Il modulo di s_prime e' "<<s_prime_x*s_prime_x+s_prime_y*s_prime_y<<std::endl;
        spin[i][j][k] = atan2(s_prime_y, s_prime_x); // DIO BENEDICA CHATGPT CHE MI HA SALVATO DA ALTRE 7 ORE DI DEBUGGING!! ATAN2() CONTEMPLA TUTTI I CASI NON COME ATAN()!!!!
        // std::cout << " s_prime e' " << spin[i][j][k] << std::endl;
        return;
    };

    void over_relaxation()
    {
        for (int k = 0; k < L; ++k)
        {
            for (int j = 0; j < L; ++j)
            {
                for (int i = 0; i < L; ++i)
                {
                    over_relaxation_singlestep(i, j, k);
                }
            }
        }

        return;
    }
    metropolis_output2 mmetropolis(double theta, double T, std::string filename)
    {
        int acc = 0;
        double eff = 0;
        metropolis_output results;
        std::vector<double> v_o_story; // energy story
        metropolis_output2 results2;
        results2.acc = 0;

        for (int i = 0; i < L * L * L; i++)
        {
            results = metropolis_singlestep(theta, T);
            results2.acc += results.b;
            v_o_story.push_back(results.v); // RIATTIVARE PER SALVARE L'ENERGIA
        }
        results2.E_i_story = v_o_story;
        return results2;
    }

    // METODO RUN BUONO

    void run(double theta, double T, int mc_moves, const std::string &filename)
    {
        metropolis_output results;

        int acc = 0;
        double eff = 0.0;
        auto start_time = std::chrono::steady_clock::now();

        std::vector<double> v_o_story(L * L * L * mc_moves, 0.0);

        for (int i = 0; i < mc_moves; ++i)
        {

            for (int j = 0; j < 10; ++j)
            {
                over_relaxation();
            }

            for (int k = 0; k < L * L * L; ++k)
            {
                results = metropolis_singlestep(theta, T);
                v_o_story[i * L * L * L + k] = results.v;

                acc += results.b;
            }

            // Calcola e mostra la barra di progresso
            auto current_time = std::chrono::steady_clock::now();
            std::chrono::duration<double> elapsed = current_time - start_time;

            double elapsed_seconds = elapsed.count();
            double estimated_total_seconds = (i > 0) ? (elapsed_seconds / i) * mc_moves : 0;
            showProgressBar(i + 1, mc_moves, elapsed_seconds, estimated_total_seconds);
        }

        eff = static_cast<double>(acc) / (mc_moves * L * L * L);

        write_vector_to_file(v_o_story, filename);

        std::cout << "\n";
        std::cout << "Mosse accettate: " << acc << std::endl;
        std::cout << "Frazione di mosse accettate: " << std::fixed << std::setprecision(3) << eff << std::endl;

        auto end_time = std::chrono::steady_clock::now();
        std::chrono::duration<double> total_elapsed = end_time - start_time;
        std::cout << "Esecuzione completata in " << total_elapsed.count() << " secondi." << std::endl;

        return;
    }

    /*
        void run(double theta, double T, int mc_moves, const std::string &filename)
        {
            // Variabile per i risultati della Metropolis
            metropolis_output results;

            // Variabili per il progresso
            int acc = 0;      // Numero di mosse accettate
            double eff = 0.0; // Efficienza delle mosse accettate
            auto start_time = std::chrono::steady_clock::now();
            int mod = 0; // Ogni quanto calcolare la media
            int counter = 0;
            double E_prev = 0;
            double E = 0;

            mod = static_cast<int>(round(0.001 * mc_moves * L * L * L)); // Calcolo di mod

            // Parametro di convergenza
            double epsilon = 1;           // Soglia per la differenza delle medie
            double previous_average = 0.0; // Media precedente
            bool converged = false;        // Flag per convergenza

            // Inizializzazione del vettore per memorizzare le medie
            std::vector<double> averages; // Salverà le medie ogni "mod" passi

            double energy_sum = 0.0; // Variabile per sommare le energie
            int count_steps = 0;     // Contatore dei passi per la media

            // Ciclo principale Monte Carlo
            for (int i = 0; i < mc_moves; ++i)
            {

                // Ciclo di over-relaxation
                for (int j = 0; j < 10; ++j)
                {
                    over_relaxation();
                }

                // Ciclo per il calcolo delle configurazioni
                for (int k = 0; k < L * L * L; ++k)
                {
                    counter++;

                    results = metropolis_singlestep(theta, T);
                    energy_sum += results.v; // Accumula l'energia
                    count_steps++;

                    // Incrementa il numero di accettazioni se la mossa è accettata
                    acc += results.b;
                }

                // Calcolo della media ogni "mod" passi
                if (counter % mod == 0)
                {
                    double average_energy = energy_sum / count_steps; // Calcolo della media
                    averages.push_back(average_energy);               // Salva la media nel vettore

                    // Controllo convergenza
                    if (abs(average_energy - previous_average) < epsilon)
                    {
                        std::cout << "Convergenza raggiunta dopo " << counter << " passi Monte Carlo.\n";
                        converged = true;
                        break; // Interrompe il ciclo Monte Carlo
                    }

                    previous_average = average_energy; // Aggiorna la media precedente
                    energy_sum = 0.0;                  // Resetta la somma
                    count_steps = 0;                   // Resetta il contatore
                }

                // Calcola e mostra la barra di progresso
                // auto current_time = std::chrono::steady_clock::now();
                // std::chrono::duration<double> elapsed = current_time - start_time;

                // double elapsed_seconds = elapsed.count();
                // double estimated_total_seconds = (i > 0) ? (elapsed_seconds / i) * mc_moves : 0;
                // showProgressBar(i + 1, mc_moves, elapsed_seconds, estimated_total_seconds);

                // Se convergenza è stata raggiunta, esci
                if (converged)
                    break;
            }

            // Calcolo dell'efficienza delle mosse accettate
            eff = static_cast<double>(acc) / (mc_moves * L * L * L);

            // Scrive le medie dei risultati su file
            write_vector_to_file(averages, filename);

            // Output dei risultati
            std::cout << "\n";
            std::cout << "Mosse accettate: " << acc << std::endl;
            std::cout << "Frazione di mosse accettate: " << std::fixed << std::setprecision(3) << eff << std::endl;

            // Tempo totale di esecuzione
            auto end_time = std::chrono::steady_clock::now();
            std::chrono::duration<double> total_elapsed = end_time - start_time;
            std::cout << "Esecuzione completata in " << total_elapsed.count() << " secondi." << std::endl;
        }
    */
};

class replica
{
public:
    int N_T;                          // Numero di sistemi
    std::vector<sys> N_T_systems;     // Sistemi
    std::vector<double> temperatures; // Temperature associate
    std::vector<double> energies;     // Energie correnti di ciascun sistema
    int L;

    // Costruttore che inizializza N_T sistemi con LxLxL spin e temperature uniformemente distribuite tra T_min e T_max
    replica(int N_T, int L_, double T_min, double T_max, bool J_fixed = 1)
        : N_T(N_T), N_T_systems(N_T, sys(L_)), temperatures(N_T), energies(N_T, 0.0), L(L_)
    {
        // Inizializza le temperature in modo uniforme tra T_min e T_max
        for (int i = 0; i < N_T; ++i)
        {
            temperatures[i] = T_min + i * (T_max - T_min) / (N_T - 1);
        }

        if (J_fixed)
        {
            N_T_systems[0].set_J();
            N_T_systems[0].set_sparse_J();
            for (int i = 0; i < N_T; ++i)
            {
                N_T_systems[i].initial_condition();
                if (i != 0)
                    N_T_systems[i].J_sparse = N_T_systems[0].J_sparse;
            }
        }
        else
        {
            for (int i = 0; i < N_T; i++)
            {
                N_T_systems[i].set_J();
                N_T_systems[i].set_sparse_J();
                N_T_systems[i].initial_condition();
            }
        }
    }

    // Metodo per calcolare l'energia di ogni sistema
    void update_energies()
    {
        double current_energy = 0;
        std::cout << "\n"
                  << std::endl;

        for (int i = 0; i < N_T; ++i)
        {
            current_energy = N_T_systems[i].total_e();

            // std::cout << "Energia del sistema " << i + 1 << ": " << current_energy << std::endl;

            energies[i] = current_energy;
        }

        return;
    }
    // Metodo per scambiare due configurazioni in caso di accettazione di Parallel Tempering tra il sistema i e i+1
    void swap_configuration(int i)
    {
        sys temp(L);
        temp = N_T_systems[i];
        N_T_systems[i] = N_T_systems[i + 1];
        N_T_systems[i + 1] = temp;
    }
    // Metodo per un passo di Parallel Tempering
    int tempering_step(int i)
    {
        int b = 0;
        double p = 0;
        double rnd = 0;
        double beta_i = 1.0 / temperatures[i];
        double beta_i_1 = 1.0 / temperatures[i + 1];
        double delta = (beta_i_1 - beta_i) * (energies[i + 1] - energies[i]);
        // std::cout << "delta beta: " << std::fixed << std::setprecision(10) << beta_i_1 - beta_i << std::endl;
        // std::cout << "delta E: " << std::fixed << std::setprecision(10) << energies[i + 1] - energies[i] << std::endl;

        if (delta <= 0)
        {
            // std::cout << "delta: " << std::fixed << std::setprecision(10) << delta << std::endl;
            p = exp(delta);
            // std::cout << "exp: " << std::fixed << std::setprecision(10) << p << std::endl;
            rnd = rng.ranf();
            // std::cout << "rnd: " << std::fixed << std::setprecision(10) << rnd << std::endl;

            if (rnd <= p)
            {
                // cout << "\n";
                // cout << "Fatto lo swap per il sistema: " << i << endl;
                swap_configuration(i);
                b = 1;
            }
            else
            {
                b = 0;
                // cout << "\n";
                // cout << "Non fatto lo swap per il sistema: " << i << endl;
            }
        }
        else
        {
            b = 1;
            swap_configuration(i);
        }

        return b;
    }
    /*

        void evolve(double theta, int mc_moves, const std::string &base_filename)
        {
            metropolis_output2 temp, results2;
            double eff;
            for (int i = 0; i < N_T; i++)
            {
                std::string filename = base_filename + "_replica_" + std::to_string(i) + ".txt";

                temp.acc = 0;
                temp.E_i_story.clear();
                results2.acc = 0;
                results2.E_i_story.clear();
                eff = 0;

                auto start_time = std::chrono::steady_clock::now();

                std::cout << "\n";
                std::cout << "Evoluzione per la replica " << i + 1 << "/" << N_T << "...\n";

                // RICORDARSI DI CANCELLARE I FILE DI ENERGIA DA PRECEDENTI SIMULAZIONI O SI AGGIUNGONO

                // Scrivi i risultati in modalità "append" per evitare di riscrivere tutto ogni volta
                std::ofstream file_out(filename, std::ios::app);
                if (!file_out.is_open())
                {
                    throw std::runtime_error("Errore nell'aprire il file: " + filename);
                }

                for (int j = 0; j < mc_moves; j++)
                {
                    for (int k = 0; k < 10; k++)
                    {
                        N_T_systems[i].over_relaxation();
                    }

                    // Chiama la funzione mmetropolis e accumula i risultati
                    temp = N_T_systems[i].mmetropolis(theta, temperatures[i], mc_moves, filename);

                    // if(i<(N_T-1))
                    //     tempering_step(i);
                    //  Accumula i risultati solo temporaneamente per il batch corrente
                    results2.acc += temp.acc;

                    // Scrivi direttamente i nuovi risultati su file
                    for (const auto &E : temp.E_i_story)
                    {
                        file_out << E << "\n";
                    }

                    // Calcolo del progresso
                    auto current_time = std::chrono::steady_clock::now();
                    std::chrono::duration<double> elapsed = current_time - start_time;
                    double elapsed_seconds = elapsed.count();
                    double estimated_total_seconds = (j > 0) ? (elapsed_seconds / j) * mc_moves : 0;

                    showProgressBar(j + 1, mc_moves, elapsed_seconds, estimated_total_seconds);
                }

                file_out.close(); // Chiudi il file dopo aver terminato la replica

                auto end_time = std::chrono::steady_clock::now();
                std::chrono::duration<double> total_elapsed = end_time - start_time;

                eff = static_cast<double>(results2.acc) / (mc_moves * L * L * L);

                std::cout << "\n";
                std::cout << "\n";
                std::cout << "Mosse accettate: " << results2.acc << std::endl;
                std::cout << "Frazione di mosse accettate: " << std::fixed << std::setprecision(3) << eff << std::endl;

                std::cout << "\nReplica " << i + 1 << "/" << N_T << " completata in " << total_elapsed.count() << " secondi.\n";
            }

            update_energies(); // Aggiorna le energie globali
        }

    */
    // Metodo per evolvere le repliche in parallelo (su diversi thread della CPU/GPU)
    std::vector<int> parallel_evolve(double theta, const std::string &base_filename)
    {
        std::vector<std::thread> threads;      // Vettore di thread
        std::vector<int> accept_moves(N_T, 0); // Contatore mosse accettate per ogni replica

        auto evolve_replica = [&](int replica_index) // Funzione lambda per evolvere una replica
        {
            std::string filename = base_filename + "_sistema_" + std::to_string(replica_index) + ".txt";
            metropolis_output2 temp, results2;
            results2.acc = 0;
            results2.E_i_story.clear();

            // LEVARE IL COMMENTO SE VOGLIO STAMPARE LE ENERGIE SU FILE ********

            std::ofstream file_out(filename, std::ios::app);

            if (!file_out.is_open())
            {
                throw std::runtime_error("Errore nell'aprire il file: " + filename);
            }

            for (int k = 0; k < 10; k++)
            {
                N_T_systems[replica_index].over_relaxation();
            }

            temp = N_T_systems[replica_index].mmetropolis(theta, temperatures[replica_index], filename);
            results2.acc += temp.acc;

            for (const auto &E : temp.E_i_story) // Salva le energie su file
            {
                file_out << E << "\n";
            }

            accept_moves[replica_index] = results2.acc; // Salva il numero di mosse accettate

            file_out.close(); // LEVARE IL COMMENTO SE VUOI STAMPARE LE ENERGIE SU FILE ********
        };

        for (int i = 0; i < N_T; ++i)
        {
            threads.emplace_back(evolve_replica, i); // Crea un thread per ogni replica chiamando la funzione lambda
        }

        // Attende che tutti i thread abbiano finito
        for (auto &t : threads)
        {
            t.join();
        }

        return accept_moves;
    }

    void run(double theta, int mc_moves, const std::string &base_filename)
    {
        std::vector<int> total_accept_moves(N_T, 0); // Contatore totale mosse accettate per ogni sistema
        auto start_time = std::chrono::steady_clock::now();
        double eff = 0;
        int acc = 0;
        for (int step = 0; step < mc_moves; ++step)
        {
            std::vector<int> step_accept_moves = parallel_evolve(theta, base_filename); // fa una mossa montecalro per ogni sistema

            // Aggiorna il totale delle mosse accettate
            for (int i = 0; i < N_T; ++i)
            {
                total_accept_moves[i] += step_accept_moves[i];
            }

            // Esegui il tempering tra le repliche
            update_energies();

            for (int i = 0; i < N_T - 1; ++i)
            {
                acc += tempering_step(i);
            }

            auto current_time = std::chrono::steady_clock::now();
            std::chrono::duration<double> elapsed = current_time - start_time;

            double elapsed_seconds = elapsed.count();
            double estimated_total_seconds = (step > 0) ? (elapsed_seconds / step) * mc_moves : 0;
            showProgressBar(step + 1, mc_moves, elapsed_seconds, estimated_total_seconds);
        }
        double efficiency = static_cast<double>(acc) / (mc_moves * N_T);
        std::cout << "Frazioni mosse tempering accetate: " << std::fixed << std::setprecision(3) << efficiency << std::endl;

        update_energies();

        for (int i = 0; i < N_T; ++i)
        {
            double efficiency = static_cast<double>(total_accept_moves[i]) / (L * L * L * mc_moves);
            {
                std::cout << "Efficienza sistema " << i + 1 << ": " << std::fixed << std::setprecision(3) << efficiency << std::endl;
            }
        }
        auto end_time = std::chrono::steady_clock::now();
        std::chrono::duration<double> total_elapsed = end_time - start_time;
        std::cout << "Esecuzione completata in " << total_elapsed.count() << " secondi." << std::endl;
    }

    // run di un passo montecarlo per ogni replica SENZA barra di progresso
    void avg_run(double theta, int mc_moves, const std::string &base_filename)
    {
        std::vector<int> total_accept_moves(N_T, 0);

        for (int step = 0; step < mc_moves; ++step)
        {
            std::vector<int> step_accept_moves = parallel_evolve(theta, base_filename);

            for (int i = 0; i < N_T; ++i)
            {
                total_accept_moves[i] += step_accept_moves[i];
            }

            for (int i = 0; i < N_T - 1; ++i)
            {
                tempering_step(i);
            }
        }
    }
};
// Classe per le repliche di N_T sistemi
class replicas2
{
public:
    std::vector<replica> replicas;
    int L;
    replicas2(int N_T, int L_, double T_min, double T_max, bool J_fixed) : replicas(2, replica(N_T, L_, T_min, T_max, J_fixed)), L(L_)
    {
    }
    void set_equal_J()
    {
        for (int i = 0; i < replicas[0].N_T; i++)
        {
            replicas[1].N_T_systems[i].J_sparse = replicas[0].N_T_systems[i].J_sparse; 
        }
    }
    void run(double theta, int mc_moves, const std::string filename)
    {
        for (int i = 0; i < 2; i++)
        {
            cout << "Replica: " << i << endl;
            std::string updated_filename = filename + "replica_" + std::to_string(i) + "_";
            replicas[i].run(theta, mc_moves, updated_filename);
        }
    }
    // metodo run per calcolare le medie di energia di ogni replica senza barra di progresso
    void avg_run(double theta, int mc_moves, const std::string filename)
    {
        for (int i = 0; i < 2; i++)
        {
            std::string updated_filename = std::to_string(i) + "_" + filename;
            replicas[i].avg_run(theta, mc_moves, updated_filename);
        }
    }
    
    std::complex<double> q_EA(int mu, int nu, double k_x, double k_y, double k_z, int sys_index) 
    {
        std::complex<double> q_EA = 0;
        std::complex<double> phase = 0;

        if (mu == 0 && nu == 0)
        {
            for (int i = 0; i < L; i++)
            {
                for (int j = 0; j < L; j++)
                {
                    for (int k = 0; k < L; k++)
                    {
                        phase = std::exp(std::complex<double>(0.0, dot_product(k_x, k_y, k_z, i, j, k)));
                        q_EA += cos(replicas[0].N_T_systems[sys_index].spin[i][j][k]) * cos(replicas[1].N_T_systems[sys_index].spin[i][j][k]) * phase;
                    }
                }
            }
        }
        if (mu == 1 && nu == 0)
        {
            for (int i = 0; i < L; i++)
            {
                for (int j = 0; j < L; j++)
                {
                    for (int k = 0; k < L; k++)
                    {
                        phase = std::exp(std::complex<double>(0.0, dot_product(k_x, k_y, k_z, i, j, k)));
                        q_EA += sin(replicas[0].N_T_systems[sys_index].spin[i][j][k]) * cos(replicas[1].N_T_systems[sys_index].spin[i][j][k]) * phase;
                    }
                }
            }
        }

        if (mu == 0 && nu == 1)
        {
            for (int i = 0; i < L; i++)
            {
                for (int j = 0; j < L; j++)
                {
                    for (int k = 0; k < L; k++)
                    {
                        phase = std::exp(std::complex<double>(0.0, dot_product(k_x, k_y, k_z, i, j, k)));
                        q_EA += cos(replicas[0].N_T_systems[sys_index].spin[i][j][k]) * sin(replicas[1].N_T_systems[sys_index].spin[i][j][k]) * phase;
                    }
                }
            }
        }

        if (mu == 1 && nu == 1)
        {
            for (int i = 0; i < L; i++)
            {
                for (int j = 0; j < L; j++)
                {
                    for (int k = 0; k < L; k++)
                    {
                        phase = std::exp(std::complex<double>(0.0, dot_product(k_x, k_y, k_z, i, j, k)));
                        q_EA += sin(replicas[0].N_T_systems[sys_index].spin[i][j][k]) * sin(replicas[1].N_T_systems[sys_index].spin[i][j][k]) * phase;
                    }
                }
            }
        }
        return q_EA / double(L * L * L);
    }

    /*
         double susceptibility(int sys_index, double k_x, double k_y, double k_z)
         {
             double susceptibility = 0;
             std::complex<double> phase = 0;
             std::complex<double> q_EA_00 = 0, q_EA_01 = 0, q_EA_10 = 0, q_EA_11 = 0;

             for (int i = 0; i < L; i++)
             {
                 for (int j = 0; j < L; j++)
                 {
                     for (int k = 0; k < L; k++)
                     {
                         phase = std::exp(std::complex<double>(0.0, dot_product(k_x, k_y, k_z, i, j, k)));
                         q_EA_00 += cos(replicas[0].N_T_systems[sys_index].spin[i][j][k]) * cos(replicas[1].N_T_systems[sys_index].spin[i][j][k]) * phase;
                         q_EA_10 += sin(replicas[0].N_T_systems[sys_index].spin[i][j][k]) * cos(replicas[1].N_T_systems[sys_index].spin[i][j][k]) * phase;
                         q_EA_01 += cos(replicas[0].N_T_systems[sys_index].spin[i][j][k]) * sin(replicas[1].N_T_systems[sys_index].spin[i][j][k]) * phase;
                         q_EA_11 += sin(replicas[0].N_T_systems[sys_index].spin[i][j][k]) * sin(replicas[1].N_T_systems[sys_index].spin[i][j][k]) * phase;
                     }
                 }
             }
             susceptibility = double(L * L * L) * (q_EA_00.real()*q_EA_00.real() +q_EA_00.imag()*q_EA_00.imag() +q_EA_10.real()*q_EA_10.real()+q_EA_10.imag()*q_EA_10.imag() + q_EA_01.real()*q_EA_01.real() + q_EA_01.imag()*q_EA_01.imag() + q_EA_11.real()*q_EA_11.real() + q_EA_11.imag()*q_EA_11.imag());
             return susceptibility;
         }

         double xi_L(int sys_index)
         {
             double xi;
             double k_min = 2 * pi / L;
             double chi_0 = 0;
             double chi_min = 0;

             chi_0 = susceptibility(sys_index, 0, 0, 0);
             chi_min = susceptibility(sys_index, 2 * pi / L, 0, 0);
             cout << "valore chi_0: " << chi_0 << endl;
             cout << "valore chi_min: " << chi_min << endl;
             xi = 1 / (2 * sin(k_min / 2)) * pow((chi_0 / chi_min) - 1, 0.5);
             return xi;
         }
        */

    double mod_sq_q(int sys_index, double k_x, double k_y, double k_z)
    {
        double susceptibility = 0;
        double k_R = 0;
        complex_n q_00, q_01, q_10, q_11; 

        q_00.a = 0;
        q_00.b = 0;
        q_10.a = 0;
        q_10.b = 0;
        q_01.a = 0;
        q_01.b = 0;
        q_11.a = 0;
        q_11.b = 0;

        for (int i = 0; i < L; i++)
        {
            for (int j = 0; j < L; j++)
            {
                for (int k = 0; k < L; k++)
                {
                    k_R = dot_product(k_x, k_y, k_z, i, j, k);

                    q_00.a += cos(k_R) * (cos(replicas[0].N_T_systems[sys_index].spin[i][j][k]) * cos(replicas[1].N_T_systems[sys_index].spin[i][j][k]));
                    q_00.b += sin(k_R) * (cos(replicas[0].N_T_systems[sys_index].spin[i][j][k]) * cos(replicas[1].N_T_systems[sys_index].spin[i][j][k]));

                    q_10.a += cos(k_R) * (sin(replicas[0].N_T_systems[sys_index].spin[i][j][k]) * cos(replicas[1].N_T_systems[sys_index].spin[i][j][k]));
                    q_10.b += sin(k_R) * (sin(replicas[0].N_T_systems[sys_index].spin[i][j][k]) * cos(replicas[1].N_T_systems[sys_index].spin[i][j][k]));

                    q_01.a += cos(k_R) * (cos(replicas[0].N_T_systems[sys_index].spin[i][j][k]) * sin(replicas[1].N_T_systems[sys_index].spin[i][j][k]));
                    q_01.b += sin(k_R) * (cos(replicas[0].N_T_systems[sys_index].spin[i][j][k]) * sin(replicas[1].N_T_systems[sys_index].spin[i][j][k]));

                    q_11.a += cos(k_R) * (sin(replicas[0].N_T_systems[sys_index].spin[i][j][k]) * sin(replicas[1].N_T_systems[sys_index].spin[i][j][k]));
                    q_11.b += sin(k_R) * (sin(replicas[0].N_T_systems[sys_index].spin[i][j][k]) * sin(replicas[1].N_T_systems[sys_index].spin[i][j][k]));
                }
            }
        }

        return (q_00.a * q_00.a + q_00.b * q_00.b + q_10.a * q_10.a + q_10.b * q_10.b + q_01.a * q_01.a + q_01.b * q_01.b + q_11.a * q_11.a + q_11.b * q_11.b);
    }

    s_chi chi(int n_mean, double delta_theta, int sys_index)
    {
        s_chi chi_ensemble;
        double xi_L = 0;
        double k_min = 2 * pi / L;
        double k_0 = 0;
        chi_ensemble._0 = 0;
        chi_ensemble.min = 0;

        for (int i = 0; i < n_mean; i++)
        {
            avg_run(delta_theta * pi, 1, "E");

            chi_ensemble._0 += mod_sq_q(sys_index, k_0, 0, 0);
            chi_ensemble.min += mod_sq_q(sys_index, k_min, 0, 0);
        }

        xi_L = 1 / (2 * L * sin(k_min / 2)) * pow((chi_ensemble._0 / chi_ensemble.min) - 1, 0.5);

        return chi_ensemble;
    }

    double xi_L(int n_mean, double delta_theta, int sys_index)
    {

        double xi_L = 0;
        double k_min = 2 * pi / L;
        double k_0 = 0;
        double chi_0 = 0;
        double chi_min = 0;

        for (int i = 0; i < n_mean; i++)
        {
            avg_run(delta_theta * pi, 1, "E");

            chi_0 += mod_sq_q(sys_index, k_0, 0, 0);
            chi_min += mod_sq_q(sys_index, k_min, 0, 0);
        }

        xi_L = 1 / (2 * L * sin(k_min / 2)) * pow((chi_0 / chi_min) - 1, 0.5);

        return xi_L;
    }
};

// Questa classe mi è servita per calcolare le funzioni di correlazioni su diversi campioni
class montecarlo
{
public:
    int N_sample;
    int L;
    std::vector<replicas2> disorder;

    montecarlo(int N_sample_, int N_T, int L_, double T_min, double T_max)
        : N_sample(N_sample_), L(L_)
    {
        for (int i = 0; i < N_sample; i++)
        {
            disorder.push_back(replicas2(N_T, L, T_min, T_max, true));
        }
        for (int i = 0; i < N_sample; i++)
        {
            disorder[i].set_equal_J();
        }
    }

    void run(double theta, int mc_moves, const std::string filename)
    {
        for (int i = 0; i < N_sample; i++)
        {
            std::string updated_filename = filename + "_sample_" + std::to_string(i) + "_";
            disorder[i].run(theta, mc_moves, updated_filename);
        }
    }

    double xi_L_mean(int n_mean, double delta_theta, int sys_index)
    {
        s_chi temp;
        s_chi chi_quenched;
        double k_min = 2 * pi / L;
        double k_0 = 0;
        double xi_L;

        for (int i = 0; i < N_sample; i++)
        {
            cout << "Sample " << i << endl;
            temp = disorder[i].chi(n_mean, delta_theta, sys_index);
            chi_quenched.min += temp.min;
            chi_quenched._0 += temp._0;
        }

        xi_L = 1 / (2 * L * sin(k_min / 2)) * pow((chi_quenched._0 / chi_quenched.min) - 1, 0.5);

        return xi_L;
    };
};

#endif

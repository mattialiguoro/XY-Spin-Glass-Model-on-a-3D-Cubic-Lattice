#include <iostream>
#include <chrono>
#include "./montecarlo.hpp"

int main(void)
{
  rng.rseed();

  config config_data;
  config_data = uploadConfig();
  printConfig(config_data);



  // DUE REPLICHE DI N_T SISTEMI X N_SAMPLE CON J DIVERSIl

  // ###############********* RICORDARSI DI CANCELLARE I FILE DI ENERGIA DA PRECEDENTI SIMULAZIONI O SI AGGIUNGONO*********###############

  montecarlo mc(config_data.N_sample  , config_data.N_T, config_data.L, config_data.T_min, config_data.T_max);
  mc.run(config_data.fraction_of_pi * pi, config_data.mc_moves, "E");

  // ###############********* RICORDARSI DI CANCELLARE I FILE DI ENERGIA DA PRECEDENTI SIMULAZIONI O SI AGGIUNGONO*********###############

  // QUESTO SERVE A CALCOLARE LE LUNGHEZZE DI CORRELAZIONE SU PIù SAMPLE

  /*
  double xi_L = 0;
  double n_mean = 30;

  for(int i = 0; i < 10; i++)
  {
    xi_L = mc.xi_L_mean(n_mean, config_data.fraction_of_pi * pi,i);
    cout << "(T , xi_L) = (" << mc.disorder[0].replicas[0].temperatures[i] << ", " << xi_L << ")" << endl;
  }
  */

  // DUE REPLICHE DI N_T SISTEMI (stesso sample)

  /*
    replicas2 mc(config_data.N_T, config_data.L, config_data.T_min, config_data.T_max, true);
    mc.run(config_data.fraction_of_pi * pi, config_data.mc_moves, "E");
    mc.set_equal_J();
    xi_L = mc.xi_L(n_mean,config_data.fraction_of_pi * pi);
    cout << "xi_L: " << xi_L << endl;


      for (int i = 0; i < n_mean; i++) {
          mc.run(config_data.fraction_of_pi * pi, 1, "E");

          chi_0 += mc.mod_sq_q(0,k_0,0,0);
          chi_min +=mc.mod_sq_q(0,k_min,0,0);
      }

      q_0 = q_0/n_mean;
      q_min = q_min/n_mean;
      xi_L = 1 / (2*config_data.L * sin(k_min / 2)) * pow((chi_0 / chi_min) - 1, 0.5);

      cout << "\n xi_L: " << xi_L << endl;
  */

  // CODICE PER UN SISTEMA ALLA VOLTA

  /*
      sys s1(config_data.L);

      s1.initial_condition();
      s1.set_J();
      //s1.visualize_J();
      s1.set_sparse_J();
        std::cout << "L'energia iniziale del sistema e': " << s1.total_e() << std::endl;
      s1.run(pi * config_data.fraction_of_pi, config_data.T_min, config_data.mc_moves, "v_o_mc.txt");
      std::cout << "L'energia finale del sistema e': " << s1.total_e() << std::endl;
  */

  // s1.read_spin_from_file("spin_data.txt");

  // s1.metropolis(1/3 * pi, T, 50000, "v_o_mc.txt");
  // s1.visualize_spin();
  // s1.write_spin();

  // s1.write_J();

  return 0;
}

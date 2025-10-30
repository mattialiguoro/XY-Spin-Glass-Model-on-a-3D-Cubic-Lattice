#ifndef _RANDNUMGEN_
#define _RANDNUMGEN_
#include<chrono>
#include<random>
// use by default mersenne-twister rng engine
template <class ntype=double, class engine=std::mt19937_64>
class randnumgen {
  std::random_device rd;
  using myclock = std::chrono::high_resolution_clock;
  engine rng;
  using ud = std::uniform_real_distribution<ntype>;
  // we use a pointer here since if a range different from [0,1) is needed
  // we can use the constructor of uniform_real_distribution class 
  // to provide a different range (see randnumgen()) 
  // If a different range is not needed, the default constructor can be used
  // and one can avoid to use a pointer to uniform_real_distribution class
  ud* unidst;
public:
  // fixed seed
  void seed(int s)
    {
      std::seed_seq seq{s+1, s+2, s+2};
      rng.seed(seq);
    }
  // random seed
  void rseed(void)
    {
      // int operator()() {}
      // ticks since 1/1/1970
      unsigned int t=myclock::now().time_since_epoch().count(); 
      // ^ is bitwise XOR
      std::seed_seq seq{rd()^t,rd()^t, rd()^t};
      rng.seed(seq);
    }
  ntype ranf()
    {
      // ntype operator()(engine mt) {}
      // uniform_real_distribution overloads the parenthesis operator 
      return (*unidst)(rng); 
    } 
  randnumgen()
    {
      unidst = new ud(0.0,1.0);
      // equivalent way to set the range of random numbers:
      // unidst->param(std::uniform_real_distribution<double>::param_type(0.0, 1.0));
    }
  ~randnumgen()
    {
      delete unidst;
    }
  
  ntype ranf_int(long int a,long int b){
    return a + static_cast<int>( (*this).ranf() * (b - a));
  }

  ntype ranf_double(ntype a,ntype b){
    return a + (b-a)* (*this).ranf();
  }

};


// PUBLIC METHODS:
// rng.seed(<integer>) : fixed seeding 
// rng.rseed();        : random seeding 
// rng.ranf();         : generate a random number in [0.0,1.0)
// rng.ranf_int(a,b)       : generate integer number between [a,b]
// rng.ranf_double(a,b)       : generate number between [a,b]

randnumgen<> rng;
#endif

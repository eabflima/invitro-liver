#include <ip_compute.h>

int main(int argc, char* argv[]){
  // Initialize QUESO environment
#ifdef QUESO_HAS_MPI
  MPI_Init(&argc,&argv);
  UQ_FATAL_TEST_MACRO(argc < 2,QUESO::UQ_UNAVAILABLE_RANK,"main()","input file must be specified in command line as argv[1], just after executable argv[0]");
  QUESO::FullEnvironment* env =  new QUESO::FullEnvironment(MPI_COMM_WORLD,argv[1],"",NULL);
#else
  QUESO::FullEnvironment* env =  new QUESO::FullEnvironment(argv[1],"",NULL);
#endif
  
  {
    // Compute
    computeIP(*env, argc, argv);
  }
  
  // Finalize QUESO environment
  delete env;
#ifdef QUESO_HAS_MPI
  MPI_Finalize();
#endif
  return 0;
}

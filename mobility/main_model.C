#include "general_libraries.h"
#include "general_functions.h"
#include "tumor.h"
#include "main_model.h"

double Compute_CCC(EquationSystems& es,const unsigned int id,double& pearson,const unsigned int n_dofs);
void init_cond_st(EquationSystems & es, const std::string & system_name);
Number exact_value (const Point & p, const Parameters & es, const std::string & sys_n, const std::string & );

void main_code(MPI_Comm lib_comm, int argc, char** argv,vector<double>& TUMOR_MASS,vector<double> Parameters){
  //********** Initialize the library **********
  LibMeshInit init (argc, argv, lib_comm);

  /////----- Check the size of the parameters vector -----/////
  if(Parameters.size()!=6){
    std::cout << " Wrong number of parameters!" << std::endl;
    throw std::exception();
  }
  
  /////----- Read argumments -----/////
  GetPot input_file("options.in");
  GetPot input_data("inp_dat.in");
  
  const unsigned int r_steps     = input_file("r_steps",0);
  const unsigned int write_files = input_file("write_files",0);
  const unsigned int file_number = input_file("model_n",0);
  const unsigned int reinit      = input_file("reinit",0);
  const double time_step         = input_file("time_step", 0.05);
  const std::string mesh_name    = input_file("mesh_name", "test.msh");
  const std::string sol_name     = input_file("ic_res", "test.xda");
  const double l_day             = input_file("l_day", 10.0);
  const double init_time         = input_file("ic_time", 10.0);
  const unsigned int n_timesteps = (l_day*100.0-init_time*100.0)/(time_step*100.0);
  const unsigned int n_sol       = input_file("n_sol",0);
  const unsigned int n_days      = input_file("n_days",0);
  std::vector<double> save_days(n_days,0);
  for(unsigned int i = 0; i < n_days; i++)
    save_days[i] = input_file("s_days",0.0,i);
  std::vector<unsigned int> save_iter(n_days,0);
  for(unsigned int i = 0; i < n_days; i++)
    save_iter[i] = (save_days[i]*100.0-init_time*100.0)/(time_step*100.0);
 
  /////----- Create/Read mesh-----/////
  SerialMesh mesh(init.comm());
  mesh.read(mesh_name);
  EquationSystems eq_sys(mesh);
  MeshRefinement mesh_refinement (mesh);
  if(write_files){
    mesh.print_info();
  }
  TransientLinearImplicitSystem &tum = eq_sys.add_system<TransientLinearImplicitSystem>("S0_Tum");
  for(unsigned int f = 0; f < n_sol; f++){
    std::string name_f;
    name_f.append(input_file("file_names", "test.dat", f));
    std::string sys_f;
    sys_f.append(input_file("fsys_names", "test.dat", f));
    eq_sys.add_system<TransientLinearImplicitSystem>(sys_f);
    sys_f.clear();
    name_f.clear();
  }
  if(reinit){
    /////----- Add systems and variables -----/////
    eq_sys.read(sol_name, libMeshEnums::READ);
    tum.update();
    for(unsigned int f = 0; f < n_sol; f++){
      std::string name_f;
      name_f.append(input_file("file_names", "test.dat", f));
      std::string sys_f;
      sys_f.append(input_file("fsys_names", "test.dat", f));
      eq_sys.get_system(sys_f).update();
      sys_f.clear();
      name_f.clear();
    }
    tum.attach_assemble_function(Assemble_Tum);
    eq_sys.reinit();
  }
  else{
    std::vector<Point> pos;
    std::vector<double> sol;
    //=========================================================================== Tumor IC
    std::string ic_tum;
    ic_tum.append(input_file("ic_tum", "test.dat"));
    std::ifstream read;
    std::string line;
    double x,y,z,t0;
    read.open(ic_tum);
    if(!read.is_open())
      std::cout << "Error opening Tum.\n";
    while(std::getline(read,line)){
      std::istringstream(line) >> x >> y >> z >> t0;
      Point a(x,y,z);
      pos.push_back(a);
      sol.push_back(t0);
    }
    read.close();
    eq_sys.parameters.set<std::vector<Point>>("r_nodes") = pos;
    eq_sys.parameters.set<std::vector<double>>("IC_tum") = sol;
    tum.add_variable("v0_tum",FIRST);
    tum.add_variable("v1_che",FIRST);
    tum.attach_assemble_function(Assemble_Tum);
    tum.attach_init_function(Initial_Condition_Tum);
    //=========================================================================== System
    for(unsigned int f = 0; f < n_sol; f++){
      std::string name_f;
      name_f.append(input_file("file_names", "test.dat", f));
      std::string sol_f;
      sol_f.append(input_file(name_f, "test.dat"));
      read.open(sol_f);
      if(!read.is_open())
	std::cout << "Error opening Sol.\n";
      for(unsigned i=0; i<sol.size(); i++){
	std::getline(read,line);
	std::istringstream(line) >> x >> y >> z >> t0;
	sol[i] = t0;
      }
      std::string sys_f;
      sys_f.append(input_file("fsys_names", "test.dat", f));
      eq_sys.parameters.set<std::vector<double>>(sys_f) = sol;
      eq_sys.get_system(sys_f).add_variable(name_f,FIRST);
      eq_sys.get_system(sys_f).attach_init_function(init_cond_st);
      read.close();
      sol_f.clear();
      sys_f.clear();
      name_f.clear();
    }
    eq_sys.init();
  }
  if(write_files){
    eq_sys.print_info();
  }
  /////----- Defining parameters -----/////
  char buffer[400];
  unsigned int t_step = 0;
  tum.time = init_time;
  const unsigned int n_nonlinear_steps = 50;//100;
  const Real nonlinear_tolerance       = 1.e-8;
  eq_sys.parameters.set<unsigned int>("linear solver maximum iterations") = 50;//250;
  eq_sys.parameters.set<Real>("time_step") = time_step;
  eq_sys.parameters.set<Real>("m_v")       = Parameters[0]; // M_V
  eq_sys.parameters.set<Real>("E_t")       = Parameters[1]; // \bar{E}_T
  eq_sys.parameters.set<Real>("eps")       = Parameters[2]; // \eps_T
  eq_sys.parameters.set<Real>("l_p")       = Parameters[3]; // l_p
  eq_sys.parameters.set<Real>("l_k")       = Parameters[4]; // l_k
  eq_sys.parameters.set<Real>("l_a")       = Parameters[5]; // l_a
  for(unsigned int i=0; i<r_steps; i++){
    mesh_refinement.uniformly_refine(1);
    eq_sys.reinit();
  }
  /////----- Saving initial condition -----/////
  if(write_files){
    sprintf(buffer, "exo%d-%05d.e",file_number,t_step);
    ExodusII_IO(mesh).write_equation_systems(buffer, eq_sys);
  }
  if(!reinit){
    eq_sys.write(sol_name, libMeshEnums::WRITE);
  }
  if(write_files){
    double CCC_v,PCC_v;
    CCC_v = Compute_CCC(eq_sys,n_sol+1,PCC_v,n_sol+2);
    cout << "Time = " << init_time << ", CCC = " << CCC_v << ", PCC = " << PCC_v << endl;
  }
  /////----- Solving the problem -----/////
  UniquePtr<NumericVector<Number> > last_nonlinear_soln_tum (tum.solution->clone());
  unsigned int n_id = 1;
  do{
    /////----- Prepare time step-----/////
    t_step++;
    tum.time = init_time+t_step*time_step;
    eq_sys.parameters.set<Real>("time") = tum.time;
    *tum.old_local_solution = *tum.current_local_solution;
    /////----- Nonlinear iteration loop -----/////
    const Real initial_linear_solver_tol = 1.e-8;
    eq_sys.parameters.set<Real> ("linear solver tolerance") = initial_linear_solver_tol;
    for (unsigned int l=0; l<n_nonlinear_steps; ++l){
      /////----- Solve systems -----/////
      last_nonlinear_soln_tum->zero();
      last_nonlinear_soln_tum->add(*tum.solution);
      tum.solve();
      last_nonlinear_soln_tum->add(-1., *tum.solution);
      last_nonlinear_soln_tum->close();
      
      /////----- Nonlinear iteration error -----/////
      const Real norm_delta = last_nonlinear_soln_tum->linfty_norm();
      const unsigned int n_linear_iterations = tum.n_linear_iterations();
      const Real final_linear_residual = tum.final_linear_residual();
      if(write_files){std::cout << "   Linear converged at step: " << n_linear_iterations << ", residual: " << final_linear_residual << " Nonlinear convergence: ||u - u_old|| = " << norm_delta << std::endl;}
      if(norm_delta < nonlinear_tolerance){
	if(write_files){std::cout << "Nonlinear converged at step: " << l << std::endl;}
	break;
      }
    }
    /////----- Post-processing -----/////
    int need_to_print = 0;
    for(unsigned int i = 0; i < n_days; i++)
      if(save_iter[i] == t_step){
	need_to_print = 1;
	n_id++;	
      }
    if(write_files){
      need_to_print=1;
      ExodusII_IO exo(mesh);
      exo.append(true);
      exo.write_timestep (buffer, eq_sys, t_step+1, tum.time);
    }
    if(need_to_print){
      double PCC_v;
      Compute_CCC(eq_sys,n_id,PCC_v,n_sol+2);
      TUMOR_MASS.push_back(PCC_v);
    }
  }while(t_step<n_timesteps);
}

/* Set initial condition */
void init_cond_st(EquationSystems & es, const std::string & system_name){
  cout << system_name << endl;
  System & system = es.get_system<System>(system_name);
  system.project_solution(exact_value, libmesh_nullptr, es.parameters);
}

/* Initial condition */
Number exact_value (const Point & p, const Parameters & es, const std::string & sys_n, const std::string & ){
  std::vector<Point> r_nodes = es.get<std::vector<Point>>("r_nodes");
  std::vector<double> IC_vals = es.get<std::vector<double>>(sys_n);
  for(unsigned int i=0; i<r_nodes.size(); i++){
    if(r_nodes[i]==p){
      return IC_vals[i];
    }
  }
  return 0.;
}

/* Calculate the Concordance Correlation Coefficient */
double Compute_CCC(EquationSystems& es,const unsigned int id,double& pearson,const unsigned int n_dofs){
  std::vector<Number> soln;
  es.build_solution_vector(soln);
  const unsigned int model_id = 0;
  const unsigned int exper_id = id;
  const unsigned int vector_n = soln.size()/n_dofs;
  double mean_x, mean_y;
  double var_x, var_y, var_xy;
  var_x = var_y = var_xy = 0.0;
  mean_x = mean_y = 0.0;
  for(unsigned int i = 0; i < vector_n; i++){
    mean_x += soln[i*(n_dofs-1)+i+model_id];
    mean_y += soln[i*(n_dofs-1)+i+exper_id];
  }
  mean_x = mean_x/vector_n;
  mean_y = mean_y/vector_n;
  for(unsigned int i = 0; i < vector_n; i++){
    var_x += std::pow(soln[i*(n_dofs-1)+i+model_id]-mean_x,2);
    var_y += std::pow(soln[i*(n_dofs-1)+i+exper_id]-mean_y,2);
    var_xy += (soln[i*(n_dofs-1)+i+model_id]-mean_x)*(soln[i*(n_dofs-1)+i+exper_id]-mean_y);
  }
  var_x = var_x/vector_n;
  var_y = var_y/vector_n;
  var_xy = var_xy/vector_n;
  
  pearson = var_xy/(std::sqrt(var_x)*std::sqrt(var_y));
  return 2.0*var_xy/(var_x+var_y+std::pow(mean_x-mean_y,2));
}

#include "general_libraries.h"
#include "general_functions.h"

/* Save output to file */
void Sprint(const EquationSystems& es,const MeshBase& mesh,string s,int file_number,int t){
  libmesh_assert_equal_to(libMesh::processor_id(), 0);
  Point p;
  std::vector<Number> soln;
  std::vector<std::string> names;
  es.build_variable_names(names);
  es.build_solution_vector(soln);
  FILE *arq; 
  const char *c = s.c_str();
  char n[100],name[200];
  sprintf(n,"%d-%05d.txt",file_number,t); 
  strcpy(name,c);
  strcat(name,n);
  arq = fopen(name,"w");
  for(unsigned int i=0;i<mesh.n_nodes();i++){
    p = mesh.point(i);
    for(unsigned int j=0;j<mesh.mesh_dimension();j++){
      fprintf(arq,"\t%e",p(j));
    }
    const unsigned int n_vars = names.size();
    for(unsigned int c=0;c<n_vars;c++){
      fprintf(arq,"\t%e",soln[i*n_vars + c]);
    }
    fprintf(arq,"\n");
  }
  fclose(arq);
}

/* Print mesh points */
void Pmesh(const MeshBase& mesh){
  Point p;
  libmesh_assert_equal_to(libMesh::processor_id(), 0);
  FILE *arq; 
  arq = fopen("mesh.txt","w");
  for(unsigned int i=0;i<mesh.n_nodes();i++){
    p = mesh.point(i);
    for(unsigned int j=0;j<mesh.mesh_dimension();j++){
      fprintf(arq,"\t%e",p(j));
    }
    fprintf(arq,"\n");
  }
  fclose(arq);
}

/* Return the Square of numbers and ignores values lower than 0 and higher than 1 */
double Square(double x){
  if(x<0){return 0.;}
  else if(x>1){return 1.;}
  else return x*x;
}

/* Ignores values lower than 0 and higher than 1 */
double Linear(double x){
  if(x<0){return 0.;}
  else if(x>1){return 1.;}
  else return x;
}

/* Regularized Heaviside function */
double Heaviside(double x){
  if(fabs(x)>=0.064){
    return (x<0) ? 0. : 1.;
  }
  else{
    return (1./(2.*0.064))*(x+0.064);
  }
}

/* Return Modulus of a vector */
double Modulus(double xg,double yg){
  return sqrt(pow(xg,2)+pow(yg,2));
}

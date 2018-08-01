#include "general_libraries.h"
#include "general_functions.h"
#include "tumor.h"

/* Assemble tumor system */
void Assemble_Tum(EquationSystems& es,const std::string& libmesh_dbg_var(system_name)){
  libmesh_assert_equal_to(system_name,"S0_Tum");
  const MeshBase& mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  TransientLinearImplicitSystem & tum = es.get_system<TransientLinearImplicitSystem>("S0_Tum");
  std::vector<unsigned int> v_tum(2);
  v_tum[0] = tum.variable_number("v0_tum");
  v_tum[1] = tum.variable_number("v1_che");
  const DofMap& tum_map = tum.get_dof_map();
  FEType fe_type = tum.variable_type(0);
  UniquePtr<FEBase> fe(FEBase::build(dim,fe_type));
  QGauss qrule(dim,fe_type.default_quadrature_order());
  fe->attach_quadrature_rule(&qrule);
  const std::vector<Real>& JxW = fe->get_JxW();
  const std::vector<std::vector<Real> >& phi = fe->get_phi();
  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();
  const Real time_size = es.parameters.get<Real>("time_step");
  const Real m_v       = es.parameters.get<Real>("m_v");
  const Real E_t       = es.parameters.get<Real>("E_t");
  const Real eps       = es.parameters.get<Real>("eps");
  const Real l_p       = es.parameters.get<Real>("l_p");
  const Real l_k       = es.parameters.get<Real>("l_k");
  const Real l_a       = es.parameters.get<Real>("l_a");
  DenseMatrix<Number> Ke;
  DenseSubMatrix<Number> Ke_var[2][2] =
    {
      {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)},
      {DenseSubMatrix<Number>(Ke), DenseSubMatrix<Number>(Ke)}
    };
  DenseVector<Number> Fe;
  DenseSubVector<Number> Fe_var[2] =
    {DenseSubVector<Number>(Fe),
     DenseSubVector<Number>(Fe)};
  std::vector<unsigned int> dof_indices;
  std::vector< std::vector<dof_id_type> > dof_indices_var(2);
  std::vector<unsigned int> dof_indices_nut;
  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end(); 
  for(;el != end_el;++el){
    const Elem* elem = *el;
    tum_map.dof_indices(elem,dof_indices);
    for (unsigned int var=0; var<2; var++)
      tum_map.dof_indices (elem, dof_indices_var[var], v_tum[var]);
    const unsigned int n_dofs = dof_indices.size();
    const unsigned int n_var_dofs = dof_indices_var[0].size();
    fe->reinit (elem);
    Ke.resize (n_dofs, n_dofs);
    for (unsigned int var_i=0; var_i<2; var_i++)
      for (unsigned int var_j=0; var_j<2; var_j++)
	Ke_var[var_i][var_j].reposition (var_i*n_var_dofs, var_j*n_var_dofs, n_var_dofs, n_var_dofs);
    Fe.resize (n_dofs);
    for (unsigned int var=0; var<2; var++)
      Fe_var[var].reposition (var*n_var_dofs, n_var_dofs);
    for(unsigned int qp=0; qp<qrule.n_points(); qp++){
      Number tum_cur = 0., tum_old = 0.;
      for(unsigned int l=0; l<phi.size(); l++){
	tum_old += phi[l][qp]*tum.old_solution(dof_indices_var[0][l]);
	tum_cur += phi[l][qp]*tum.current_solution(dof_indices_var[0][l]);
      }
      double m_t = m_v;//*pow(Linear(tum_cur),2);//*pow(1.0-Linear(tum_cur),2);
      for(unsigned int i=0; i<phi.size(); i++){
	//-- Tumor --//
	Fe_var[0](i) += JxW[qp]*(tum_old
				 )*phi[i][qp];
	//-- Chemical Potential --//
	Fe_var[1](i) += JxW[qp]*2.0*E_t*tum_old*(2.0*std::pow(tum_old,2)-3.0*tum_old)*phi[i][qp];
	for (unsigned int j=0; j<phi.size(); j++){
	  //-- Ka Tumor --//
	  Ke_var[0][0](i,j) += JxW[qp]*(1.0+
					time_size*(l_a-l_p*(1.0-(tum_cur/l_k))))*phi[j][qp]*phi[i][qp];
	  Ke_var[0][1](i,j) += JxW[qp]*time_size*m_t*dphi[j][qp]*dphi[i][qp];
	  //-- Kb Chemical_tumor --//
	  Ke_var[1][0](i,j) -= JxW[qp]*2.0*E_t*phi[j][qp]*phi[i][qp];
	  Ke_var[1][0](i,j) -= JxW[qp]*std::pow(eps,2)*dphi[j][qp]*dphi[i][qp];
	  Ke_var[1][1](i,j) += JxW[qp]*phi[j][qp]*phi[i][qp];
	}
      }
    }
    tum_map.heterogenously_constrain_element_matrix_and_vector(Ke,Fe,dof_indices);
    tum.matrix->add_matrix(Ke, dof_indices);
    tum.rhs->add_vector(Fe, dof_indices);
  }
  return;
}

/* Calculate the Mass of the Tumor */
void Mass_Tum(EquationSystems& es,const std::string& libmesh_dbg_var(system_name),double &value_mass){
  libmesh_assert_equal_to(system_name,"S0_Tum");
  const MeshBase& mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  Number tumor_mass = 0.; 
  TransientLinearImplicitSystem & tum = es.get_system<TransientLinearImplicitSystem>("S0_Tum");
  const unsigned int v_tum = tum.variable_number("v0_tum");
  FEType fe_type = tum.variable_type(0);
  UniquePtr<FEBase> fe(FEBase::build(dim,fe_type));
  QGauss qrule(dim,fe_type.default_quadrature_order());
  fe->attach_quadrature_rule(&qrule);
  const std::vector<Real>& JxW = fe->get_JxW();
  const std::vector<std::vector<Real> >& phi = fe->get_phi();
  const DofMap& tum_map = tum.get_dof_map();
  std::vector<unsigned int> dof_indices;
  std::vector<unsigned int> dof_indices_tum;
  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end(); 
  for(;el != end_el;++el){
    const Elem* elem = *el;
    tum_map.dof_indices(elem,dof_indices);
    tum_map.dof_indices(elem,dof_indices_tum,v_tum);
    fe->reinit(elem);
    for(unsigned int qp=0; qp<qrule.n_points(); qp++){
      Number partial_tumor_mass = 0.;
      for(unsigned int l=0; l<phi.size(); l++){
	double tumor_val = tum.current_solution(dof_indices_tum[l]);
	if(tumor_val > 1.0)
	  tumor_val = 1.0;
	else if(tumor_val < 0.0)
	  tumor_val = 0.0;
	partial_tumor_mass += JxW[qp]*phi[l][qp]*tumor_val;
      }
      tumor_mass += partial_tumor_mass;
    }
  }
  value_mass = tumor_mass;
  return;
}

/* Initial condition */
Number exact_value_Tum(const Point& p,const Parameters& es,const std::string& ,const std::string& var_name){
  std::vector<Point> r_nodes = es.get<std::vector<Point>>("r_nodes");
  std::vector<double> IC_vals = es.get<std::vector<double>>("IC_tum");
  if (var_name == "v0_tum"){
    for(unsigned int i=0; i<r_nodes.size(); i++){
      if(r_nodes[i]==p){
	return IC_vals[i];
      }
    }
  }
  return 0.;
}

/* Set initial condition */
void Initial_Condition_Tum(EquationSystems& es,const std::string& libmesh_dbg_var(system_name)){
  libmesh_assert_equal_to(system_name,"S0_Tum");
  TransientLinearImplicitSystem & sys = es.get_system<TransientLinearImplicitSystem>("S0_Tum");
  es.parameters.set<Real>("time") = sys.time = 0;
  sys.project_solution(exact_value_Tum,NULL,es.parameters);
}

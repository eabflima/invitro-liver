/* Assemble tumor system */
void Assemble_Tum(EquationSystems& es,const std::string& libmesh_dbg_var(system_name));

/* Calculate the Mass of the Tumor */
void Mass_Tum(EquationSystems& es,const std::string& libmesh_dbg_var(system_name),double &value_mass);

/* Initial condition */
Number exact_value_Tum(const Point& p,const Parameters& es,const std::string& ,const std::string& var_name);

/* Set initial condition */
void Initial_Condition_Tum(EquationSystems& es,const std::string& libmesh_dbg_var(system_name));

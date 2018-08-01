/////----- General Functions -----/////

/* Save output to file */
void Sprint(const EquationSystems& es,const MeshBase& mesh,string s,int file_number,int t);

/* Print mesh points */
void Pmesh(const MeshBase& mesh);

/* Return the Square of numbers and ignores values lower than 0 and higher than 1 */
double Square(double x);

/* Ignores values lower than 0 and higher than 1 */
double Linear(double x);

/* Regularized Heaviside function */
double Heaviside(double x);

/* Return Modulus of a vector */
double Modulus(double xg,double yg);

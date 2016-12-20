/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2007 - 2015 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------

 *
 * Author: David Neckels, Boulder, Colorado, 2007, 2008
 */


// @sect3{Include files}

// First a standard set of deal.II includes. Nothing special to comment on
// here:
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/function_parser.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/std_cxx11/array.h>

#include <deal.II/lac/vector.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_in.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/solution_transfer.h>

#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>


// Sacado is the automatic differentiation package within Trilinos, which is
// used to find the Jacobian for a fully implicit Newton iteration:
// Trilinos::Sacado (at least until version 11.10.2) package will trigger
// warnings when compiling this file. Since we are not responsible for this,
// we just suppress the warning by wrapping the <code>#include</code>
// directive into a pair of macros that simply suppress these warnings:
DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#include <Sacado.hpp>
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS

// And this again is C++:
#include <iostream>
#include <fstream>
#include <vector>
#include <memory>
#include <cmath>


namespace Step33
{
  using namespace dealii;


  // @sect3{Euler equation specifics}

  // Here we define the flux function for this particular system of conservation laws, as well as pretty much everything else that's specific
  // to the Euler equations for gas dynamics, for reasons discussed in the introduction. We group all this into a structure that defines everything
  // that has to do with the flux. All members of this structure are static, i.e. the structure has no actual state specified by instance member
  // variables. The better way to do this, rather than a structure with all static members would be to use a namespace -- but namespaces can't be
  // templatized and we want some of the member variables of the structure to depend on the space dimension, which we in our usual way introduce using
  // a template parameter.
  template <int dim>
  struct PTsolver
  {
	  // PAGE 8
	  // COMPONENT DESCRIPTION AND NAMES (VECTOR/SCALAR)
    static const unsigned int n_components             = 2*dim + 4;
    static const unsigned int VES_component      	   = 0;
    static const unsigned int pf_component             = 1;
    static const unsigned int vf_first_component       = 2;
    static const unsigned int vr_first_component       = 2+dim;
    static const unsigned int poro_component           = 2+2*dim;
    static const unsigned int temp_component           = 2+2*dim+1;

    static
    std::vector<std::string>
    component_names ()
    {
    	std::vector<std::string> names (1, "VES-rock-fluid_pressure");
    	names.push_back ("fluid_pressure");
    	std::vector<std::string> vf (dim, "fluid_velocity");
    	std::vector<std::string> vr (dim, "rock_velocity");
    	names.insert(names.end(), vf.begin(), vf.end());
    	names.insert(names.end(), vr.begin(), vr.end());
    	names.push_back ("porosity");
    	names.push_back ("temperature");

      return names;
    }

    static
    std::vector<DataComponentInterpretation::DataComponentInterpretation>
    component_interpretation ()
    {

    	std::vector<DataComponentInterpretation::DataComponentInterpretation>
    	data_component_interpretation (1, DataComponentInterpretation::component_is_scalar); //p_r-p_f - initial data interpretation
        data_component_interpretation.push_back (DataComponentInterpretation::component_is_scalar); // add another scalar for p_f
    	std::vector<DataComponentInterpretation::DataComponentInterpretation> vfDCI (dim, DataComponentInterpretation::component_is_part_of_vector);
    	std::vector<DataComponentInterpretation::DataComponentInterpretation> vrDCI (dim, DataComponentInterpretation::component_is_part_of_vector);
    	data_component_interpretation.insert(data_component_interpretation.end(), vfDCI.begin(), vfDCI.end()); // two vector push backs for dim-dimensional amount of data interpretation
    	data_component_interpretation.insert(data_component_interpretation.end(), vrDCI.begin(), vrDCI.end());
    	data_component_interpretation.push_back (DataComponentInterpretation::component_is_scalar);
    	data_component_interpretation.push_back (DataComponentInterpretation::component_is_scalar);

      return data_component_interpretation;
    }


    // PARAMETER DECLATAION AND DEFINE FUNCTIONS
    static const double lambda;
    static const double gamma;
    static const double shear_modulus;
    static const double nusselt;
    static const double mu_0;
    static const double T_0;
    static const double rho_r;
    static const double E;
    static const double R;
    static const double c_f;
    static const double C_tilde; // C*Re*Ri*phi_0¹.1
    static const double phi_0;

    //RHO_F
    template <typename InputVector>
    static
    typename InputVector::value_type
    compute_rho_f (const InputVector &W)
    {
      return -9.9e-6*W[temp_component]*W[temp_component]*W[temp_component] + 0.002712*W[temp_component]*W[temp_component] -
    		  	  	  0.8176*W[temp_component] + 1032.0;
    }
    //C_R
    template <typename InputVector>
       static
       typename InputVector::value_type
       compute_c_r (const InputVector &W)
       {
         return -0.001361*W[temp_component]*W[temp_component] + 1.478*W[temp_component] +737.0;
       }
    //KAPPA
    template <typename InputVector>
         static
         typename InputVector::value_type
         compute_kappa (const InputVector &W)
         {
           typename InputVector::value_type kappar = 0.0;
           typename InputVector::value_type kappaf = 0.0;
           typename InputVector::value_type kappa = 0.0;

           kappaf = -6.0e-6*W[temp_component]*W[temp_component] + 0.00162*W[temp_component] + 0.5643;
           kappar = 5.0e-6*W[temp_component]*W[temp_component] - 0.00585*W[temp_component] + 3.515;
           kappa = pow(kappaf,W[poro_component]) * pow(kappar,(1.0-W[poro_component]));

           return kappa;
         }
    //MU_R
    template <typename InputVector>
        static
        typename InputVector::value_type
        compute_mu_r (const InputVector &W)
        {
          return mu_0*(exp((E/R) * (1.0/ W[temp_component] - 1.0/T_0)));
        }
    //K, PERMEABILITY
    template <typename InputVector>
    	static
		typename InputVector::value_type
		compute_perm (const InputVector &W)
    	{
    		return pow(( W[poro_component]/phi_0 ),8.0);
    	}
    //CHEMICAL COMPACTION COEFFICIENT, EVERYTHING UP TO (PR-PF)
    template <typename InputVector>
    	static
		typename InputVector::value_type
		compute_chem_coeff (const InputVector &W)
    	{
    		return ( C_tilde * pow(W[poro_component],0.1) * pow((1.0-W[poro_component]),2.0) * (1.0/compute_mu_r(W)) );

    	}
    // mech compac coeff where = gamma*exp(-gamma*ves))
    template <typename InputVector>
    	static
		typename InputVector::value_type
		compute_mech_coeff (const InputVector &W)
    	{
    		return (gamma * exp(-gamma*W[VES_component]));
    	}
    //BULK DENSITY USED IN EQUATIONS 4 (FORCING)
    template <typename InputVector>
        	static
    		typename InputVector::value_type
    		compute_bulk_density (const InputVector &W)
        	{
        		return ( (1.0-W[poro_component])*rho_r + W[poro_component]*compute_rho_f(W) );

        	}
    //BULK DENSITY WITH RELATIVE HEAT CAPACITY for rock and fluid USED IN TIME FOR TEMP EQUATION: BE CAREFUL - THIS IS RETURNING A NUMBER
    template <typename InputVector>
            	static
        		typename InputVector::value_type
        		compute_bulk_c_r (const InputVector &W)
            	{
            		return ( (1.0-W[poro_component])*rho_r*compute_c_r(W) ); // (1-phi)*rho_r*c_r

            	}
    template <typename InputVector>
              	static
          		typename InputVector::value_type
          		compute_bulk_c_f (const InputVector &W)
              	{
              		return ( W[poro_component]*compute_rho_f(W)*c_f ); // phi*rho_f*c_f

              	}



    template <typename InputVector>
    static
    void compute_flux_matrix (const InputVector &W,
                              std_cxx11::array <std_cxx11::array
                              <typename InputVector::value_type, dim>, // a row of (dim columns) and n_components rows of them stacked below. called 'flux'
                              PTsolver<dim>::n_components > &flux)
    {

      const typename InputVector::value_type bulk_c_r = compute_bulk_c_r(W);
      const typename InputVector::value_type bulk_c_f = compute_bulk_c_f(W);
      const typename InputVector::value_type rho_f = compute_rho_f(W);

      // equations for first two continuity equations

      for (unsigned int d=0; d<dim; ++d)
        flux[VES_component][d] = (1.0-W[poro_component])*rho_r*W[vr_first_component+d];

      for (unsigned int d=0; d<dim; ++d)
        flux[pf_component][d] = W[poro_component]*rho_f*W[vf_first_component+d];

      // Momentum equations portion

      for (unsigned int d=0; d<dim; ++d)
        {
          for (unsigned int e=vf_first_component; e<dim; ++e)
            flux[vf_first_component+d][e] = 0.0; // initialise all values to 0.

          flux[vf_first_component+d][d] += W[pf_component]; // diagonals equal to p_f

        }
      for (unsigned int d=0; d<dim; ++d)
             {
               for (unsigned int e=vr_first_component; e<dim; ++e)
                 flux[vr_first_component+d][e] = 0.0; // initialise all values to 0.

               flux[vr_first_component+d][d] += (1.0-W[poro_component])*W[VES_component] + W[pf_component];

             }

      // components for porosity and temperature equations

      for (unsigned int d=0; d<dim; ++d)
             flux[poro_component][d] = 0.0; // no straight forward div terms in porosity equation


      for (unsigned int d=0; d<dim; ++d)
             flux[temp_component][d] = ( bulk_c_f*W[vf_first_component+d] + bulk_c_r*W[vr_first_component+d] ) * W[temp_component];


    }


    template <typename InputVector>
    static
    void numerical_normal_flux (const Tensor<1,dim>                &normal,
                                const InputVector                  &Wplus,
                                const InputVector                  &Wminus,
                                const double                        alpha,
                                std_cxx11::array
                                <typename InputVector::value_type, n_components>
                                &normal_flux)
    {
      std_cxx11::array
      <std_cxx11::array <typename InputVector::value_type, dim>,
      PTsolver<dim>::n_components > iflux, oflux;

      compute_flux_matrix (Wplus, iflux);
      compute_flux_matrix (Wminus, oflux);

      for (unsigned int di=0; di<n_components; ++di)
        {
          normal_flux[di] = 0;
          for (unsigned int d=0; d<dim; ++d)
            normal_flux[di] += 0.5*(iflux[di][d] + oflux[di][d]) * normal[d];

          normal_flux[di] += 0.5*alpha*(Wplus[di] - Wminus[di]);
        }
    }



    template <typename InputVector>
    static
    void compute_forcing_vector (const InputVector &W, // NOW THE MATRIX IS CALLED FORCING
                                 std_cxx11::array
                                 <typename InputVector::value_type, n_components> // JUST N COMPONENTS
                                 &forcing)
    {
    	const typename InputVector::value_type permeability = compute_perm(W);
    	const typename InputVector::value_type chem_coeff = compute_chem_coeff(W);
    	const typename InputVector::value_type rho_f = compute_rho_f(W);
    	const typename InputVector::value_type bulk_density = compute_bulk_density(W);

      for (unsigned int c=0; c < n_components; ++c)
        switch (c) // whatever the dimension, the z direction and porosity component are filled and everything else is 0.0
          {
          case vf_first_component+dim-1:
        	forcing[c] = 1./(lambda*permeability) * (W[vf_first_component+dim-1] - W[vr_first_component+dim-1]) - rho_f;
        	break;
          case vr_first_component+dim-1:
            forcing[c] = bulk_density;
            break;
          case poro_component:
            forcing[c] = chem_coeff * W[VES_component];
            break;
          default:
            forcing[c] = 0; // ZERO FOR EVERYTHING ELSE
        }

      switch (dim) // add in teh extra components that come up in darcyś for 2d and 3d. Note we dont need to do anything with bulk density here anymore for rock mom eqn
      	  {
      	  case 2:
      		  forcing[vf_first_component] = 1./(lambda*permeability) * (W[vf_first_component] - W[vr_first_component]);
      		break;
      	  case 3:
      		  forcing[vf_first_component] = 1./(lambda*permeability) * (W[vf_first_component] - W[vr_first_component]);
      		  forcing[vf_first_component+1] = 1./(lambda*permeability) * (W[vf_first_component+1] - W[vr_first_component+1]);
      		  break;
      	  }


    }

    template <typename InputVector>
     static
     void compute_coeff_vector (const InputVector &W,
                                  std_cxx11::array
                                  <typename InputVector::value_type, n_components>
                                  &coeff)
     {
    	const typename InputVector::value_type kappa_total = compute_kappa(W);

       for (unsigned int c=0; c<n_components; ++c)
         switch (c)
           {
           case vr_first_component+dim-1:
             coeff[c] = -(1.0-W[poro_component])*2.0*shear_modulus;
             break;
           case temp_component:
             coeff[c] = -1.0/(nusselt*phi_0) * kappa_total;
             break;
           default:
             coeff[c] = 0;
           }
     }

    template <typename InputVector>
    static
	void compute_extraporo1_vector (const InputVector &W,
									std_cxx11::array
			                        <typename InputVector::value_type, n_components>
    								&extraporo1)
    {

    		for (unsigned int c=0; c< n_components; ++c)
    			switch(c)
    			{
    			case poro_component:
    				extraporo1[c] = W[vr_first_component+dim-1];
    				break;
    			default:
    				extraporo1[c] = 0;
    			}
    }

    template <typename InputVector>
    static
	void compute_extraporo2_vector (const InputVector &W,
									std_cxx11::array
			                        <typename InputVector::value_type, n_components>
    								&extraporo2)
    {
    		const typename InputVector::value_type mech_coeff = compute_mech_coeff(W);

    		for (unsigned int c=0; c< n_components; ++c)
    			switch(c)
    			{
    			case poro_component:
    				extraporo2[c] = mech_coeff*W[vr_first_component+dim-1];
    				break;
    			default:
    				extraporo2[c] = 0;
    			}
    }

    template <typename InputVector>
     static
 	void compute_time_vector (const InputVector &W,
 									std_cxx11::array
 			                        <typename InputVector::value_type, n_components>
     								&time)
     {
    	const typename InputVector::value_type rho_f = compute_rho_f(W);
    	const typename InputVector::value_type bulk_c_r = compute_bulk_c_r(W);
    	const typename InputVector::value_type bulk_c_f = compute_bulk_c_f(W);

     		for (unsigned int c=0; c< n_components; ++c)
     			switch(c)
     			{
     			case VES_component:
     				time[c] = rho_r*(1.0-W[poro_component]);
     				break;
     			case pf_component:
     				time[c] = rho_f*W[poro_component];
     				break;
     			case poro_component:
     				time[c] = W[poro_component] - exp(-gamma* (W[VES_component] ));
     				break;
     			case temp_component:
     				time[c] = (bulk_c_r + bulk_c_f)*W[temp_component];
     				break;
     			default:
     				time[c] = 0;
     			}
     }


    };




  template <int dim>
  class PTsolver
  {
  public:
		  PTsolver (const char *input_filename); // For input of parameters
		  void run ();

  private:
		 void setup_system ();

		 void assemble_system (); // 1 of 3: main loop over all cells

		 void assemble_cell_term (const FEValues<dim>             &fe_v,    		// 2 of 3: calls other two function for integrals over cells
								  const std::vector<types::global_dof_index> &dofs);

		 void assemble_face_term (const unsigned int               face_no,			// 3 of 3: and over faces
								  const FEFaceValuesBase<dim>     &fe_v,
								  const FEFaceValuesBase<dim>     &fe_v_neighbor,
								  const std::vector<types::global_dof_index> &dofs,
								  const std::vector<types::global_dof_index> &dofs_neighbor,
								  const bool                       external_face,
								  const unsigned int               boundary_id,
								  const double                     face_diameter);

		 std::pair<unsigned int, double> solve (Vector<double> &solution);

		 void compute_refinement_indicators (Vector<double> &indicator) const;
		 void refine_grid (const Vector<double> &indicator);

		 void output_results () const;

		 Triangulation<dim>   triangulation;
		 const MappingQ1<dim> mapping;

		 const FESystem<dim>  fe;
		 DoFHandler<dim>      dof_handler;

		 const QGauss<dim>    quadrature;
		 const QGauss<dim-1>  face_quadrature;

		 Vector<double>       old_solution;			// solution of previous time step
		 Vector<double>       current_solution;		// guess of current solution (NOT NECESSARILY CONVERGED)
		 Vector<double>       predictor;			// predictor for the solution at next time step

		 Vector<double>       right_hand_side;

		 TrilinosWrappers::SparseMatrix system_matrix;

		 Parameters::AllParameters<dim>  parameters;
		 ConditionalOStream              verbose_cout;
  };

// PTsolver::PTsolver here

  template <int dim>
  void
  ConservationLaw<dim>::
  assemble_cell_term (const FEValues<dim>             &fe_v,
                      const std::vector<types::global_dof_index> &dof_indices)
  {
    const unsigned int dofs_per_cell = fe_v.dofs_per_cell;
    const unsigned int n_q_points    = fe_v.n_quadrature_points;

    Table<2,Sacado::Fad::DFad<double> >
    W (n_q_points, PTsolver<dim>::n_components);

    Table<2,double>
    W_old (n_q_points, PTsolver<dim>::n_components);

    Table<3,Sacado::Fad::DFad<double> >
    grad_W (n_q_points, PTsolver<dim>::n_components, dim);

    Table<3,double>
    grad_W_old(n_q_points, PTsolver<dim>::n_components, dim);

    std::vector<double> residual_derivatives (dofs_per_cell);

    std::vector<Sacado::Fad::DFad<double> > independent_local_dof_values(dofs_per_cell);
    for (unsigned int i=0; i<dofs_per_cell; ++i)
      independent_local_dof_values[i] = current_solution(dof_indices[i]);

    for (unsigned int i=0; i<dofs_per_cell; ++i)
    	independent_local_dof_values[i].diff (i, dofs_per_cell);

    ///////////////////////////////////////////////////////////////////////////
    // initialise the values of W etc then actually calculate

    for (unsigned int q=0; q<n_q_points; ++q)
      for (unsigned int c=0; c<PTsolver<dim>::n_components; ++c)
        {
          W[q][c]       = 0;
          W_old[q][c]   = 0;
          for (unsigned int d=0; d<dim; ++d)
            {
              grad_W[q][c][d] = 0;
              grad_W_old[q][c][d] = 0;
            }
        }

    for (unsigned int q=0; q<n_q_points; ++q)
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
          const unsigned int c = fe_v.get_fe().system_to_component_index(i).first;

          W[q][c] += independent_local_dof_values[i] *
                     fe_v.shape_value_component(i, q, c);
          W_old[q][c] += old_solution(dof_indices[i]) *
                         fe_v.shape_value_component(i, q, c);

          for (unsigned int d = 0; d < dim; d++)
            {
              grad_W[q][c][d] += independent_local_dof_values[i] *
                                 fe_v.shape_grad_component(i, q, c)[d];
              grad_W_old[q][c][d] += old_solution(dof_indices[i]) *
                                     fe_v.shape_grad_component(i, q, c)[d];
            }
        }

    // allocate and declare
    std::vector < std_cxx11::array <std_cxx11::array <Sacado::Fad::DFad<double>, dim>, PTsolver<dim>::n_components > > flux(n_q_points); // allocate flux
    std::vector <std_cxx11::array <std_cxx11::array <double, dim>, PTsolver<dim>::n_components > > flux_old(n_q_points); // allocate flux_old

    std::vector < std_cxx11::array< Sacado::Fad::DFad<double>, PTsolver<dim>::n_components> > forcing(n_q_points);
    std::vector < std_cxx11::array< double, PTsolver<dim>::n_components> > forcing_old(n_q_points);

    std::vector < std_cxx11::array< Sacado::Fad::DFad<double>, PTsolver<dim>::n_components> > time(n_q_points);
    std::vector < std_cxx11::array< double, PTsolver<dim>::n_components> > time_old(n_q_points);

    std::vector < std_cxx11::array< Sacado::Fad::DFad<double>, PTsolver<dim>::n_components> > coeff(n_q_points);
    std::vector < std_cxx11::array< double, PTsolver<dim>::n_components> > coeff_old(n_q_points);

    std::vector < std_cxx11::array< Sacado::Fad::DFad<double>, PTsolver<dim>::n_components> > extraporo1(n_q_points);
    std::vector < std_cxx11::array< double, PTsolver<dim>::n_components> > extraporo1_old(n_q_points);
    std::vector < std_cxx11::array< Sacado::Fad::DFad<double>, PTsolver<dim>::n_components> > extraporo2(n_q_points);
    std::vector < std_cxx11::array< double, PTsolver<dim>::n_components> > extraporo2_old(n_q_points);

    for (unsigned int q=0; q<n_q_points; ++q)
      {
        PTsolver<dim>::compute_flux_matrix (W_old[q], flux_old[q]);
        PTsolver<dim>::compute_forcing_vector (W_old[q], forcing_old[q]);
        PTsolver<dim>::compute_coeff_vector (W_old[q], coeff_old[q]);
        PTsolver<dim>::compute_extraporo1_vector (W_old[q], extraporo1_old[q]);
        PTsolver<dim>::compute_extraporo2_vector (W_old[q], extraporo2_old[q]);
        PTsolver<dim>::compute_time_vector (W_old[q], time_old[q]);

        PTsolver<dim>::compute_flux_matrix (W[q], flux[q]);
        PTsolver<dim>::compute_forcing_vector (W[q], forcing[q]);
        PTsolver<dim>::compute_coeff_vector (W[q], coeff[q]);
        PTsolver<dim>::compute_extraporo1_vector (W[q], extraporo1[q]);
        PTsolver<dim>::compute_extraporo2_vector (W[q], extraporo2[q]);
        PTsolver<dim>::compute_time_vector (W[q], time[q]);

      }




    for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
      {
        Sacado::Fad::DFad<double> R_i = 0;

        const unsigned int
        component_i = fe_v.get_fe().system_to_component_index(i).first;


        for (unsigned int point=0; point<fe_v.n_quadrature_points; ++point)
          {
            if (parameters.is_stationary == false)  // time stuff
            	R_i += 1./ parameters.time_step *
						(time[point][component_i] - time_old[point][component_id]) *
						fe_v.shape_value_component(i, point, component_i) *
						fe_v.JxW(point);

            for (unsigned int d=0; d<dim; d++) // Flux components
            	R_i += ( parameters.theta * flux[point][component_i][d] +
                       (1.0-parameters.theta) * flux_old[point][component_i][d] ) *
                     fe_v.shape_grad_component(i, point,
                    		 component_i)[d] *
                     fe_v.JxW(point);

            for (unsigned int d=0; d<dim; d++)
            	R_i +=  ( parameters.theta * coeff[point][component_i] * grad_W[point][component_i][d] +
										 (1.0-parameters.theta) * coeff_old[point][component_i] * grad_W_old[point][component_i][d] ) *
										  fe_v.shape_grad_component(i, point, component_i)[d] *
										  fe_v.JxW(point);

            for (unsigned int d=0; d<dim; d++) // got rid of if statement as extraporo[component_i] is 0 for non-porosity component ones
                R_i += ( parameters.theta * ( extraporo1[point][component_i]*grad_W[point][component_i][d]
																		+ extraporo2[point][component_i]*grad_W[point][PTsolver<dim>::VES_component][d] ) +
                                   (1.0-parameters.theta) * ( extraporo1_old[point][component_i]*grad_W_old[point][component_i][d]
																		+ extraporo2_old[point][component_i]*grad_W_old[point][PTsolver<dim>::VES_component][d] ) ) *
                                   		fe_v.shape_value_component(i, point, component_i) *
										fe_v.JxW(point);

          }


        for (unsigned int k=0; k<dofs_per_cell; ++k)
          residual_derivatives[k] = R_i.fastAccessDx(k);
        system_matrix.add(dof_indices[i], dof_indices, residual_derivatives);
        right_hand_side(dof_indices[i]) -= R_i.val();


  }

}








int main (int argc, char *argv[])
{

      using namespace dealii;
      using namespace Step33;



  return 0;
}

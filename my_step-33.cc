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

  template <int dim>
  struct specifics
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
                              specifics<dim>::n_components > &flux)
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
      std_cxx11::array <std_cxx11::array <typename InputVector::value_type, dim>, specifics<dim>::n_components > iflux, oflux;

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

		 void move_mesh();
		 void create_grid();
  };


  template <int dim>
  void PTsolver<dim>::create_grid()
  {
	  GridGenerator::hyper_cube (triangulation, 0, 1);
	  // see step-7 for how to check if a boundary equals the boundary id

	    for (typename Triangulation<dim>::active_cell_iterator
	         cell=triangulation.begin_active();
	         cell!=triangulation.end(); ++cell)
	      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
	        if (cell->face(f)->at_boundary())
	          {
	            const Point<dim> face_center = cell->face(f)->center();

	            if (face_center[dim] == 0)  // z coordinate of mid point of a face is 0
	              cell->face(f)->set_boundary_id (0);

	            else if (face_center[dim] == 1)
	              cell->face(f)->set_boundary_id (1);

	            else
	              cell->face(f)->set_boundary_id (2);
	          }

	    triangulation.refine_global (3);
  }


  template <int dim>
  PTsolver<dim>::PTsolver (const char *input_filename)
    :
    mapping (),
    fe (FE_Q<dim>(1), specifics<dim>::n_components),
    dof_handler (triangulation),
    quadrature (2),
    face_quadrature (2),
    verbose_cout (std::cout, false)
  {
    ParameterHandler prm;
    Parameters::AllParameters<dim>::declare_parameters (prm);

    prm.read_input (input_filename);
    parameters.parse_parameters (prm);

    verbose_cout.set_condition (parameters.output == Parameters::Solver::verbose);
  }



  template <int dim>
  void PTsolver<dim>::setup_system ()
  {
    DynamicSparsityPattern dsp (dof_handler.n_dofs(),
                                dof_handler.n_dofs());
    DoFTools::make_sparsity_pattern (dof_handler, dsp);

    system_matrix.reinit (dsp);
  }


  template <int dim>
  void PTsolver<dim>::assemble_system ()
  {
    const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;

    std::vector<types::global_dof_index> dof_indices (dofs_per_cell);
    std::vector<types::global_dof_index> dof_indices_neighbor (dofs_per_cell);

    const UpdateFlags update_flags               = update_values
                                                   | update_gradients
                                                   | update_q_points
                                                   | update_JxW_values,
                                                   face_update_flags          = update_values
                                                       | update_q_points
                                                       | update_JxW_values
                                                       | update_normal_vectors,
                                                       neighbor_face_update_flags = update_values;

    FEValues<dim>        fe_v                  (mapping, fe, quadrature,
                                                update_flags);
    FEFaceValues<dim>    fe_v_face             (mapping, fe, face_quadrature,
                                                face_update_flags);
    FESubfaceValues<dim> fe_v_subface          (mapping, fe, face_quadrature,
                                                face_update_flags);
    FEFaceValues<dim>    fe_v_face_neighbor    (mapping, fe, face_quadrature,
                                                neighbor_face_update_flags);
    FESubfaceValues<dim> fe_v_subface_neighbor (mapping, fe, face_quadrature,
                                                neighbor_face_update_flags);

    typename DoFHandler<dim>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    for (; cell!=endc; ++cell)
      {
        fe_v.reinit (cell);
        cell->get_dof_indices (dof_indices);

        assemble_cell_term(fe_v, dof_indices);

        for (unsigned int face_no=0; face_no<GeometryInfo<dim>::faces_per_cell;
             ++face_no)
          if (cell->at_boundary(face_no))
            {
              fe_v_face.reinit (cell, face_no);
              assemble_face_term (face_no, fe_v_face,
                                  fe_v_face,
                                  dof_indices,
                                  std::vector<types::global_dof_index>(),
                                  true,
                                  cell->face(face_no)->boundary_id(),
                                  cell->face(face_no)->diameter());
            }

          else
            {
              if (cell->neighbor(face_no)->has_children())
                {
                  const unsigned int neighbor2=
                    cell->neighbor_of_neighbor(face_no);

                  for (unsigned int subface_no=0;
                       subface_no < cell->face(face_no)->n_children();
                       ++subface_no)
                    {
                      const typename DoFHandler<dim>::active_cell_iterator
                      neighbor_child
                        = cell->neighbor_child_on_subface (face_no, subface_no);

                      Assert (neighbor_child->face(neighbor2) ==
                              cell->face(face_no)->child(subface_no),
                              ExcInternalError());
                      Assert (neighbor_child->has_children() == false,
                              ExcInternalError());

                      fe_v_subface.reinit (cell, face_no, subface_no);
                      fe_v_face_neighbor.reinit (neighbor_child, neighbor2);

                      neighbor_child->get_dof_indices (dof_indices_neighbor);

                      assemble_face_term (face_no, fe_v_subface,
                                          fe_v_face_neighbor,
                                          dof_indices,
                                          dof_indices_neighbor,
                                          false,
                                          numbers::invalid_unsigned_int,
                                          neighbor_child->face(neighbor2)->diameter());
                    }
                }

              else if (cell->neighbor(face_no)->level() != cell->level())
                {
                  const typename DoFHandler<dim>::cell_iterator
                  neighbor = cell->neighbor(face_no);
                  Assert(neighbor->level() == cell->level()-1,
                         ExcInternalError());

                  neighbor->get_dof_indices (dof_indices_neighbor);

                  const std::pair<unsigned int, unsigned int>
                  faceno_subfaceno = cell->neighbor_of_coarser_neighbor(face_no);
                  const unsigned int neighbor_face_no    = faceno_subfaceno.first,
                                     neighbor_subface_no = faceno_subfaceno.second;

                  Assert (neighbor->neighbor_child_on_subface (neighbor_face_no,
                                                               neighbor_subface_no)
                          == cell,
                          ExcInternalError());

                  fe_v_face.reinit (cell, face_no);
                  fe_v_subface_neighbor.reinit (neighbor,
                                                neighbor_face_no,
                                                neighbor_subface_no);

                  assemble_face_term (face_no, fe_v_face,
                                      fe_v_subface_neighbor,
                                      dof_indices,
                                      dof_indices_neighbor,
                                      false,
                                      numbers::invalid_unsigned_int,
                                      cell->face(face_no)->diameter());
                }
            }
      }
  }



// PTsolver::PTsolver here

  template <int dim>
  void
  PTsolver<dim>::
  assemble_cell_term (const FEValues<dim>             &fe_v,
                      const std::vector<types::global_dof_index> &dof_indices)
  {
    const unsigned int dofs_per_cell = fe_v.dofs_per_cell;
    const unsigned int n_q_points    = fe_v.n_quadrature_points;

    Table<2,Sacado::Fad::DFad<double> >
    W (n_q_points, specifics<dim>::n_components);

    Table<2,double>
    W_old (n_q_points, specifics<dim>::n_components);

    Table<3,Sacado::Fad::DFad<double> >
    grad_W (n_q_points, specifics<dim>::n_components, dim);

    Table<3,double>
    grad_W_old(n_q_points, specifics<dim>::n_components, dim);

    std::vector<double> residual_derivatives (dofs_per_cell);

    std::vector<Sacado::Fad::DFad<double> > independent_local_dof_values(dofs_per_cell);
    for (unsigned int i=0; i<dofs_per_cell; ++i)
      independent_local_dof_values[i] = current_solution(dof_indices[i]);

    for (unsigned int i=0; i<dofs_per_cell; ++i)
    	independent_local_dof_values[i].diff (i, dofs_per_cell);

    ///////////////////////////////////////////////////////////////////////////
    // initialise the values of W etc then actually calculate

    for (unsigned int q=0; q<n_q_points; ++q)
      for (unsigned int c=0; c<specifics<dim>::n_components; ++c)
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
    std::vector < std_cxx11::array <std_cxx11::array <Sacado::Fad::DFad<double>, dim>, specifics<dim>::n_components > > flux(n_q_points); // allocate flux
    std::vector <std_cxx11::array <std_cxx11::array <double, dim>, specifics<dim>::n_components > > flux_old(n_q_points); // allocate flux_old

    std::vector < std_cxx11::array< Sacado::Fad::DFad<double>, specifics<dim>::n_components> > forcing(n_q_points);
    std::vector < std_cxx11::array< double, specifics<dim>::n_components> > forcing_old(n_q_points);

    std::vector < std_cxx11::array< Sacado::Fad::DFad<double>, specifics<dim>::n_components> > time(n_q_points);
    std::vector < std_cxx11::array< double, specifics<dim>::n_components> > time_old(n_q_points);

    std::vector < std_cxx11::array< Sacado::Fad::DFad<double>, specifics<dim>::n_components> > coeff(n_q_points);
    std::vector < std_cxx11::array< double, specifics<dim>::n_components> > coeff_old(n_q_points);

    std::vector < std_cxx11::array< Sacado::Fad::DFad<double>, specifics<dim>::n_components> > extraporo1(n_q_points);
    std::vector < std_cxx11::array< double, specifics<dim>::n_components> > extraporo1_old(n_q_points);
    std::vector < std_cxx11::array< Sacado::Fad::DFad<double>, specifics<dim>::n_components> > extraporo2(n_q_points);
    std::vector < std_cxx11::array< double, specifics<dim>::n_components> > extraporo2_old(n_q_points);

    for (unsigned int q=0; q<n_q_points; ++q)
      {
        specifics<dim>::compute_flux_matrix (W_old[q], flux_old[q]);
        specifics<dim>::compute_forcing_vector (W_old[q], forcing_old[q]);
        specifics<dim>::compute_coeff_vector (W_old[q], coeff_old[q]);
        specifics<dim>::compute_extraporo1_vector (W_old[q], extraporo1_old[q]);
        specifics<dim>::compute_extraporo2_vector (W_old[q], extraporo2_old[q]);
        specifics<dim>::compute_time_vector (W_old[q], time_old[q]);

        specifics<dim>::compute_flux_matrix (W[q], flux[q]);
        specifics<dim>::compute_forcing_vector (W[q], forcing[q]);
        specifics<dim>::compute_coeff_vector (W[q], coeff[q]);
        specifics<dim>::compute_extraporo1_vector (W[q], extraporo1[q]);
        specifics<dim>::compute_extraporo2_vector (W[q], extraporo2[q]);
        specifics<dim>::compute_time_vector (W[q], time[q]);

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
																		+ extraporo2[point][component_i]*grad_W[point][specifics<dim>::VES_component][d] ) +
                                   (1.0-parameters.theta) * ( extraporo1_old[point][component_i]*grad_W_old[point][component_i][d]
																		+ extraporo2_old[point][component_i]*grad_W_old[point][specifics<dim>::VES_component][d] ) ) *
                                   		fe_v.shape_value_component(i, point, component_i) *
										fe_v.JxW(point);

          }


        for (unsigned int k=0; k<dofs_per_cell; ++k)
          residual_derivatives[k] = R_i.fastAccessDx(k);
        system_matrix.add(dof_indices[i], dof_indices, residual_derivatives);
        right_hand_side(dof_indices[i]) -= R_i.val();


  }

}

  template <int dim>
  void
  PTsolver<dim>::assemble_face_term(const unsigned int           face_no,
                                           const FEFaceValuesBase<dim> &fe_v,
                                           const FEFaceValuesBase<dim> &fe_v_neighbor,
                                           const std::vector<types::global_dof_index>   &dof_indices,
                                           const std::vector<types::global_dof_index>   &dof_indices_neighbor,
                                           const bool                   external_face,
                                           const unsigned int           boundary_id,
                                           const double                 face_diameter)
  {
    const unsigned int n_q_points = fe_v.n_quadrature_points;
    const unsigned int dofs_per_cell = fe_v.dofs_per_cell;

    std::vector<Sacado::Fad::DFad<double> >
    independent_local_dof_values (dofs_per_cell),
                                 independent_neighbor_dof_values (external_face == false ?
                                     dofs_per_cell :
                                     0);

    const unsigned int n_independent_variables = (external_face == false ?
                                                  2 * dofs_per_cell :
                                                  dofs_per_cell);

    for (unsigned int i = 0; i < dofs_per_cell; i++)
      {
        independent_local_dof_values[i] = current_solution(dof_indices[i]);
        independent_local_dof_values[i].diff(i, n_independent_variables);
      }

    if (external_face == false)
      for (unsigned int i = 0; i < dofs_per_cell; i++)
        {
          independent_neighbor_dof_values[i]
            = current_solution(dof_indices_neighbor[i]);
          independent_neighbor_dof_values[i]
          .diff(i+dofs_per_cell, n_independent_variables);
        }


    Table<2,Sacado::Fad::DFad<double> >
    Wplus (n_q_points, specifics<dim>::n_components),
          Wminus (n_q_points, specifics<dim>::n_components);
    Table<2,double>
    Wplus_old(n_q_points, specifics<dim>::n_components),
              Wminus_old(n_q_points, specifics<dim>::n_components);

    for (unsigned int q=0; q<n_q_points; ++q)
      for (unsigned int i=0; i<dofs_per_cell; ++i)
        {
          const unsigned int component_i = fe_v.get_fe().system_to_component_index(i).first;
          Wplus[q][component_i] +=  independent_local_dof_values[i] *
                                    fe_v.shape_value_component(i, q, component_i);
          Wplus_old[q][component_i] +=  old_solution(dof_indices[i]) *
                                        fe_v.shape_value_component(i, q, component_i);
        }

    if (external_face == false)
      {
        for (unsigned int q=0; q<n_q_points; ++q)
          for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
              const unsigned int component_i = fe_v_neighbor.get_fe().
                                               system_to_component_index(i).first;
              Wminus[q][component_i] += independent_neighbor_dof_values[i] *
                                        fe_v_neighbor.shape_value_component(i, q, component_i);
              Wminus_old[q][component_i] += old_solution(dof_indices_neighbor[i])*
                                            fe_v_neighbor.shape_value_component(i, q, component_i);
            }
      }
    else
      {
        Assert (boundary_id < Parameters::AllParameters<dim>::max_n_boundaries,
                ExcIndexRange (boundary_id, 0,
                               Parameters::AllParameters<dim>::max_n_boundaries));

        std::vector<Vector<double> >
        boundary_values(n_q_points, Vector<double>(specifics<dim>::n_components));
        parameters.boundary_conditions[boundary_id]
        .values.vector_value_list(fe_v.get_quadrature_points(),
                                  boundary_values);

        for (unsigned int q = 0; q < n_q_points; q++)
          {
            specifics<dim>::compute_Wminus (parameters.boundary_conditions[boundary_id].kind,
                                                 fe_v.normal_vector(q),
                                                 Wplus[q],
                                                 boundary_values[q],
                                                 Wminus[q]);
            specifics<dim>::compute_Wminus (parameters.boundary_conditions[boundary_id].kind,
                                                 fe_v.normal_vector(q),
                                                 Wplus_old[q],
                                                 boundary_values[q],
                                                 Wminus_old[q]);
          }
      }



    std::vector< std_cxx11::array < Sacado::Fad::DFad<double>, specifics<dim>::n_components> >  normal_fluxes(n_q_points);
    std::vector< std_cxx11::array < double, specifics<dim>::n_components> >  normal_fluxes_old(n_q_points);

    double alpha;

    switch (parameters.stabilization_kind)
      {
      case Parameters::Flux::constant:
        alpha = parameters.stabilization_value;
        break;
      case Parameters::Flux::mesh_dependent:
        alpha = face_diameter/(2.0*parameters.time_step);
        break;
      default:
        Assert (false, ExcNotImplemented());
        alpha = 1;
      }

    for (unsigned int q=0; q<n_q_points; ++q)
      {
        specifics<dim>::numerical_normal_flux(fe_v.normal_vector(q),
                                                   Wplus[q], Wminus[q], alpha,
                                                   normal_fluxes[q]);
        specifics<dim>::numerical_normal_flux(fe_v.normal_vector(q),
                                                   Wplus_old[q], Wminus_old[q], alpha,
                                                   normal_fluxes_old[q]);
      }

    std::vector<double> residual_derivatives (dofs_per_cell);
    for (unsigned int i=0; i<fe_v.dofs_per_cell; ++i)
      if (fe_v.get_fe().has_support_on_face(i, face_no) == true)
        {
          Sacado::Fad::DFad<double> R_i = 0;

          for (unsigned int point=0; point<n_q_points; ++point)
            {
              const unsigned int
              component_i = fe_v.get_fe().system_to_component_index(i).first;

              R_i += ( parameters.theta * normal_fluxes[point][component_i] +
                       (1.0 - parameters.theta) * normal_fluxes_old[point][component_i] ) *
                     fe_v.shape_value_component(i, point, component_i) *
                     fe_v.JxW(point);
            }

          for (unsigned int k=0; k<dofs_per_cell; ++k)
            residual_derivatives[k] = R_i.fastAccessDx(k);
          system_matrix.add(dof_indices[i], dof_indices, residual_derivatives);

          if (external_face == false)
            {
              for (unsigned int k=0; k<dofs_per_cell; ++k)
                residual_derivatives[k] = R_i.fastAccessDx(dofs_per_cell+k);
              system_matrix.add (dof_indices[i], dof_indices_neighbor,
                                 residual_derivatives);
            }

          right_hand_side(dof_indices[i]) -= R_i.val();
        }
  }



  template <int dim>
  void PTsolver<dim>::run ()
  {
     // create grid. not reading in in this case.

    create_grid();

    dof_handler.clear();
    dof_handler.distribute_dofs (fe);

    old_solution.reinit (dof_handler.n_dofs());
    current_solution.reinit (dof_handler.n_dofs());
    predictor.reinit (dof_handler.n_dofs());
    right_hand_side.reinit (dof_handler.n_dofs());

    setup_system();

    VectorTools::interpolate(dof_handler,
                             parameters.initial_conditions, old_solution);
    current_solution = old_solution;
    predictor = old_solution;

    if (parameters.do_refine == true)
      for (unsigned int i=0; i<parameters.shock_levels; ++i)
        {
          Vector<double> refinement_indicators (triangulation.n_active_cells());

          compute_refinement_indicators(refinement_indicators);
          refine_grid(refinement_indicators);

          setup_system();

          VectorTools::interpolate(dof_handler,
                                   parameters.initial_conditions, old_solution);
          current_solution = old_solution;
          predictor = old_solution;
        }

    output_results ();

    Vector<double> newton_update (dof_handler.n_dofs());

    double time = 0;
    double next_output = time + parameters.output_step;

    predictor = old_solution;
    while (time < parameters.final_time)
      {
        std::cout << "T=" << time << std::endl
                  << "   Number of active cells:       "
                  << triangulation.n_active_cells()
                  << std::endl
                  << "   Number of degrees of freedom: "
                  << dof_handler.n_dofs()
                  << std::endl
                  << std::endl;

        std::cout << "   NonLin Res     Lin Iter       Lin Res" << std::endl
                  << "   _____________________________________" << std::endl;

        unsigned int nonlin_iter = 0;
        current_solution = predictor;
        while (true)
          {
            system_matrix = 0;

            right_hand_side = 0;
            assemble_system ();

            const double res_norm = right_hand_side.l2_norm();
            if (std::fabs(res_norm) < 1e-10)
              {
                std::printf("   %-16.3e (converged)\n\n", res_norm);
                break;
              }
            else
              {
                newton_update = 0;

                std::pair<unsigned int, double> convergence
                  = solve (newton_update);

                current_solution += newton_update;

                std::printf("   %-16.3e %04d        %-5.2e\n",
                            res_norm, convergence.first, convergence.second);
              }

            ++nonlin_iter;
            AssertThrow (nonlin_iter <= 10,
                         ExcMessage ("No convergence in nonlinear solver"));
          }

        time += parameters.time_step;

        if (parameters.output_step < 0)
          output_results ();
        else if (time >= next_output)
          {
            output_results ();
            next_output += parameters.output_step;
          }

        predictor = current_solution;
        predictor.sadd (2.0, -1.0, old_solution);

        old_solution = current_solution;

        if (parameters.do_refine == true)
          {
            Vector<double> refinement_indicators (triangulation.n_active_cells());
            compute_refinement_indicators(refinement_indicators);

            refine_grid(refinement_indicators);
            setup_system();

            newton_update.reinit (dof_handler.n_dofs());
          }
      }
  }





int main (int argc, char *argv[])
{

      using namespace dealii;
      using namespace Step33;



  return 0;
}

/* ---------------------------------------------------------------------
 *
 * Copyright (C) 1999 - 2015 by the deal.II authors
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

 */


#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_in.h>

#include <iostream>
#include <fstream>
#include <cmath>

using namespace dealii;

template <int dim>
void print_mesh_info(const Triangulation<dim> &tria,
                     const std::string        &filename)
{
  std::cout << "Mesh info:" << std::endl
            << " dimension: " << dim << std::endl
            << " no. of cells: " << tria.n_active_cells() << std::endl;
  {
    std::map<unsigned int, unsigned int> boundary_count;
    typename Triangulation<dim>::active_cell_iterator
    cell = tria.begin_active(),
    endc = tria.end();
    for (; cell!=endc; ++cell)
      {
        for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
          {
            if (cell->face(face)->at_boundary())
              boundary_count[cell->face(face)->boundary_id()]++;
          }
      }
    std::cout << " boundary indicators: ";
    for (std::map<unsigned int, unsigned int>::iterator it=boundary_count.begin();
         it!=boundary_count.end();
         ++it)
      {
        std::cout << it->first << "(" << it->second << " times) ";
      }
    std::cout << std::endl;
  }
  std::ofstream out (filename.c_str());
  GridOut grid_out;
  grid_out.write_eps (tria, out);
  std::cout << " written to " << filename
            << std::endl
            << std::endl;
}

template <int dim>
class gridtest
{
public:

	void run ();

private:
	void grid ();
	void grid_2 ();
	void trap();

};

template <int dim>
void gridtest<dim>::grid ()
{
  Triangulation<dim> triangulation1;
  Triangulation<dim> triangulation2;

  const Point<dim> p1;
  const Point<dim> p2;

  GridGenerator::hyper_cube (triangulation1, 0,1);
  triangulation1.refine_global (4);

  std::ofstream out1 ("grid-1.eps");
  GridOut grid_out1;
  grid_out1.write_eps (triangulation1, out1);
  std::cout << "Grid written to grid-1.eps" << std::endl;


  GridGenerator::truncated_cone	(triangulation2, 1.0, 0.5, 1.0 );
  triangulation2.refine_global (4);

  std::ofstream out2 ("cone.eps");
  GridOut grid_out2;
  grid_out1.write_eps (triangulation2, out2);
  std::cout << "Grid written to cone.eps" << std::endl;
}

template <int dim>
void gridtest<dim>::grid_2 ()
{
  Triangulation<dim> hole1;
  GridGenerator::hyper_cube_with_cylindrical_hole (hole1, 0.25, 1.0);

  Triangulation<dim> hole2;
  std::vector< unsigned int > repetitions(2);
  repetitions[0]=3;
  repetitions[1]=2;
  GridGenerator::subdivided_hyper_rectangle (hole2, repetitions,
                                             Point<dim>(1.0,-1.0),
                                             Point<dim>(4.0,1.0));
  Triangulation<dim> hole;
  GridGenerator::merge_triangulations (hole1, hole2, hole);
  print_mesh_info(hole, "grid-hole.eps");
}

template<int dim>
void gridtest<dim>::trap ()
{
	Triangulation<dim> tri;
	Triangulation<dim> rec;

	Point<dim> v1 (0.0, -0.5);
	Point<dim> v2 (0.0, 0.5);
	Point<dim> v3 (-1.0, 0.5);
	std::vector< Point<dim> > v = {v1,v2,v3};

	GridGenerator::simplex	(tri, v);
	std::vector< unsigned int > repetitions(2);
	repetitions[0]=2;
	repetitions[1]=2;
    GridGenerator::subdivided_hyper_rectangle (rec, repetitions,
	                                             Point<dim>(0.0, -0.5),
	                                             Point<dim>(1.0, 0.5));


	Triangulation<dim> trap;
	GridGenerator::merge_triangulations (tri, rec, trap);

	typename Triangulation<dim>::cell_iterator
	          cell = trap.begin (),
	          endc = trap.end();
	          for (; cell!=endc; ++cell)
	            for (unsigned int face_number=0;
	                 face_number<GeometryInfo<dim>::faces_per_cell;
	                 ++face_number)
	              if ((std::fabs(cell->face(face_number)->center()(0) - (1.0)) < 1e-12)) // if x=1.0, RHS boundary
	                cell->face(face_number)->set_boundary_id (2);
	              else if ((std::fabs(cell->face(face_number)->center()(1) - (-0.5)) < 1e-12)) // if y=-0.5, bottom boundary
	            	  cell->face(face_number)->set_boundary_id (3);
	              else if ((std::fabs(cell->face(face_number)->center()(1) - (0.5)) < 1e-12)) // if y=0.5, top boundary
	            	  cell->face(face_number)->set_boundary_id (1);


	trap.refine_global (2);
	print_mesh_info(trap, "trap.eps");


}




template <int dim>
void gridtest<dim>::run ()
{
	//grid ();
	//grid_2 ();
	trap ();


}



int main ()
{


	gridtest<2> test;
	test.run ();

	return 0;

}

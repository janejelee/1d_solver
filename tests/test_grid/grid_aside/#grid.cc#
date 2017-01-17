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

#include <iostream>
#include <fstream>
#include <cmath>

using namespace dealii;

template <int dim>
class gridtest
{
public:

	void run ();

private:
	void grid ();

	Triangulation<dim> triangulation;
};

template <int dim>
void gridtest<dim>::grid ()
{
  Triangulation<dim> triangulation;
  const Point<dim> p1;
  const Point<dim> p2;

  GridGenerator::hyper_cube (triangulation, 0,1);
  triangulation.refine_global (4);

  std::ofstream out ("grid-1.eps");
  GridOut grid_out;
  grid_out.write_eps (triangulation, out);
  std::cout << "Grid written to grid-1.eps" << std::endl;
}

template <int dim>
void gridtest<dim>::run ()
{
	grid ();
}




int main ()
{


	gridtest<1> test;
	test.run ();

	return 0;

}

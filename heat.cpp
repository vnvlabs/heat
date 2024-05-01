#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include "VnV.h"

using namespace std;

int main(int argc, char *argv[]);


/**
 *
 *    @title Solving the Steady State Heat Equation using Finite Differences. 
 * 
 *    This code solves the steady state heat equation on a rectangular region.
 *
 *    The sequential version of this program needs approximately
 *    18/epsilon iterations to complete.
 *
 *    The steady state solution to the discrete heat equation satisfies the
 *    following condition at an interior grid point:
 *
 *      W[Central] = (1/4) * ( W[North] + W[South] + W[East] + W[West] )
 *
 *    where "Central" is the index of the grid point, "North" is the index
 *    of its immediate neighbor to the "north", and so on.
 *
 *    Given an approximate solution of the steady state heat equation, a
 *    "better" solution is given by replacing each interior point by the
 *    average of its 4 neighbors - in other words, by using the condition
 *    as an ASSIGNMENT statement:
 *
 *      W[Central]  <=  (1/4) * ( W[North] + W[South] + W[East] + W[West] )
 *
 *    If this process is repeated often enough, the difference between
 *    successive estimates of the solution will go to zero.
 *
 *    This program carries out such an iteration, using a tolerance specified by
 *    the user, and writes the final estimate of the solution to a file that can
 *    be used for graphic processing.  
 *
 */
INJECTION_EXECUTABLE(Heat,"{}");

class OptionsStruct {
public:  
  double eps, l,r,t,b;
  int N,M, max_it;
};


/**
 * @title Heated Plate Code with VnV 
 * 
 * Version: 1.0.0
 * Author: Ben O'Neill
 * Licence: M.I.T
 * 
 */
INJECTION_OPTIONS(Heat, R"({
  "type" : "object",
  "properties" : { 
     "epsilon" : {"type" : "number" , "default" : 0.001 , "description" : "The simulation will run until the difference between two iteration is less than this value" },
     "max_it" : {"type" : "integer" , "default" : 100000 , "description" : "The Maximum number of iterations to take"},
     "width" : {"type" : "integer" , "default" : 500 , "description" : "The height of the rectangular domain", "min" : 10 , "max" : 1000 },
     "height" : {"type" : "integer" , "default" : 100 , "description" : "The width of the rectangular domain", "min" : 10 , "max" : 1000 },
     "left" : {"type" : "number" , "default" : 100.0 , "description" : "The Heat on the left boundary" },
     "right" : {"type" : "number" , "default" : 100.0 , "description" : "The Heat on the left boundary" },
     "top" : {"type" : "number" , "default" : 0.0 , "description" : "The Heat on the left boundary" },
     "bottom" : {"type" : "number" , "default" : 100.0 , "description" : "The Heat on the left boundary"}
  }
})",OptionsStruct) {

  OptionsStruct *c  = new OptionsStruct();
  c->l = config.value("left",100.0);
  c->r = config.value("right",100.0);
  c->t = config.value("top",  0.0);
  c->b = config.value("bottom",100.0);
  c->eps = config.value("epsilon",0.1);
  c->M = config.value("width", 500);
  c->N = config.value("height",100);
  c->max_it = config.value("max_it", 1000000);
  return c;

}

class Plate {
  int M,N;
  std::vector<double> vals;
public:

  enum class BoundaryType { TOP, LEFT, BOTTOM, RIGHT };

  Plate(int M_, int N_) : M(M_), N(N_), vals(M_*N_) {}

  int width() {return M;}
  int height() {return N;}

  std::vector<double>& data() { return vals; }

  double& get(int i, int j ) {
    return get(i,j,vals);
  }

  double& get(int i, int j, std::vector<double>& vec ) {
    if ( i < 0 || i > N) throw INJECTION_EXCEPTION("Row Index %d is out of bounds [0,%d]", i, N);
    if ( j < 0 || j > M) throw INJECTION_EXCEPTION("Column Index %d is out of bounds [0,%d]", j, M);
    return vec[i*M + j];
  }

  void set(int i, int j, const double& val) {
    set(i,j,val,vals);
  }

  void set(int i, int j, const double& val, std::vector<double>& vec) {
    if ( i < 0 || i > N) throw INJECTION_EXCEPTION("Row Index %d is out of bounds [0,%d]", i, N);
    if ( j < 0 || j > M) throw INJECTION_EXCEPTION("Column Index %d is out of bounds [0,%d]", j, M);
    vec[i*M + j] = val;
  }

  void setBoundary(BoundaryType type, double val, bool include_edges) {
     
      int start = include_edges ? 0 : 1;

      if (type == BoundaryType::TOP || type == BoundaryType::BOTTOM) {
        int i =  type == BoundaryType::TOP ? 0 : (N-1);
        for ( int j = start; j < M - start ; j++) { set(i,j,val); }    
      } else {
        int j = type == BoundaryType::LEFT ? 0 : M-1 ;
        for ( int i = start; i < N - start ; i++) { set( i, j, val); }    
      }
  }
  
  void setBoundary(double left, double right, double top, double bottom, bool prioritize_top) {
     setBoundary(BoundaryType::LEFT, left, !prioritize_top);
     setBoundary(BoundaryType::RIGHT, right, !prioritize_top);
     setBoundary(BoundaryType::TOP, top, prioritize_top);
     setBoundary(BoundaryType::BOTTOM, bottom, prioritize_top);
  }

  void randomizeInterior(int lower = 0, int upper = 100) {
    if (lower >= upper) throw INJECTION_EXCEPTION("Lower (%d) must be less than upper (%d)", lower, upper);

    for (int i = 1; i < N - 1; i++ ) {
      for (int j = 1; j < M -1; j++ ) {
          set(i,j,(((float) rand()) / (float) RAND_MAX)*(upper-lower) + lower);
      }
    }

  }

  void dump() {

    for (int i = 0; i < N; i++ ) {
      for (int j = 0; j < M; j++ ) {
          std::cout << get(i,j) << " " ;
      }
      std::cout << std::endl;
    }
  }

  double iterate() {

    double diff = 0; 
    
    std::vector<double> tmp = vals; // copy. (Could remove this by toggling between two vectors)
    for (int i = 1; i < N - 1; i++) {
      for (int j = 1; j < M - 1; j++) {
         set(i,j, ( get(i-1,j,tmp) + get(i+1,j,tmp) + get(i,j - 1,tmp) + get(i,j+1,tmp) ) / 4.0 ) ;

         double change = fabs( get(i,j) - get(i,j,tmp) );
         if (diff < change) {
           diff = change;
         }
      }
    }
    return diff;
  }

};

/**
 * @title Animating the Solution
 * @shortTitle Anitmation
 * 
 * The figure below shows an anitmation of the solution process. Images are 
 * drawn on iterations by power of 2. . 
 * 
 * .. vnv-animation::
 *      :trace.sol: contour
 *      :layout.title.text: Solution 
 *      :layout.yaxis.autorange: reversed
 *      :labels: {{iterations}}
 *      :prefix: Iteration
 *      :values: {{solution}}
 *      :sol.z: ${i}
 *      
 * .. vnv-animation::
 *      :trace.sol: surface
 *      :layout.title.text: Solution 
 *      :labels: {{iterations}}
 *      :prefix: Iteration
 *      :values: {{solution}}
 *      :sol.z: ${i}
 * 
 */
INJECTION_TEST(Heat, Animate) {

  int freq = m_config.getAdditionalParameters().value("freq",-1);
  int iter = GetRef("iterations", int);
  double err = GetRef("diff", double);

  if (type == InjectionPointType::Begin && m_config.getAdditionalParameters().value("surface",false)) {
   engine->Put("surface", "surface");
  }

  if (type == InjectionPointType::End || iter == 0 || !(iter & (iter-1)) ) {

    auto p = GetRef("plate", Plate);
  
    auto g = std::make_pair(p.height(),p.width());
    auto o = std::make_pair(0,0);


   engine->Put_Matrix("solution", p.height(), p.data(), g,o);
   engine->Put("iterations", iter);
   engine->Put("error", err);
  }
  return SUCCESS;
}


int main(int argc, char *argv[]) {

  double diff;
  int i,j;
  int iterations;
  int iterations_print;
  int success;
  
 /**
 *   @title Steady State Heat Distribution for a heated plate. 
 * 
 *   This code solves the steady state heat equation for a retangular plate
 *   with heated sides. The code assumes the sequential version of this program needs approximately
 *   18/epsilon iterations to complete.
 *
 * 
 */
 INJECTION_INITIALIZE(Heat, &argc, &argv, (argc > 1) ? argv[1] : "vv-input.json");

  OptionsStruct* opts = (OptionsStruct*) INJECTION_GET_CONFIG(Heat);
  
  int M = opts->M;
  int N = opts->N;
  
  // Create the plate;
  Plate plate(opts->M,opts->N);
  
  //Set the bounary
  plate.setBoundary(opts->l, opts->r, opts->t, opts->b, true);
  
  // Randomize the interior. 
  plate.randomizeInterior();

  /**
   * @title Defining the Problem
   * @shortTitle Problem Configuration
   * 
   * The goal of this code was to solve the steady state heat equation for 
   * a rectangular plate being heated uniformly along each edge. The boundary 
   * conditions along each edge were:
   * 
   *   -  Left: :vnv:`L` units 
   *   -  Right: :vnv:`R` units
   *   -  Top: :vnv:`T` units
   *   -  Bottom: :vnv:`B` units
   * 
   * The region was discretized into a [ :vnv:`M` , :vnv:`N` ] grid with 
   * random values being used for initialize all interior points. After applying the
   * boundary conditions and randomizing the interior points, the inital guess was
   * as shown below
   *     
   * .. vnv-plotly::
   *     :trace.sol: contour
   *     :sol.z: {{as_json(init)}}
   *     :layout.title.text: Basic contour plot
   *     :layout.yaxis.autorange: reversed
   *  
   */
  INJECTION_POINT(Heat, VWORLD,Setup, VNV_CALLBACK{
     data.engine->Put("L", opts->l);
     data.engine->Put("R", opts->r);
     data.engine->Put("T", opts->t);
     data.engine->Put("B", opts->b);
     data.engine->Put("M", opts->M);
     data.engine->Put("N", opts->N);
      auto g = std::make_pair(opts->N,opts->M);
      auto o = std::make_pair(0,0);
      
     data.engine->Put_Matrix("init", opts->N, plate.data(), g,o);
      
  }, opts);

  

 /**
 *   @title Iterative Solve:
 *   @shortTitle Solve 
 *  
 *   The steady state solution to the discrete heat equation satisfies the
 *   following condition at an interior grid point:
 *
 *   .. vnv-math::
 *      
 *        W[i,j] = (1/4) * ( W[i+1,j] + W[i-1,j] + W[i,j+1] + W[i,j-1] )
 *
 *   The solution is found by iterativly applying this condition at every 
 *   iterior point. This code uses an assignment operation to do this. All
 *   assignments are completed using the solution at the previous iteration. That is,
 *   a copy of the solution is made at each step and the final solution is returned.
 *
 *  
 *   .. note:: 
 * 
 *       Convergence will also stop if the iteration count reaches the maximum number
 *       of iterations selected by the user (in this case max_it was :vnv:`max_it`)
 *
 *   Convergence:
 *   ------------
 * 
 *   The figure below plots the maximum difference against the number of overall iterations
 *   in the solve.  
 *   
 *   .. vnv-plotly::
 *      :trace.conv: scatter
 *      :conv.y: {{convergence}}
 *      :layout.xaxis.type: log
 *      :layout.xaxis.autorange: true
 *      :layout.yaxis.type: log
 *      :layout.yaxis.autorange: true
 *
 */
  INJECTION_LOOP_BEGIN(Heat,VWORLD,Solve, VNV_CALLBACK {
       data.engine->Put("epsilon", opts->eps);
       data.engine->Put("max_it", opts->max_it);
  }, plate, opts, diff, iterations);
  
  
  iterations = 0;
  diff = opts->eps + 100;
  while (opts->eps <= diff && ++iterations < opts->max_it ) {
    diff = plate.iterate();
    
    INJECTION_LOOP_ITER(Heat,Solve,"Iteration",VNV_CALLBACK {
           data.engine->Put("convergence", diff);
    });

    VnV_Info(Heat, "Error at iteration %d: %f", iterations-1, diff);
  }

  INJECTION_LOOP_END(Heat,Solve,VNV_NOCALLBACK);
/**  @title Steady State Solution:
 *   @shortTitle Solution
 *  
 *   The total number of iterations performed was :vnv:`total`. At the end, the tolerance 
 *   was :vnv:`diff`. 
 * 
 *   The steady state solution is
 * 
 *   .. vnv-plotly::
 *       :trace.sol: contour
 *       :sol.z: {{as_json(solution)}}
 *       :layout.title.text: Steady State solution for Heated Plate.
 *       :layout.yaxis.autorange: reversed 
 * 
 * 
 */  

  INJECTION_POINT(Heat, VWORLD, Solution, VNV_CALLBACK {
    auto g = std::make_pair(opts->N,opts->M);
    auto o = std::make_pair(0,0);
    data.engine->Put_Matrix("solution", opts->N, plate.data(), g,o);
    data.engine->Put("total", iterations-1);  
    data.engine->Put("diff", diff);
  }, plate, opts, diff, iterations)


  INJECTION_FINALIZE(Heat)

  return 0;

}


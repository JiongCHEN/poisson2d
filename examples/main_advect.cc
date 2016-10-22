#include <iostream>
#include <fstream>

#include "src/grid.h"

using namespace std;

// class advection_equation
// {
// public:
//   advection_equation() {
//     // init T, u and v
//   }
//   int solve() {
//     grid2d<double> next_u, next_v, half_u, half_v;
    
//     for (size_t frm = 0; frm < 100; ++frm) {
//       const double curr_t = frm*h_;
      
//       // time invariant velocity field
//       next_u = u_;             
//       next_v = v_;

//       // assumption
//       half_u = (u_+next_u)/2.;
//       half_v = (v_+next_v)/2.;

//       grid2d<double> next_T;
      
//       // solve the internal points
//       for (size_t i = 1; i < T_->nx()-1; ++i) {
//         for (size_t j = 1; j < T_->ny()-1; ++j) {          
//           Vec2 curr_p = T->pos(i, j), prev_p;
//           Vec2 curr_v = Vec2(u_(i, j), v_(i, j));
          
//           int iter = 5;
//           while ( iter-- ) {
//             prev_p = curr_p-0.5*h_*curr_v;
//             curr_vel(0) = half_u.interplinear(prev_p);
//             curr_vel(1) = half_v.interplinear(prev_p);
//           }
          
//           prev_p = curr_p-h_*curr_v;

//           next_T(i, j) = T_->interpcubic(prev_p);
//         }
//       }

//       // assign boundary
//       next_T = ;
      
//       *T_ = next_T;
//     }
//     return 0;
//   }
// private:
//   double h_;
//   std::shared_ptr<grid2d<double>> T_, u_, v_;
// };

int main(int argc, char *argv[])
{
  return 0;
}

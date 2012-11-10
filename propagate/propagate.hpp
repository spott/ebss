#pragma once

namespace propagate
{

//functions that we need:  step (takes the "matrix finding function" 
//                         solve (which just steps through a set number of times),
//

PetscErrorCode
step(PetscReal time, PetscReal dt, std::function<Mat> 

}

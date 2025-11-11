# Simplified Laplacian Mesh Deformation — Implementation Prompt Template

This markdown file defines a **modular prompt** for instructing an LLM to implement a single-linear-solve Laplacian mesh deformation system.  
It follows a clear six-module structure so both human and LLM can collaborate efficiently.

---

## 🧩 Overview

We perform smooth mesh deformation after moving handle vertices, using one linear solve of the Laplacian system (no rotations, no ARAP).  
Language: **Python (numpy + scipy.sparse)**

---

## 🧱 Prompt for LLM

**Copy-paste the following prompt into any code-capable LLM:**

> **Overview:**  
> Implement the *Simplified Laplacian Mesh Deformation* using the six modules below.  
> Follow the module names and input/output structure exactly, one module = one function.  
> Look through the current codebase to understand it structure and the provided API. You should prioritize the provided API than writing a new function.  
> Each module must be small, commented, and independently editable.  
> Use cotangent weights by default. Add small diagonal regularization (1e-8) if needed.  
>
> **Modules:**  
>
> 1. `build_laplacian(p, faces, weight_type="cotangent")`  
> 2. `compute_laplacian_coordinates(L, p)`  
> 3. `apply_constraints(L, δ, handles, u_c)`  
> 4. `solve_positions(L_ff, b_f)`  
> 5. `merge_results(p, free_idx, cons_idx, p_f_solved, u_c)`  
> 6. `run_test()`  
>
> ###### GUI
>
> 1. there should be a "fix" button, when I click the button, it enters "drag points mode" in which I can drag one single vertex and the mesh should deform as I move the points.
> 2. there should be a "move" button. after I click the button, it enters "fix points mode" in which I can select points to fix their positions. These points should be highlighted. 
>
> ###### Terminal output
>
> 1. when I do the select / fix / drag operation, output information in the terminal.
>
> ###### LLM output
>
> 1. ONLY put the terminal command to run the file and some import information. No need to output other conclusion content.
> 2. don't output any other information.
>
> ###### Test
>
> 1. the test mesh for GUI is stored at ../meshes/woody_lo.obj
> 2. **Re-compile the file** until there is a up-to-date runnable exe file. Remind me if you need access to run a terminal command.
>
> ###### Solver:
>
> 1. You should use the following solver for speed up.
>
>    Eigen::VectorXd b;
>      Eigen::VectorXd x;
>      Eigen::SparseMatrix<double> A(m, n);
>
>      // populate b
>      ...
>
>      // populate A
>      std::vector<Eigen::Triplet<double>> A_elem;
>      for (non-zero elements) {
>          A_elem.push_back( Eigen::Triplet<double>(row, column, value) );
>      }
>
>      A.setFromTriplets(A_elem.begin(), A_elem.end());
>
>      // solve using your prefered solver.
>      // For symmetric systems use LDLT or LLT
>      // For non-symmetric systems use PartialPivLU or FullPivLU
>      Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> > solver;
>
>      // Call this only once for each value of A
>      solver.compute(A);
>
>      // Call this as many times as you want for one or different values of b,
>      // as long as A is not changed.
>      x = solver.solve(b);
>
>      // Remember that calling solver.compute() is much more expensive that solver.solve().
>      // Therefore, if the matrix A has not changed, be careful not to call solver.compute()
>      // more than required.
>
> At the end, print: constraint error, residual norm, and final positions.  
> Do not implement ARAP or iterative rotations — only single linear solve.

---



based on the current implementation without rotation, implement rotation following the instructions below:

###### Updata problem:



there is a piece of information that might help. in our basic setting for solver without rotation, it says:
The method should perform a single linear solve step to compute new vertex positions in response to change in handle vertex positions. Initialize the the rotation matrices per vertex using a heuristic of your choice (identity is fine but maybe you can do better?).
not sure whether this can help.

###### Rotation:

1. Once this method works, use it as a basis for a full implementation of ARAP that includes update of the rotation matrices within a loop that uses solve and update steps to reach the desired output.

###### LLM output: 

1. ONLY put the terminal command to run the file.
2. don't output any other information unless it's SUPER important to let me know.

###### GUI:

1. Do not change the current GUI. Only optimize the deforming algorithm.



Solver:

```
Eigen::VectorXd b;
  Eigen::VectorXd x;
  Eigen::SparseMatrix<double> A(m, n);

  // populate b
  ...
  
  // populate A
  std::vector<Eigen::Triplet<double>> A_elem;
  for (non-zero elements) {
      A_elem.push_back( Eigen::Triplet<double>(row, column, value) );
  }

  A.setFromTriplets(A_elem.begin(), A_elem.end());

  // solve using your prefered solver.
  // For symmetric systems use LDLT or LLT
  // For non-symmetric systems use PartialPivLU or FullPivLU
  Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> > solver;

  // Call this only once for each value of A
  solver.compute(A);

  // Call this as many times as you want for one or different values of b,
  // as long as A is not changed.
  x = solver.solve(b);

  // Remember that calling solver.compute() is much more expensive that solver.solve().
  // Therefore, if the matrix A has not changed, be careful not to call solver.compute()
  // more than required.
```


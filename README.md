# Optimization Methods (UMA-035) - MATLAB Lab Implementations

**Thapar Institute of Engineering and Technology, Patiala**  
**School of Mathematics**

This repository contains comprehensive MATLAB implementations for all 9 lab experiments in the Optimization Methods course.

## 📋 Table of Contents

1. [Lab 1 - Graphical Method](#lab-1---graphical-method)
2. [Lab 2 - Basic Feasible Solutions](#lab-2---basic-feasible-solutions)
3. [Lab 3 - Simplex Method](#lab-3---simplex-method)
4. [Lab 4 - Big-M Method](#lab-4---big-m-method)
5. [Lab 5 - Two-Phase Method](#lab-5---two-phase-method)
6. [Lab 6 - Dual Simplex Method](#lab-6---dual-simplex-method)
7. [Lab 7 - Least Cost Method (Transportation)](#lab-7---least-cost-method)
8. [Lab 8 - Multi-Objective LPP](#lab-8---multi-objective-lpp)
9. [Lab 9 - Fibonacci Search](#lab-9---fibonacci-search)

## 🚀 Quick Start

### Running Individual Labs

Each lab can be run independently by executing its MATLAB file:

```matlab
% Example: Run Lab 1
Lab1_Graphical_Method()

% Example: Run Lab 3
Lab3_Simplex_Method()
```

### Running All Labs

To run all labs sequentially:

```matlab
Run_All_Labs()
```

## 📖 Lab Descriptions

### Lab 1 - Graphical Method

**File:** `Lab1_Graphical_Method.m`

**Purpose:** Solve 2-variable Linear Programming Problems using graphical method

**Features:**
- Visual representation of constraints and feasible region
- Automatic detection of corner points
- Support for both maximization and minimization
- Plots objective function contours

**Problems Included:**
1. Max/Min (3x₁ + 2x₂) s.t. 2x₁+4x₂≤8, 3x₁+5x₂≥15
2. Max/Min (3x₁ + 2x₂) s.t. 2x₁+4x₂≥8, 3x₁+5x₂≥15
3. Max/Min (3x₁ + 2x₂) s.t. 2x₁+4x₂≤8, 3x₁+5x₂≤15

**Usage:**
```matlab
Lab1_Graphical_Method()
% Select problem number (1-3)
% Select optimization type (1=Max, 0=Min)
```

### Lab 2 - Basic Feasible Solutions

**File:** `Lab2_Basic_Feasible_Solutions.m`

**Purpose:** Find and classify all Basic Feasible Solutions (BFS) of LPPs

**Features:**
- Systematic enumeration of all possible bases
- Classification: Non-degenerate BFS, Degenerate BFS, Not BFS
- Determinant checking for basic solutions
- Feasibility verification

**Problems Included:**
1. Find all BFS for standard form problem
2. Check if specific variable sets form valid bases
3. Solve LPP by examining all BFS
4. Detect degenerate BFS

**Usage:**
```matlab
Lab2_Basic_Feasible_Solutions()
% Select problem number (1-4)
```

### Lab 3 - Simplex Method

**File:** `Lab3_Simplex_Method.m`

**Purpose:** Solve LPPs with ≤ type constraints using the Simplex Method

**Features:**
- Step-by-step iteration display
- Automatic slack variable handling
- Optimality checking via reduced costs
- Minimum ratio test for pivot selection

**Problems Included:**
1. Max z = x₁ + 2x₂ (2 constraints)
2. Max z = 4x₁ + 6x₂ + 3x₃ + x₄ (3 constraints)
3. Min z with equality constraints

**Usage:**
```matlab
Lab3_Simplex_Method()
% Select problem number (1-3)
```

### Lab 4 - Big-M Method

**File:** `Lab4_BigM_Method.m`

**Purpose:** Solve LPPs with mixed constraints (≤, ≥, =) using Big-M Method

**Features:**
- Automatic artificial variable introduction
- Large penalty coefficient (M = 10⁶)
- Infeasibility detection
- Handles equality and ≥ constraints

**Problems Included:**
1. Min z = 3x₁ + 5x₂ (2 ≥ constraints)
2. Min z = 12x₁ + 10x₂ (3 ≥ constraints)
3. Max z = 3x₁ + 2x₂ (mixed constraints with equality)

**Usage:**
```matlab
Lab4_BigM_Method()
% Select problem number (1-3)
```

### Lab 5 - Two-Phase Method

**File:** `Lab5_TwoPhase_Method.m`

**Purpose:** Solve LPPs using Two-Phase Simplex Method

**Features:**
- **Phase I:** Find initial BFS by minimizing sum of artificial variables
- **Phase II:** Optimize original objective using Phase I solution
- Automatic infeasibility detection
- More efficient than Big-M for some problems

**Problems Included:**
1. Min z = 3x₁ + 5x₂
2. Min z = 12x₁ + 10x₂
3. Max z = 3x₁ + 2x₂ (with equality constraint)

**Usage:**
```matlab
Lab5_TwoPhase_Method()
% Select problem number (1-3)
```

### Lab 6 - Dual Simplex Method

**File:** `Lab6_Dual_Simplex_Method.m`

**Purpose:** Solve LPPs that are initially optimal but infeasible

**Features:**
- Starts with optimal but infeasible solution
- Maintains optimality while gaining feasibility
- Useful for sensitivity analysis
- Minimum ratio test on dual variables

**Problems Included:**
1. Min z = 3x₁ + 5x₂
2. Min z = 12x₁ + 10x₂
3. Min z = 3x₁ + 2x₂ (mixed constraints)

**Usage:**
```matlab
Lab6_Dual_Simplex_Method()
% Select problem number (1-3)
```

### Lab 7 - Least Cost Method

**File:** `Lab7_Least_Cost_Method.m`

**Purpose:** Find initial BFS for Transportation Problems using Least Cost Method

**Features:**
- Cost-based allocation strategy
- Automatic balancing for unbalanced problems
- Step-by-step allocation display
- Total transportation cost calculation

**Problems Included:**
1. 3 Sources × 4 Destinations
2. 4 Sources × 5 Destinations

**Usage:**
```matlab
Lab7_Least_Cost_Method()
% Select problem number (1-2)
```

### Lab 8 - Multi-Objective LPP

**File:** `Lab8_MultiObjective_LPP.m`

**Purpose:** Solve Multi-Objective Linear Programming using Weighted Sum Method

**Features:**
- Combines multiple objectives into single objective
- Equal weighting strategy (arithmetic mean)
- Reports values for all individual objectives
- Uses Big-M method for solving

**Problems Included:**
1. Two objectives with 3 variables
2. Two objectives with 3 variables (different constraints)

**Usage:**
```matlab
Lab8_MultiObjective_LPP()
% Select problem number (1-2)
```

### Lab 9 - Fibonacci Search

**File:** `Lab9_Fibonacci_Search.m`

**Purpose:** Find minimum/maximum of univariate functions using Fibonacci Search

**Features:**
- Fibonacci number generation
- Automatic interval reduction
- Visual plotting of objective function and optimal point
- Achieves target interval of uncertainty

**Example Problem:**
- Minimize f(x) = x(x-2) on [0, 1.5]
- Target: Final interval = 0.25 × Initial interval

**Usage:**
```matlab
Lab9_Fibonacci_Search()
```

## 📊 Algorithm Implementations

### Key Algorithms Implemented:

1. **Simplex Algorithm:** Standard tableau method with pivot operations
2. **Big-M Method:** Artificial variables with penalty coefficients
3. **Two-Phase Method:** 
   - Phase I: Minimize artificial variables
   - Phase II: Optimize original objective
4. **Dual Simplex:** Maintains optimality, gains feasibility
5. **Least Cost:** Greedy allocation based on minimum cost cells
6. **Fibonacci Search:** Golden-ratio-based interval reduction

## 🔧 Code Features

### Common Features Across All Labs:

- **Interactive:** User-friendly menu-driven interface
- **Verbose Output:** Detailed iteration-by-iteration progress
- **Error Handling:** Checks for infeasibility, unboundedness, singularity
- **Validation:** Automatic constraint and solution verification
- **Visualization:** Graphical outputs where applicable (Labs 1, 9)

### Code Quality:

- Well-commented and documented
- Follows MATLAB best practices
- Modular design with reusable functions
- Clear variable naming conventions
- Comprehensive test problems

## 📝 Input/Output Format

### Typical Input Format:

```matlab
% For constraints Ax ≤ b:
A = [coefficient matrix];
b = [right-hand side vector];
c = [objective coefficients];
basis = [initial basic variables];
```

### Typical Output:

```
========================================
OPTIMAL SOLUTION
========================================
x1 = 1.5000
x2 = 0.5000
Optimal value: z = 7.5000
Iterations: 3
```

## 🎯 Learning Objectives

By working through these labs, students will:

1. Understand fundamental LP solution methods
2. Implement optimization algorithms from scratch
3. Recognize when to use different methods
4. Interpret algorithm iterations and convergence
5. Handle special cases (infeasibility, unboundedness, degeneracy)

## 📚 References

- Course: Optimization Methods (UMA-035)
- Institution: Thapar Institute of Engineering and Technology, Patiala
- School: School of Mathematics

## 🤝 Usage Guidelines

### For Students:

1. Run each lab before the corresponding class
2. Try different problem variations
3. Observe iteration details to understand the algorithm
4. Experiment with your own problem instances

### For Instructors:

- Use as teaching demonstrations
- Assign as programming exercises
- Modify problems for assessments
- Extend for advanced topics

## ⚠️ Important Notes

1. **MATLAB Version:** Tested on MATLAB R2019b and later
2. **Numerical Precision:** Uses tolerance of 10⁻⁶ for comparisons
3. **Big-M Value:** Set to 10⁶ (adjust if needed for specific problems)
4. **Maximum Iterations:** Default 100 (prevents infinite loops)

## 🐛 Troubleshooting

**Issue:** "Singular matrix" error  
**Solution:** Check constraint linear independence; problem may be degenerate

**Issue:** "Maximum iterations reached"  
**Solution:** Increase max_iterations or check problem formulation

**Issue:** Unexpected infeasibility  
**Solution:** Verify constraint equations and signs

## 📧 Support

For issues or questions related to these implementations:
- Review the inline comments in each file
- Check the algorithm descriptions in the lab manual
- Verify input data format

## 📄 License

Educational use only. Developed for Optimization Methods (UMA-035) course.

---

**Version:** 1.0  
**Last Updated:** January 2026  
**Author:** Course Material - Thapar Institute

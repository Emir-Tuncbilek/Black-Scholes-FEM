# 1D Finite Element Solver  
### FEM Spatial Discretization + Crank–Nicolson Time Integration

This repository contains a **1D time-dependent numerical solver** implemented in **MATLAB**.  
The solution is obtained using the **Finite Element Method (FEM)** for spatial discretization
(with **P1 and P2 Lagrange elements**) and the **Crank–Nicolson finite-difference scheme**
for time integration.

> **Important note on provenance**  
> This codebase is **derived from a 1D FEM solver framework originally provided by  
> Prof. Serge Prudhomme (Polytechnique Montréal)** in the context of the course  
> **MTH8207 – Finite Element Methods**.
>
> The original solver served as a pedagogical foundation.  
> The present repository **extends and modifies** this framework substantially by introducing:
> - time-dependent formulations,
> - Crank–Nicolson time stepping,
> - Method of Manufactured Solutions (MMS),
> - error estimation in \(L^2\) and \(H^1\) norms,
> - application to the Black–Scholes equation,
> - and quantity-of-interest (QoI) convergence studies.

---

## Numerical Methods

- **Spatial discretization:** Finite Element Method (FEM)
- **Element types:** Lagrange P1 (linear), P2 (quadratic)
- **Time discretization:** Crank–Nicolson scheme
- **Verification:** Method of Manufactured Solutions (MMS)
- **Error norms:** \(L^2\), \(H^1\)
- **Applications:**  
  - Time-dependent PDEs  
  - Black–Scholes equation  
  - Quantity-of-Interest (QoI) convergence analysis

---

## Academic Context

This project was developed as part of the course:

- **MTH8207 – Finite Element Methods**  
- **Instructor:** Prof. Serge Prudhomme  
- **Institution:** Polytechnique Montréal (PolyMTL)

The base 1D FEM solver was provided by the instructor for instructional purposes.  
All additional features—including time stepping, boundary-condition handling,
manufactured solutions, convergence studies, and QoI estimation—were implemented
by the authors.

---

## Authors

- **Samatar Aberkane**
- **Emir Tuncbilek**

---

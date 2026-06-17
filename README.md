# Finite Difference Method for Hyperbolic Equations

A numerical implementation of an explicit finite-difference method for the one-dimensional wave equation. The repository derives the centered finite-difference scheme, expresses the time-stepping recurrence in matrix form, improves the first time-step approximation, and compares the numerical result with an analytical solution.

## Problem statement

The method solves the one-dimensional wave equation

$$\frac{\partial^2 u}{\partial t^2}(x,t)-\alpha^2\frac{\partial^2 u}{\partial x^2}(x,t)=0,\qquad 0<x<l,\quad t>0,$$

subject to homogeneous Dirichlet boundary conditions

$$u(0,t)=u(l,t)=0,$$

and the initial conditions

$$u(x,0)=f(x),\qquad \frac{\partial u}{\partial t}(x,0)=g(x).$$

The spatial and temporal grids are defined by

$$x_i=ih,\qquad t_j=jk,\qquad h=\frac{l}{m},\qquad k=\frac{T}{N}.$$

## Numerical method

Using centered differences in both space and time gives

$$\frac{w_{i,j+1}-2w_{i,j}+w_{i,j-1}}{k^2}-\alpha^2\frac{w_{i+1,j}-2w_{i,j}+w_{i-1,j}}{h^2}=0.$$

With the Courant number

$$\lambda=\frac{\alpha k}{h},$$

the explicit recurrence becomes

$$w_{i,j+1}=2(1-\lambda^2)w_{i,j}+\lambda^2\left(w_{i+1,j}+w_{i-1,j}\right)-w_{i,j-1}.$$

The first time level is initialized with the second-order approximation

$$w_{i,1}=(1-\lambda^2)f(x_i)+\frac{\lambda^2}{2}\left[f(x_{i+1})+f(x_{i-1})\right]+kg(x_i).$$

For the standard explicit wave-equation scheme, stability requires the Courant–Friedrichs–Lewy condition

$$\lambda=\frac{\alpha k}{h}\leq 1.$$

## Matrix formulation

Let

$$\mathbf{w}^{(j)}=\begin{bmatrix}w_{1,j}&w_{2,j}&\cdots&w_{m-1,j}\end{bmatrix}^{\mathsf T}.$$

The update can be written as

$$\mathbf{w}^{(j+1)}=A\mathbf{w}^{(j)}-\mathbf{w}^{(j-1)},$$

where the tridiagonal matrix is represented below using an HTML table so that it renders consistently on GitHub.

<p align="center"><strong>A =</strong></p>

<table align="center">
  <tbody>
    <tr>
      <td>2(1 − λ<sup>2</sup>)</td>
      <td>λ<sup>2</sup></td>
      <td>0</td>
      <td>⋯</td>
      <td>0</td>
    </tr>
    <tr>
      <td>λ<sup>2</sup></td>
      <td>2(1 − λ<sup>2</sup>)</td>
      <td>λ<sup>2</sup></td>
      <td>⋱</td>
      <td>⋮</td>
    </tr>
    <tr>
      <td>0</td>
      <td>⋱</td>
      <td>⋱</td>
      <td>⋱</td>
      <td>0</td>
    </tr>
    <tr>
      <td>⋮</td>
      <td>⋱</td>
      <td>λ<sup>2</sup></td>
      <td>2(1 − λ<sup>2</sup>)</td>
      <td>λ<sup>2</sup></td>
    </tr>
    <tr>
      <td>0</td>
      <td>⋯</td>
      <td>0</td>
      <td>λ<sup>2</sup></td>
      <td>2(1 − λ<sup>2</sup>)</td>
    </tr>
  </tbody>
</table>

## Included example

The Python script evaluates the problem with

- $l=1$
- $T=1$
- $\alpha=2$
- $m=10$, so $h=0.1$
- $N=20$, so $k=0.05$
- $f(x)=\sin(\pi x)$
- $g(x)=0$

The analytical solution is

$$u(x,t)=\sin(\pi x)\cos(2\pi t).$$

For these parameters, $\lambda=1$. The selected spatial eigenmode is propagated exactly by the discrete dispersion relation, up to floating-point roundoff, so the error at $t=1$ is approximately machine precision.

## Repository structure

```text
.
├── Backward_Difference_Hyperbolic.tex       # Original Spanish report
├── Finite_Difference_Method_Hyperbolic_EN.tex # Corrected English report
├── finite_difference_method_hyperbolic.py   # NumPy/Pandas implementation
└── README.md
```

## Requirements

- Python 3.8 or newer
- NumPy
- pandas

Install the dependencies with

```bash
python -m pip install numpy pandas
```

## Usage

Run the numerical example from the repository root:

```bash
python finite_difference_method_hyperbolic.py
```

The script prints the numerical solution at each time level and creates

```text
error_table.csv
```

containing the grid points, the numerical approximation at $t=1$, the analytical solution, and the absolute error.

## Main functions

### `finite_difference_method_hyperbolic`

```python
finite_difference_method_hyperbolic(l, T, alpha, N, m, f, g)
```

Constructs the tridiagonal update matrix, initializes the first two time levels, and advances the solution through the temporal grid.

### `error_table`

```python
error_table(m, x, w, u)
```

Compares the final numerical vector with the analytical solution and exports the results to CSV.

## Notes

The source filename refers to a “backward difference,” but the implemented algorithm is an explicit centered-difference scheme: it uses a centered second derivative in space and a centered second derivative in time. The English report uses terminology consistent with the actual method.

**Author:** BSc. Julio Medina

## Reference

Richard L. Burden and J. Douglas Faires, *Numerical Analysis*, 9th edition, Brooks/Cole, Cengage Learning.

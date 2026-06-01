# Microplastic Simulation code

Simulation code for the Master's Thesis:
**Mathematical Modeling of Microplastics Spread  in the Marine Environment**
$\def\ref{\text{ref}}$

## Note

### Normalization

As the PDE system is stiff, the simulation code normalizes the system. 
For this we define the referenced parameters,
$$
\begin{align*}
\tilde{P}_m &= \frac{P_m}{m_\ref} \\
\tilde{P}_s &= \frac{P_s}{m_\ref} \\
\tilde{P}_l &= \frac{P_l}{m_\ref}
\end{align*}
$$
and 
$$
\begin{align*}
\tilde{N}_s = \frac{N_s}{{N}_{s, \ref}} \\
\tilde{N}_l = \frac{N_l}{{N}_{l, \ref}}.
\end{align*}
$$

Inserting these values into the system, yields

$$
\begin{align*}
\frac{\partial \tilde{P}_m}{\partial t} &= D_m \Delta \tilde{P}_m - (c_s \cdot {N}_{s, \ref}) \tilde{N}_s \tilde{P}_m - (c_l \cdot {N}_{l, \ref}) \tilde{N}_l \tilde{P}_m + d_s \tilde{P}_s + d_l \tilde{P}_l \\
\frac{\partial \tilde{P}_s}{\partial t} &= D_s \Delta \tilde{P}_s + (c_s \cdot {N}_{s, \ref}) \tilde{N}_s \tilde{P}_m - d_s \tilde{P}_s - (\beta \cdot {N}_{l, \ref}) \tilde{N}_l \tilde{P}_s \\
\frac{\partial \tilde{P}_l}{\partial t} &= D_l \Delta \tilde{P}_l + (c_l \cdot {N}_{l, \ref}) \tilde{N}_l \tilde{P}_m - d_l \tilde{P}_l + (\beta \cdot {N}_{l, \ref}) \tilde{N}_l \tilde{P}_s \\
\frac{\partial \tilde{N}_s}{\partial t} &= D_s \Delta \tilde{N}_s \\
\frac{\partial \tilde{N}_l}{\partial t} &= D_l \Delta \tilde{N}_l. 
\end{align*}
$$

This normalized system then gets solved.
As the reference values, the peak concentration of the Gaussian is used.

### FiPy Behavior

To solve the complete system implicitly, FiPy requires the terms to be inside `Implicit...` Terms, otherwise 
terms would be treated explicitly. When solving implicitly, the terms must also be **added** to the 
complete equation and if the term is negative, the sign must be inside the `Implicit...` term to have correct behavior.
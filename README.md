# Microplastic Simulation code

Simulation code for the Master's Thesis:
**Mathematical Modeling of Microplastics Spread  in the Marine Environment**

## Note

### Normalization

As the PDE system is stiff, the simulation code normalizes the system. 
For this we define the referenced parameters,

![equation](https://latex.codecogs.com/svg.image?\large&space;\begin{align*}\tilde{P}_m&space;&=\frac{P_m}{m_\text{ref}}\\&space;\tilde{P}_s&space;&=\frac{P_s}{m_\text{ref}}\\&space;\tilde{P}_l&space;&=\frac{P_l}{m_\text{ref}}\end{align*})

and 

![equation](https://latex.codecogs.com/svg.image?\large&space;\begin{align*}\tilde{N}_s=\frac{N_s}{{N}_{s,\text{ref}}}\\&space;\tilde{N}_l=\frac{N_l}{{N}_{l,\text{ref}}}.&space;\end{align*})

Inserting these values into the system, yields

![equation](https://latex.codecogs.com/svg.image?\large&space;\begin{align*}\frac{\partial&space;\tilde{P}_m}{\partial&space;t}&=D_m&space;\Delta&space;\tilde{P}_m-(c_s&space;\cdot{N}_{s,\text{ref}})\tilde{N}_s&space;\tilde{P}_m-(c_l&space;\cdot{N}_{l,\text{ref}})\tilde{N}_l&space;\tilde{P}_m&plus;d_s&space;\tilde{P}_s&plus;d_l&space;\tilde{P}_l&space;\\&space;\frac{\partial&space;\tilde{P}_s}{\partial&space;t}&=D_s&space;\Delta&space;\tilde{P}_s&plus;(c_s&space;\cdot{N}_{s,\text{ref}})\tilde{N}_s&space;\tilde{P}_m-d_s&space;\tilde{P}_s-(\beta&space;\cdot{N}_{l,\text{ref}})\tilde{N}_l&space;\tilde{P}_s&space;\\&space;\frac{\partial&space;\tilde{P}_l}{\partial&space;t}&=D_l&space;\Delta&space;\tilde{P}_l&plus;(c_l&space;\cdot{N}_{l,\text{ref}})\tilde{N}_l&space;\tilde{P}_m-d_l&space;\tilde{P}_l&plus;(\beta&space;\cdot{N}_{l,\text{ref}})\tilde{N}_l&space;\tilde{P}_s&space;\\&space;\frac{\partial&space;\tilde{N}_s}{\partial&space;t}&=D_s&space;\Delta&space;\tilde{N}_s&space;\\&space;\frac{\partial&space;\tilde{N}_l}{\partial&space;t}&=D_l&space;\Delta&space;\tilde{N}_l.&space;&space;\end{align*})

This normalized system then gets solved.
As the reference values, the peak concentration of the Gaussian is used.

### FiPy Behavior

To solve the complete system implicitly, FiPy requires the terms to be inside `Implicit...` Terms, otherwise 
terms would be treated explicitly. When solving implicitly, the terms must also be **added** to the 
complete equation and if the term is negative, the sign must be inside the `Implicit...` term to have correct behavior.
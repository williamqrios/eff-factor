### Effectiveness factors for reversible reactions
1. Required data
    - Effective diffusivity of the components in the fluid mixture
    - Catalyst properties ($\varepsilon$, $\tau$, $L$ or $R_p$, $\rho_p$)
    - Thiele modulus for the irreversible reaction (can be computed inside the function based on geometry data)
    - Equilibrium concentration of the limiting reactant $C_{A,eq}$ 
    - Concentration of the components in the mixture 
    $$ C_A, C_B, C_C, ... $$
    - Equilibrium constant in terms of concentrations
    $$ K_c = \prod_{i=1}^n C_{i, \mathrm{eq}}^{\nu_i} $$
    - Kinetic constant, $k$

$$ \phi_g = \frac{\ell \sqrt{ \frac{\rho_p}{D_{\mathrm{ef, A}}}   } r_s  }{ \sqrt{2 k} \ F} = \frac{\ell \sqrt{k \frac{\rho_p }{D_{\mathrm{ef, A}}}   } r_s  }{ \sqrt{2}\ k \ F} = \frac{\phi\  r_s}{ \sqrt{2}\ k \ F} $$ 

where $\ell$ is a characteristic dimension.

2. Computation based on catalyst geometry (sphere, cylinder, or semi-infinite slab)
$$ \eta = \frac{\tanh \phi_g}{\phi_g} \quad \mathrm{slab} $$
$$ \eta = \frac{3}{\phi_g} \left( \frac{1}{\tanh \phi_g} - \frac{1}{\phi_g} \right) \quad \mathrm{sphere} $$
$$ \eta = \frac{2}{\phi_g} \frac{I_1 (\phi_g)}{I_0 (\phi_g)} \quad \mathrm{cylinder} $$
where $I_n$ ($n = 0, 1$) are modified Bessel functions of the first kind. 

For $n=0$ the Bessel function can be approximated as the series:
$$ I_0(z) = \sum_{k=0}^\infty \frac{\left(\frac{1}{4}  z^2\right)^k}{\left(k!\right)^2} $$
For a real number $n$ the function is generically given by:
$$ I_n(z) = \left( \frac{1}{2} z \right)^n\  \sum_{k=0}^\infty \frac{\left(\frac{1}{4}  z^2\right)^k}{k!\ \Gamma (n + k + 1) }   $$
where $\Gamma (n)$ is the gamma function, which is an extension of the factorial to complex and real numbers, defined as:
$$ \Gamma(n) = (n-1)! $$ 
$$ \Gamma(n+k+1)\big|_{n=1}  =  \Gamma(k+2) = (k+2-1)!=(k+1)!$$
So the $I_1$ function is:
$$ I_1(z) = \frac{1}{2} z \  \sum_{k=0}^\infty \frac{\left(\frac{1}{4}  z^2\right)^k}{k!\ (k+1)! }   $$



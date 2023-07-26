# Electronic-Relaxation-Coordinates

This code reads the energy of two electronic states as a function of time in an excited state simulation, 
along with the coordinates of the nuclei.
In order to elucidate the nuclear rearrangements involved in the S$_1 \to$ S$_0$ relaxation, this code performs a variance analysis combining the nuclear 
coordinate fluctuations ($x-\overline{x}$) along the ab initio excited state trajectories and the energy difference between the diabatic 
electronic states ($\Delta E^D$).  
If $N$ is the number of atoms in the system, we define the $3N-$dimensional S$_1 \to$ S$_0$ relaxation pathway vector $\mathbf{c}$ as:
 
\begin{equation}
\label{eq:PCA1}
\mathbf{c_i} = \dfrac{\langle[x_i(t)-\overline{x_i}]\text{Sign}[-\Delta E^D(t)]\text{exp}^{-\frac{|\Delta E^D(t)|}{\alpha kT}} \rangle}{\sqrt{\langle [x_i(t)-\overline{x_i}]^2 \rangle \langle [\text{exp}^{-\frac{|\Delta E^D(t)|}{\alpha kT}}]^2} \rangle},
\end{equation}

where the index $i$ spans over the $3N$ Cartesian coordinates of the system, 
$x_i(t)$ represents the $i$-th component of the cartesian position vector at time $t$, 
$\Delta E^D$ is the diabatic energy difference between S$_0$ and S$_1$, 
and the angular brackets as well as the over-bar represent a time average. 
The diabatic energies were approximated by the adiabatic ones before the CoIn crossing and swapping the S$_1$ and S$_0$ 
identities after the CoIn passage.  The second term in the numerator, $\text{Sign}[-\Delta E^D(t)]$, sets the direction of the relaxation 
pathway vector from the S$_1$ to the S$_0$ configurations. The third term, $\text{exp}^{-\frac{|\Delta E^D(t)|}{\alpha kT}}$, 
is an Arrhenius-like factor, where T is the room temperature, 
$k$ is the Boltzmann constant, and $\alpha$ is an adjustable parameter that controls the width of the exponential term with respect to 
$\Delta E^D(t)$ (throughout this work $\alpha$ was fixed to $100$, see Figure S10 in the SI). 
This factor enables disentangling the thermal fluctuations from the relaxation process.  
As the nuclear configurations get closer to the CoIn, the $|\Delta E^D(t)|$ tends to zero and the Arrhenius-like term peaks for these 
configurations, increasing their weight in the ensemble average. 
In this way, the vector $\mathbf{c}$ is a linear estimator of the nuclear fluctuations in the S$_1$ state that lead to the S$_1$-S$_0$ crossing. 
It is important to note that usually the excited state landscape is characterized by many accessible CoIns, 
and several different decay pathways can be accessible. 
In these cases, the relaxation pathway vector $\mathbf{c}$ represents a statistical average of all the decay motions.

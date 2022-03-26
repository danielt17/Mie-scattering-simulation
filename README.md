# Mie-scattering-simulation
A simulation of mie scattering from a homogenous sphere written in python


Mie scattering:
Mie scattering/solution is a solution to Maxwell's equations which describes the scattering of an electromagnetic plane wave by homogeneous sphere. The solution takes the form of an infinite series of spherical multipole partial waves.  The Mie scattering formulation is most useful when the size of the scattering particles is comparable to the wavelength of light, rather than much smaller or much larger.
Mathematics of Mie Scattering:
We start by considering the following setup, we have spherical nanoparticle (about the same size of our wavelength), which want to scatter an EM wave from. We consider a plane wave propagating along the z-axis polarized along the x-axis. Dielectric and magnetic permeabilities of a particle are ϵ_1 and μ_1 while the environment has constant are ϵ and μ, Figure 1 describes the experimental setup.
 
Figure 1 - Scattering of the plane wave incidence direction is parallel to the z-axis, polarization is parallel to the x-axis, nanoparticle's radius is a.
In  order to solve the scattering problem described here, we write first the solutions of the vector Helmholtz equation in spherical coordinates, since the fields inside and outside the particles must satisfy it.
We first derive the Helmholtz equations from maxwell equations.
∇⋅D=ρ_F
∇×E+∂B/∂t=0
∇⋅B=0
∇×H=J_F+∂D/∂t
Where E is the electric field B the magnetic induction. The electric displacement D and magnetic field H are defined by
D=ϵ_0 E+P
H=B/μ_0 -M
Where P is the electric polarization (average electric dipole moment per unit volume), M the magnetization (average magnetic dipole per unit volume), ϵ_0 the permittivity, and μ_0 the permeability of free space. The charge density ρ_F and the current density J_F are associated with so-called "free" charges. The terms "free" and "bound" re sometimes set in quotation marks, which indicates that they are slightly suspect. In addition, the equations above are not sufficient in themselves they must be supplemented with constitutive relations, which are assumed to have the form
J_F=σE
B=μH
P=ϵ_0 χE
Where σ is the conductivity, μ the permeability, and χ the electric susceptibility. The phenomenological coefficients σ,μ,χ depend on the medium under consideration, but will be assumed to be independent of the fields (thus, the medium is linear), independent of position (the medium is homogeneous), and independent of direction (the medium is isotropic). There are generalization which deal with many classes of materials for which these assumptions are not valid. Thus, the free relations above are not universal laws of nature, but merely describe a particular class of materials which, fortunately, has a large number of members.
We will go to frequency domain now, as it is more useful in solving maxwells equations. One can assume that all fields in maxwells equation have time dependence e^(-iωt), therefore we obtain
∇⋅D=ρ_F→∇⋅(ϵ_0 E+P)=ρ_F→∇⋅(ϵ_0 E+ϵ_0 χE)=ρ_F
∇×E+∂B/∂t=0→∇×E+μ ∂H/∂t=0
∇⋅B=0→1/μ ∇⋅H=0→∇⋅H=0
∇×H=J_F+∂D/∂t→∇×H=σE+∂(ϵ_0 E+P)/∂t→∇×H=σE+∂(ϵ_0 E+ϵ_0 χE)/∂t
Assuming free space we obtain (ρ_F=0)
∇⋅(ϵ_0 (1+χ)E)=0→∇⋅(ϵ_0 (1+χ) E_c e^(-iωt) )=0→∇⋅(ϵ_0 (1+χ) E_c )=0
∇×E+1/μ  ∂H/∂t=0→∇×E_c e^(-iωt)-iωμ (∂H_c)/∂t e^(-iωt)=0→∇×E_c=iωμH_c
∇⋅H=0→∇⋅H_c e^(-iωt)=0→∇⋅H_c=0
∇×H=σE+∂(ϵ_0 (1+χ)E)/∂t→∇×H_c e^(-iωt)=σE_c e^(-iωt)-iω(ϵ_0 (1+χ) E_c ) e^(-iωt)




After some work we get:
∇⋅(ϵE_c )=0
∇×E_c=iωμH_c
∇⋅H_c=0
∇×H_c=(σ-iω(ϵ_0 (1+χ))) E_c→∇×H_c=-iω(ϵ_0 (1+χ)+iσ/ω) E_c
Finally
∇⋅(ϵE_c )=0
∇×E_c=iωμH_c
∇⋅H_c=0
∇×H_c=-iωϵE_c
Where the complex permittivity is
ϵ=ϵ_0 (1+χ)+iσ/ω
From now on we switch the letters with subscript c to not have it and derive the Helmholtz equation (and we use the Fourier transform to derive this)
∇⋅E=0
∇⋅H=0
∇×E=iωμH
∇×H=-iωϵE
At all points where ϵ and μ are continuous. The curl of the equation results in:
∇×(∇×E)=∇×(iωμH)=ω^2 ϵμE
∇×(∇×H)=-∇×(iωϵE)=ω^2 ϵμH
We use the vector identity:
∇×(∇×A)=∇(∇⋅A)-∇ ⋅(∇A)→using zero divergance=∇⋅(∇A)
We get:
∇^2 E+k^2 E=0,    ∇^2 H+k^2 H=0
Where k^2=ω^2 ϵμ and ∇^2 A=∇⋅(∇A).
Now we have to remember that inorder to solve this setup we must satisfy divergence free vector fields, and the fact that E and H are not independent. 
To solve this we start by assuming we can construct a vector function M such that it is built using a scalar function ψ and an arbitrary constant vector c.
M=∇×(cψ)  
The divergence of the curl of any vector function vanishes (theorem):
∇⋅M=∇⋅∇×(cψ)=0
one can show after some work that solving the equations above is the same as solving:
∇^2 M+k^2 M=∇×[c(∇^2 ψ+k^2 ψ)]
If ψ is a solution to the scalar wave equation
∇^2 ψ+k^2 ψ=0
One may write (M is orthogonal to c)
M=-c×∇ψ
Finally after constructing another perpendicular vector N which is defined as
N=(∇×M)/k
We also have
∇^2 N+k^2 N=0,∇×N=kM
We work in spherical coordinates as the problem has spherical symmetry therefore we may write the arbitrary vector c as r.
M=∇×(rψ)
We need to solve the scalar wave equation in spherical polar coordinates:
1/r^2   ∂/∂r (r^2  ∂ψ/∂r)+1/(r^2 sinθ)  ∂/∂θ (sinθ ∂ψ/∂θ)+1/(r^2 sinθ)  (∂^2 ψ)/(∂ϕ^2 )+k^2 ψ=0 
The solution to this equation is spherical harmonics.
ψ_emn=cos⁡(mϕ) P_n^m (cosθ) z_n (kr)
ψ_emn=sin⁡(mϕ) P_n^m (cosθ) z_n (kr)
Where P_n^m (x) are the associated Legendre polynomials and z_n (x) are spherical Bessel functions.
Using the relations written below we get M,N 
M=∇×(rψ),N=(∇×M)/k






Finally, the full solution of scattering types:
Incident plane wave (superscript (1) means that in the radial part of the function ψ are spherical Bessel functions of the first kind):
E_inc=E_0 ∑_(n=1)^∞▒〖i^n⋅(2n+1)/n(n+1)  (M_o1n^((1) ) (k,r)-iN_e1n^((1) ) (k,r)) 〗
H_inc=-k/ωμ E_0 ∑_(n=1)^∞▒〖i^n⋅(2n+1)/n(n+1)  (M_e1n^((1) ) (k,r)+iN_o1n^((1) ) (k,r)) 〗
The scattered fields (superscript (3) means that in the radial part of the function ψ are spherical Hankel functions of the first kind):
E_s=E_0 ∑_(n=1)^∞▒〖i^n⋅(2n+1)/n(n+1)  (ia_n N_e1n^((3) ) (k,r)-b_n M_o1n^((3) ) (k,r)) 〗
H_s=k/ωμ E_0 ∑_(n=1)^∞▒〖i^n⋅(2n+1)/n(n+1)  (a_n M_e1n^((3) ) (k,r)+ib_n N_o1n^((3) ) (k,r)) 〗
The internal fields:
E_int=E_0 ∑_(n=1)^∞▒〖i^n⋅(2n+1)/n(n+1)  (-id_n N_e1n^((1) ) (k_1,r)+c_n M_o1n^((1) ) (k_1,r)) 〗
H_int=-k_1/(ωμ_1 ) E_0 ∑_(n=1)^∞▒〖i^n⋅(2n+1)/n(n+1)  (-id_n N_e1n^((1) ) (k_1,r)+c_n M_o1n^((1) ) (k_1,r)) 〗
Where k=ω/c n and k_1=ω/c n_1 where the former is the wave vector outside the particle and the latter is the wave vector inside the particle. n and n_1 are the refractive indices of the medium and the particle. Applying the boundary conditions and finite field inside the center of the particle we get.
a_n (ω)=(μn_1^2 [ρj_n (ρ)]^' j_n (ρ_1 )-μ_1 n^2 [ρ_1 j_n (ρ_1 )]^' j_n (ρ))/(μn_1^2 [ρj_n (ρ)]^' j_n (ρ_1 )-μ_1 n^2 [ρ_1 j_n (ρ_1 )]^' h_n (ρ) )
b_n (ω)=(μ_1 [ρj_n (ρ)]^' j_n (ρ_1 )-μ[ρ_1 j_n (ρ_1 )]^' j_n (ρ))/(μ_1 [ρj_n (ρ)]^' j_n (ρ_1 )-μ[ρ_1 j_n (ρ_1 )]^' h_n (ρ) )
c_n (ω)=(μ_1 [ρh_n (ρ)]^' j_n (ρ)-μ_1 [ρ_1 j_n (ρ_1 )]^' h_n (ρ))/(μ_1 [ρh_n (ρ)]^' j_n (ρ_1 )-μ[ρ_1 j_n (ρ_1 )]^' h_n (ρ) )
d_n (ω)=(μ_1 n_1 n[ρh_n (ρ)]^' j_n (ρ)-μ_1 n_1 n[ρj_n (ρ)]^' h_n (ρ))/(μ_1 n_1^2 [ρh_n (ρ)]^' j_n (ρ_1 )-μ_1 n_1^2 [ρ_1 j_n (ρ_1 )]^' h_n (ρ) )
Where ρ=ka and ρ=k_1 a, and a is the radius of the sphere. j_n and h_n represent the spherical functions of Bessel and Hankel of the first kind, respectively.
This is the complete solution of Mie scattering problem for a homogenous sphere.

We don’t measure the EM waves, but measure different values of cross sections of which we have three:  Extinction (Scattered + absorbed) Q_e, scattering Q_s and absorption Q_a. Each of the efficiencies are defined by their respective cross sections Q_i=σ_i/(πa^2 ).
Each of them is defined in the following way:
W_i=-∫_A▒〖S_i⋅e ̂_r  dA〗→C_i=W_i/I_i 
Where S_i is the poynting vector.
C_sca=2π/k^2  ∑_(n=1)^∞▒(2n+1)(|a_n |^2+|b_n |^2 ) 
C_ext=2π/k^2  ∑_(n=1)^∞▒(2n+1)Re(a_n+b_n ) 


![Mie scattering](https://user-images.githubusercontent.com/60748408/160252081-679808dd-4a3c-4e34-8ee5-5e46f70c17a0.png)


![Asymmetry coefficent](https://user-images.githubusercontent.com/60748408/160252082-92271390-8dfe-4def-8f79-ce9726f60b4c.png)

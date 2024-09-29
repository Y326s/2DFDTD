# 2DFDTD
An electric current source in the strike direction
$\bold{J}=\hat{y}J_y=\hat{y}\delta(r-r_s)s(t)$
, where 

$$
s(t)=\left\{
\begin{array}{ll}
\sum_{n=0}^3 a_n cos(2n\pi t/T) &\text{if } 0\leq t \leq T, \\ 
0 &\text{otherwise},
\end{array} 
\right.
$$
and $a_0=0.35322222$, $a_1=-0.448$, $a_2=0.145$, $a_3=-0.01022222$.

PEC condition for all the boundaries.
## Case I: Homogeneous media (air) 

<p>Run <kbd>main_homogeneous.m</kbd>.</p>

## Case II: 3-layer media

<p>Run <kbd>main_layer.m</kbd>.</p>

## Case III: Homogeneous media (air) with irregular shape boundaries inside the conputational domain

<p>**See mask.png for the visualized default boundary (mask).</p>
<p>Run <kbd>main_homogeneous_maskVer.m</kbd> first, then run <kbd>plotEz.m</kbd> to visualize the wave propagation.</p>
<p></p>



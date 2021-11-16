# CrowLongWavelengthTheory
## CrowTheory.m

A MATLAB code to calculate the growth rate predicted by Crow's theory. 

The self-induction function and its asymptotic format proposed in this theory valid only for $ka\ll 1$, and lead to the prediction of spurious unstable wavelength (short-wave) band out of this range of validity.

## CalSelfInductionFun.m

A MATLAB code to calculate the self-induction function (self-induction angular frequency) by exact method, crow theory, long-wavelength asymptotic method and numerical fitting. 

### The dispersion relation (exact method)

$$
\frac{1}{\beta a}\frac{J_m^{\prime}(\beta a)}{J_m(\beta a)}=-\frac{K_m^{\prime}(ka)}{kaK_m(ka)}-\frac{sm\sqrt{(\beta a)^2+(ka)^2}}{ka(\beta a)^2}
$$

The expression at $ka = 0$ has singularity, but its solution can be obtained by its limiting behavior.

> From the properties of the modified Bessel function, the first term at the right-hand side, $-\frac{K_m^{\prime}(ka)}{kaK_m(ka)}$, is a positive monotonically decreasing function of $ka$, decreasing from $\infin$ to $0$ as $ka$ increases from $0$ to $\infin$.
>
> ...
>
> There are an infinite number of roots, both retrograde ($s=1$) and co-grade ($s=-1$). It follows from the expansions for small $ka$ and $\beta a$ that there is no co-grade root for small $\beta a$. For $ka \rightarrow \infin$, $\beta a$ is found to be bounded $\sigma \rightarrow (2s-m)\Omega$.
>
> As $ka \rightarrow 0$, $\beta a \rightarrow j_{mn}$ (the $n$th root of $J_m(x)=0$) and $\sigma \rightarrow -\Omega$, except for the smallest retrograde root.
>
> Excerpts from *Vortex Dynamics* by *P. G. Saffman*

Now know that, the roots of the dispersion relation locate at the interval $\beta a \in [j_{mn},j_{mn+1})$. The dispersion relation can be considered as a non-linear black box , for each $ka$, a $\beta a$ can be obtained. If $ka=0$, $\beta=j_{mn}$, otherwise, $\beta a \in (j_{mn},j_{mn+1})$. This helps to find the roots of the nonlinear equation.


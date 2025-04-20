### A Pluto.jl notebook ###
# v0.20.4

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
    #! format: on
end

# ╔═╡ dab5f085-ac90-4342-b2eb-60b097320de9
begin
	using CSV
	using DataFrames
	using PlutoUI
	using Plots
	using Statistics
	using Interact
	using Interpolations
	using Distributions
	using Random
	using StatsPlots 
	using LaTeXStrings
	using ColorSchemes
	using Optim
	using Loess, Polynomials, GLM, LsqFit

	plotly()  # Set the plotting backend separately, outside of `using`
end


# ╔═╡ f269ff32-3e80-4548-97c6-d337827db7aa
md"""
#### Parag Chettri ~ Santiago Berumen ~ Isaac Whitson
# An Analysis of Milky Way's Rotation Curve: A Julia Approach
## Data from Gaia DR 2
"""

# ╔═╡ 25699170-0fcd-11f0-2202-772042e9fd47
begin
    md"""
    ### QUERY

    ```sql
    SELECT TOP 100000 
        source_id, ra, dec, l, b, parallax, radial_velocity, 
        pmra, pmdec, parallax_error, radial_velocity_error, 
        pmra_error, pmdec_error
    FROM gaiadr2.gaia_source
    WHERE abs(b) < 2.5
        AND radial_velocity IS NOT NULL
        AND ra IS NOT NULL
        AND dec IS NOT NULL
        AND l IS NOT NULL
        AND b IS NOT NULL
        AND pmra IS NOT NULL
        AND pmdec IS NOT NULL
        AND source_id IS NOT NULL
        AND parallax > 0.01
        AND parallax < 99
        AND (
            (l > -5 AND l < 5) 
            OR 
            (l > 175 AND l < 185)
        )
    ```
    """
end


# ╔═╡ 80ec79c4-9971-47a2-a884-5417254a31fd
begin
	#hideall

	file_path = "DataSets/1744209227335O-result.csv"

	# Load the Gaia data
	dfInitial = CSV.read(file_path, DataFrame)

	# Randomly selecting 5000 rows
	df1 = dfInitial[randperm(nrow(dfInitial))[1:5000], :];
end;


# ╔═╡ ee84b436-12a5-4f87-8187-e064574632bb
begin
	md"""
	### Conversion from (RA, Dec) to Galactic (l, b)
	
	The relationship between Equatorial and Galactic coordinates is given by a rotation transformation. The Galactic north pole is at:
	
	
	$\alpha_{GP} = 192.8595^\circ, \quad \delta_{GP} = 27.1284^\circ$
	
	The Galactic center is at:
	
	$\alpha_{GC} = 266.4051^\circ, \quad \delta_{GC} = -28.9362^\circ$
	
	The inclination angle is:
	

	$i = 62.8717^\circ$
	
	Using these, we can compute the Galactic latitude $b$ for each star. The formula for $b$ is given by:
	
	$\sin(b) = \sin(\delta) \sin(\delta_{GP}) + \cos(\delta) \cos(\delta_{GP}) \cos(\alpha - \alpha_{GP})$
	
	
	Where: $\alpha$ is the Right Ascension $(RA)$ of the star in degrees
	
	Where: $\delta$ is the Declination $(Dec)$ of the star in degrees
	
	Where: $\alpha_{GP}$ and $\delta_{GP}$ are the coordinates of the Galactic north pole
	
	Where: $b$ is the Galactic latitude

	We can now restrict our data by only allowing a certain maximum deviation from 0 of $b$. 
	
	As deviation can be positive or negative, we will use the absolute value for $b$ to limit our values.
	
	We will choose a value of $b = 2.5^\circ$. 
	
	(This is reflected in our query as $abs(b) < 2.5$)
	"""
end

# ╔═╡ d1af481b-5a8e-469f-982a-401fbffd4590
begin
	md"""
	### Plotting the positions of our Stars (within our Galactic Coordinate Parameters)

	Our parameters were 
	"""
end
	

# ╔═╡ 72bfcd0b-47a4-4535-8a9d-2e059a6e86a1
begin
	md"""	
	The equation for the tangential velocity $V_{\text{tan}}$ is given by:
	
	
	$V_{\text{tan}} = 4.74 \times D \times \sqrt{\mu_l^2 + \mu_b^2}$
	
	
	Where $4.74$ is the conversion factor from milliarcseconds per year $(mas/yr)$ to ($km/s$)
	
	Where $D$ is the distance in kiloparsecs $(kpc)$
	
	Where $\mu_l$ and $\mu_b$ are the proper motions in the sky (in $mas/yr$).

	The array of $V_{\text{tan}}$ in $km/s$ is shown below.
	"""
end

# ╔═╡ 30c74cb8-8da9-4521-8f4a-48319390f5e6
begin
	md"""
	
	The true velocity $V_{\text{true}}$ is calculated using the radial velocity $V_{\text{radial}}$ and tangential velocity $V_{\text{tangential}}$:
	
	$V_{\text{true}} = \sqrt{V_{\text{radial}}^2 + V_{\text{tangential}}^2}$

	
	where: $V_{\text{true}}$ is the true velocity in $km/s$ 
	
	where: $V_{\text{radial}}$ is the radial velocity (along the line of sight) 


	
	where: $V_{\text{tangential}}$ is the tangential velocity (perpendicular to the line of sight)

	For simplicity's sake, we will ignore the fact that our sun and solar system are also orbiting the Galactic center.
	
	The array of $V_{\text{true}}$ in $km/s$ is shown below.
	
	"""
end

# ╔═╡ b07f485a-7168-4d87-acc2-026c1e2b4ad2
begin
	md"""
	Now there is something still wrong with our data. 
	
	Gaia, our source, does not correct for the orbital velocity of our sun and the solar system. This means, in our frame of reference, all the stars within our neighborhood with similar velocities seem to move with near-zero velocities. 

	We must now find a way to correct our data, particularly by shifting the rest frame to the Galactic Center and have our independant value of Distance be the Galactocentric Radius instead!
	
	"""
end

# ╔═╡ 34394d4f-5c13-4fa5-977e-1f0aee78dc2c
begin
	md"""
	# Translating our data to a Galactocentric Rest Frame
	"""
end

# ╔═╡ 8c1e2c97-381d-4c3b-af76-e269e46feb8c

begin
	md""" ## **Transformation Function**

The function `compute_galactocentric_kinematics` computes a star's Galactocentric kinematics using Galactic coordinates directly.

### **1. Galactic Coordinate Inputs**

Given:
- Galactic longitude ($l$) and latitude ($b$) in degrees
- Distance ($d$) in kiloparsecs (kpc)
- Proper motions in RA and Dec (pmRA, pmDec) in mas/yr
- Radial velocity $v_r$ in km/s

We compute the Cartesian position relative to the Galactic Center:

$x = d \cdot \cos(b) \cdot \cos(l)$  
$y = d \cdot \cos(b) \cdot \sin(l)$  
$z = d \cdot \sin(b)$

Then, we shift to Galactocentric coordinates by subtracting the Sun’s location (assuming $R_0 = 8.5 \, \text{kpc}$):

$X_{GC} = R_0 - x$  
$Y_{GC} = -y$  
$Z_{GC} = z$

The Galactocentric radius is:

$R_{GC} = \sqrt{X_{GC}^2 + Y_{GC}^2 + Z_{GC}^2}$

### **2. Velocity Components in Galactic Coordinates**

Proper motions are scaled using:

$\mu_l = \frac{\text{pmRA}}{\cos(b)}$  
$\mu_b = \text{pmDec}$

Velocities in the Galactic coordinate frame are computed using:

$v_l = 4.74 \cdot d \cdot \mu_l$  
$v_b = 4.74 \cdot d \cdot \mu_b$

Then, the space velocity components $U, V, W$ in the Galactic frame are:

$U = v_r \cdot \cos(b) \cdot \cos(l) - v_l \cdot \sin(l) - v_b \cdot \cos(l) \cdot \sin(b)$  
$V = v_r \cdot \cos(b) \cdot \sin(l) + v_l \cdot \cos(l) - v_b \cdot \sin(l) \cdot \sin(b)$  
$W = v_r \cdot \sin(b) + v_b \cdot \cos(b)$

### **3. Solar Motion and Galactic Rotation Correction**

We correct for the solar peculiar motion $(U_\odot, V_\odot, W_\odot) = (11.1, 12.24, 7.25)\,\text{km/s}$, and the Galactic circular speed at the Sun ($\Theta_0 = 220 \, \text{km/s}$):

$U_{\text{corr}} = U + U_\odot$  
$V_{\text{corr}} = V + V_\odot + \Theta_0$  
$W_{\text{corr}} = W + W_\odot$

### **4. Orbital (Azimuthal) Velocity**

The azimuthal velocity $V_\phi$, projected onto the Galactic plane, is given by:

$V_\phi = \frac{X_{GC} \cdot V_{\text{corr}} - Y_{GC} \cdot U_{\text{corr}}}{\sqrt{X_{GC}^2 + Y_{GC}^2}}$


Source: [(Stellar Kinematics)](https://en.wikipedia.org/wiki/Stellar_kinematics)

	
### **5. Applying the Transformation**

The function is applied to each row of the dataset `N` using the 'map' function from DataFrames.
	
	"""
end
	


# ╔═╡ 95432f3f-e4e5-48e4-abc2-522a1af4d402
begin
	md"""
	### The Keplerian Orbit
	"""
end

# ╔═╡ 7d85f4cd-e4e8-4076-afc5-33daf90c08fe
begin
	md"""

	As we know, Keplerian orbits are a fundamental part of Newtonian Physics, which we can overlay to check againt in this case. 

	$V_{\text{Keplerian}}(R) = \sqrt{\frac{GM_{\text{enclosed}}}{R}}$

	Where: $V_{\text{Keplerian}}$ is the orbital velocity in $km/s$

	Where: $M_{\text{enclosed}}$ is the mass within the galactic center that other galactic objects orbit.

	Where: $R$ is the distance from that mass.

	"""
end
	

# ╔═╡ 3ae662a1-2c66-4810-9770-514cc93897ff
begin
	md"""
	**Adjust X-Axis Scale (View Window)  
	Max R (kpc):** $(@bind xmax Slider(0:1:50; default=25, show_value=true))
	
	**Galactic Enclosed Mass / Solar Masses:** $(@bind mass_slider Slider(1e10:1e10:1e12, show_value=true, default=1e11))
		"""
end

# ╔═╡ 077c4e59-79e9-4108-85c8-cf6e06458394
begin
	md"""
	Needless to say, the rotation of the galaxy is clearly not bound just to Keplarian parameters.
	This has many reasons, but one of the most obvious ones being that there are many objects in this case and varying distances which increase the complexity with many interactions.
	"""
end

# ╔═╡ 63e51762-c84e-4839-9160-fb2dcf82e8ac
begin
	md"""
	**Adjust X-Axis Scale (View Window)  
	Max R (kpc):** $(@bind x_max5 Slider(1:1:50; default=25, show_value=true))
	"""
end

# ╔═╡ bab951b8-227f-4c89-a4a3-08f65a16fc35
begin
	md"""
	**Adjust Average Velocity:** $(@bind scaler1 Slider(100:5:300, show_value=true, default=200))
	"""
end

# ╔═╡ 627e2f9c-65d9-4645-a1a4-48c612876bd5
begin
	md"""
	**Adjust X-Axis Scale (View Window)  
	Max R (kpc):** $(@bind x_max4 Slider(1:1:50; default=25, show_value=true))
	"""
end

# ╔═╡ 711fdcf5-333a-44b9-a821-30cea343d89c
begin
    md"""# A Simple Model
	
Instead, we can understand this as a mathematicaL curve with minimal parameters. One of the possible approximations for this odd flat curve are as follows:
	
$$f(x) = V_{x \to \infty} \left( 1 - e^{-0.6 x} \right) - 0.15$$



    """
end

# ╔═╡ 57fa7cfd-8521-46ca-b2fd-70dc506e2593
begin
	md"""
	**Adjust Average Velocity:** $(@bind scaler Slider(180:5:270, show_value=true, default=200))
	"""
end

# ╔═╡ 9a619a6a-0814-4c21-a3be-50041c3e198f
begin
	md"""
	**Adjust X-Axis Scale (View Window)  
	Max R (kpc):** $(@bind x_max Slider(1:1:50; default=25, show_value=true))
	"""
end

# ╔═╡ 4d4b096d-a87a-482f-a19c-7918e7af47fa
begin
	md"""
	# A Complex Model - The Dark Matter-corrected Milky Way Rotation Curve
	
	The total rotation curve of the Milky Way is the quadrature sum of multiple contributions from the bulge, stellar disk, HI layer, H2 layer, and dark halo (a grand total of 9 parameters make up this model):
	
	
	$V_{\text{total}}(r) = \sqrt{V_{\text{bulge}}^2 + V_{\text{disk}}^2 + V_{\text{HI}}^2 + V_{\text{H2}}^2 + V_{\text{halo}}^2}$
	
	#### **1. Bulge (Blue Line)**
	The bulge follows a steep Keplerian-like rise and then quickly levels off:
	
	$V_{\text{bulge}}(r) = V_0 \frac{r}{(r^2 + a^2)^{3/4}}$
	
	where: $V_0$ is a scaling factor $(~260 km/s)$,
	where: $a$ is a scale radius $(~0.5–1 kpc)$.
	
	#### **2. Stellar Disk (Green Line)**
	The stellar disk follows an exponential profile:
	
	$V_{\text{disk}}(r) = V_0 \left(\frac{r}{R_d}\right) e^{-r/(2 R_d)}$
	
	where: $R_d$ is the disk scale length $(~3 kpc)$,
	where: $V_0$ is a normalization constant $(~150 km/s)$.
	
	#### **3. HI Layer (Yellow Line)**
	Neutral hydrogen contributes at intermediate radii:
	
	$V_{\text{HI}}(r) = V_0 \frac{r}{(r^2 + b^2)^{1/2}}$
	
	where: $b$ is a scale parameter $(~5 kpc)$,
	where: $V_0$ is a normalization constant $(~100 km/s)$.
	
	#### **4. H2 Layer (Pink Line)**
	Molecular hydrogen has a localized peak in the inner regions:
	
	$V_{\text{H2}}(r) = V_0 r e^{-r/R_{\text{H2}}}$
	
	where: $R_{\text{H2}}$ is a characteristic radius $(~4 kpc)$,
	where: $V_0$ is a normalization factor $(~80 km/s)$.
	
	#### **5. Dark Halo (Grey Line)**
	The dark matter halo is best described using an isothermal sphere model:
	
	$V_{\text{halo}}(r) = V_0 \frac{r}{\sqrt{r^2 + c^2}}$
	
	or for a nearly flat rotation curve:
	
	$V_{\text{halo}}(r) = V_{\infty} \left(1 - e^{-r/R_H}\right)$
	
	where: $V_{\infty}$ is the asymptotic velocity $(~205 km/s)$,
	where: $R_H$ is the halo scale radius $(~20 kpc)$,
	where: $c$ is a core radius $(~10 kpc)$.

	#### **4. Final Estimation (Black Line)**
	
	$V_{\text{total}}(r)$

	This is our final estimation as sum of the given distributions.

	Source: [(Ueshima, 2010)](https://www-sk.icrr.u-tokyo.ac.jp/xmass/publist/ueshima_PhD.pdf)
	"""
end


# ╔═╡ ab4b6e35-bbd6-4c5d-ae6c-cf3cb2092a85
begin
	md"""
	**Adjust X-Axis Scale (View Window)  
	Max R (kpc):** $(@bind xmax_2 Slider(1:1:50; default=25, show_value=true))
	"""
end

# ╔═╡ d7213a5c-f663-420d-99fd-1ea15758c4bc
begin
	md"""
	# F Tests and P-value

	tests

	"""
end

# ╔═╡ 4d5dbdaa-7859-45fe-ba01-5d783400543c


# ╔═╡ 93763c84-79f2-4cca-9a97-a871361c92a3
# ╠═╡ disabled = true
#=╠═╡
begin
	rss_simple = RSSSimple(N_New, r, V_total_curve)
	rss_complex = RSS(N_New, r, V_total_curve)
	
	k1 = 2   # number of parameters in simple model (example)
	k2 = 8   # number of parameters in complex model (example)
	n = length(N_New.Galactocentric_Radius)  # or after filtering
	
	F = f_test_rss(rss_simple, rss_complex, k1, k2, n)
	println("F-statistic = ", F)
end

  ╠═╡ =#

# ╔═╡ f7c435ef-deb8-49b5-8406-ec0b36116aac
# ╠═╡ disabled = true
#=╠═╡
function f_test_rss(RSS1, RSS2, k1, k2, n)
    numerator = (RSS1 - RSS2) / (k2 - k1)
    denominator = RSS2 / (n - k2)
    F_stat = numerator / denominator
    return F_stat
end

  ╠═╡ =#

# ╔═╡ bc5edd76-7aa4-494c-b78e-2949a032b0a9
begin
	md"""
	# Supporting Functions
	"""
end

# ╔═╡ f0807783-6915-49e5-a49f-eb770695b0ce
begin
	function plot_residuals(N_New, r, V_total_curve)
	    # Step 1: Interpolation with linear extrapolation
	    V_total_interp = interpolate((r,), V_total_curve, Gridded(Linear()))
	    V_total_interp = extrapolate(V_total_interp, Line())
	
	    # Step 2: Compute residuals
	    residuals = (N_New.V_phi) .- V_total_interp.(N_New.Galactocentric_Radius)
	
	    # Step 3: Sort for smooth plotting
	    sorted_indices = sortperm(N_New.Galactocentric_Radius)
	    sorted_radii = N_New.Galactocentric_Radius[sorted_indices]
	    sorted_residuals = residuals[sorted_indices]
	
	    # Step 4: Outlier removal (Option A: 0.25σ clipping)
	    residual_mean = mean(sorted_residuals)
	    residual_std = std(sorted_residuals)
	    threshold = 0.25 * residual_std
	    valid_indices = abs.(sorted_residuals .- residual_mean) .< threshold
	
	    ## Option B: radius threshold (uncomment to use)
	    # valid_indices = sorted_radii .< 1000
	
	    # Step 5: Filtering data
	    filtered_radii = sorted_radii[valid_indices]
	    filtered_residuals = sorted_residuals[valid_indices]
	
	    # Step 6: Ploting with or without error bars
	    if :V_phi_error in propertynames(N_New)
	        sorted_errors = N_New.V_phi_error[sorted_indices]
	        filtered_errors = sorted_errors[valid_indices]
	
	        plot(
	            filtered_radii, filtered_residuals, yerror=filtered_errors,
	            label="Residuals", marker=:circle, color=:blue, linewidth=2,
	            xlabel="Galactocentric Radius (kpc)", ylabel="Residual (km/s)",
	            title="Residuals: Observed - Model", alpha = 0.4
	        )
	    else
	        plot(
	            filtered_radii, filtered_residuals,
	            label="Residuals", marker=:circle, color=:blue, linewidth=2,
	            xlabel="Galactocentric Radius (kpc)", ylabel="Residual (km/s)",
	            title="Residuals: Observed - Model", alpha = 0.4
	        )
	    end
	
	    # Step 7: Adding horizontal line at zero
	    hline!([0], color=:black, linestyle=:dash, linewidth=2, label="Zero Residual")
	end
end


# ╔═╡ 7f1f1fc5-1077-4f56-94de-6a56521ad0e1
begin
	function RSS(N_New, r, V_total_curve)
	    # Step 1: Interpolation with linear extrapolation
	    V_total_interp = interpolate((r,), V_total_curve, Gridded(Linear()))
	    V_total_interp = extrapolate(V_total_interp, Line())
	
	    # Step 2: Compute residuals
	    residuals = (N_New.V_phi .+ 220) .- V_total_interp.(N_New.Galactocentric_Radius)
	
	    # Step 3: Sort for smooth plotting
	    sorted_indices = sortperm(N_New.Galactocentric_Radius)
	    sorted_radii = N_New.Galactocentric_Radius[sorted_indices]
	    sorted_residuals = residuals[sorted_indices]
	
	    # Step 4: Outlier removal (Option A: 0.25σ clipping)
	    residual_mean = mean(sorted_residuals)
	    residual_std = std(sorted_residuals)
	    threshold = 0.25 * residual_std
	    valid_indices = abs.(sorted_residuals .- residual_mean) .< threshold
	
	    ## Option B: radius threshold (uncomment to use)
	    # valid_indices = sorted_radii .< 1000
	
	    # Step 5: Filtering data
	    filtered_radii = sorted_radii[valid_indices]
	    filtered_residuals = sorted_residuals[valid_indices]
		
		# Step 6: Computing RSS
    RSSsimple = sum(filtered_residuals .^ 2)
	end
end


# ╔═╡ f5902359-7f0b-4244-8ebd-191358e3a227
function plot_residualsSimple(N_New, r, V_total_curve)
    # Step 1: Interpolation with linear extrapolation
    V_vals = V_total_curve.(r)
    V_total_interp = interpolate((r,), V_vals, Gridded(Linear()))
    V_total_interp = extrapolate(V_total_interp, Line())

    # Step 2: Computing residuals
    residuals = (N_New.V_phi) .- V_total_interp.(N_New.Galactocentric_Radius)

    # Step 3: Sorting for smooth plotting
    sorted_indices = sortperm(N_New.Galactocentric_Radius)
    sorted_radii = N_New.Galactocentric_Radius[sorted_indices]
    sorted_residuals = residuals[sorted_indices]

    # Step 4: Outlier removal (0.25σ clipping)
    residual_mean = mean(sorted_residuals)
    residual_std = std(sorted_residuals)
    threshold = 0.25 * residual_std
    valid_indices = abs.(sorted_residuals .- residual_mean) .< threshold

    # Step 5: Filtering data
    filtered_radii = sorted_radii[valid_indices]
    filtered_residuals = sorted_residuals[valid_indices]

    # Step 6: Plotting
    if :V_phi_error in propertynames(N_New)
        sorted_errors = N_New.V_phi_error[sorted_indices]
        filtered_errors = sorted_errors[valid_indices]

        plot(
            filtered_radii, filtered_residuals, yerror=filtered_errors,
            label="Residuals", marker=:circle, color=:blue, linewidth=2,
            xlabel="Galactocentric Radius (kpc)", ylabel="Residual (km/s)",
            title="Residuals: Observed - Model", xlims = (0,x_max)
        )
    else
        plot(
            filtered_radii, filtered_residuals,
            label="Residuals", marker=:circle, color=:blue, linewidth=2,
            xlabel="Galactocentric Radius (kpc)", ylabel="Residual (km/s)",
            title="Residuals: Observed - Model", xlims = (0,x_max)
        )
    end

    # Step 7: Adding horizontal line at zero
    hline!([0], color=:black, linestyle=:dash, linewidth=2, label="Zero Residual")

end


# ╔═╡ de2c68cb-640c-4fbd-83f8-eb8ef48104c0
function RSSSimple(N_New, r, V_total_curve)
    # Step 1: Interpolation with linear extrapolation
    V_vals = V_total_curve.(r)
    V_total_interp = interpolate((r,), V_vals, Gridded(Linear()))
    V_total_interp = extrapolate(V_total_interp, Line())

    # Step 2: Computing residuals
    residuals = (N_New.V_phi) .- V_total_interp.(N_New.Galactocentric_Radius)

    # Step 3: Sorting for smooth plotting
    sorted_indices = sortperm(N_New.Galactocentric_Radius)
    sorted_radii = N_New.Galactocentric_Radius[sorted_indices]
    sorted_residuals = residuals[sorted_indices]

    # Step 4: Outlier removal (0.25σ clipping)
    residual_mean = mean(sorted_residuals)
    residual_std = std(sorted_residuals)
    threshold = 0.25 * residual_std
    valid_indices = abs.(sorted_residuals .- residual_mean) .< threshold

    # Step 5: Filtering data
    filtered_radii = sorted_radii[valid_indices]
    filtered_residuals = sorted_residuals[valid_indices]

	# Step 6: Computing RSS
    RSSsimple = sum(filtered_residuals .^ 2)
end


# ╔═╡ efa3f064-bc0f-491b-b06a-d8cd6d42a174
function plot_residualsFlat(N_New, r, V_total_curve)
    # Step 1: Interpolation with linear extrapolation
    V_vals = V_total_curve.(r)
    V_total_interp = interpolate((r,), V_vals, Gridded(Linear()))
    V_total_interp = extrapolate(V_total_interp, Line())

    # Step 2: Computing residuals
    residuals = (N_New.V_phi) .- V_total_interp.(N_New.Galactocentric_Radius)

    # Step 3: Sorting for smooth plotting
    sorted_indices = sortperm(N_New.Galactocentric_Radius)
    sorted_radii = N_New.Galactocentric_Radius[sorted_indices]
    sorted_residuals = residuals[sorted_indices]

    # Step 4: Outlier removal (0.25σ clipping)
    residual_mean = mean(sorted_residuals)
    residual_std = std(sorted_residuals)
    threshold = 0.25 * residual_std
    valid_indices = abs.(sorted_residuals .- residual_mean) .< threshold

    # Step 5: Filtering data
    filtered_radii = sorted_radii[valid_indices]
    filtered_residuals = sorted_residuals[valid_indices]

    # Step 6: Plotting
    if :V_phi_error in propertynames(N_New)
        sorted_errors = N_New.V_phi_error[sorted_indices]
        filtered_errors = sorted_errors[valid_indices]

        plot(
            filtered_radii, filtered_residuals, yerror=filtered_errors,
            label="Residuals", marker=:circle, color=:blue, linewidth=2,
            xlabel="Galactocentric Radius (kpc)", ylabel="Residual (km/s)",
            title="Residuals: Observed - Model", xlims = (0,x_max4)
        )
    else
        plot(
            filtered_radii, filtered_residuals,
            label="Residuals", marker=:circle, color=:blue, linewidth=2,
            xlabel="Galactocentric Radius (kpc)", ylabel="Residual (km/s)",
            title="Residuals: Observed - Model", xlims = (0,x_max4)
        )
    end

    # Step 7: Adding horizontal line at zero
    hline!([0], color=:black, linestyle=:dash, linewidth=2, label="Zero Residual")

end


# ╔═╡ 54ce63de-9640-4c00-b685-f8163f5f4bd8
begin
		
		function V_disk(r, V0, R_d)
		    return V0 * (r ./ R_d) .* exp.(-r ./ (2 * R_d))
		end
		
		function V_HI(r, V0, b)
		    return V0 * (r ./ (r.^2 .+ b^2).^(1/2))
		end
		
		function V_H2(r, V0, R_H2)
		    return V0 * r .* exp.(-r ./ R_H2)
		end
		
		function V_halo(r, V0, c)
		    return V0 * (r ./ sqrt.(r.^2 .+ c^2))
		end
	
	 V0_bulge = 260
		a_bulge = 1
		
		V0_disk = 150
		R_d = 3
		
		V0_HI = 100
		b_HI = 5
		
		V0_H2 = 80
		R_H2 = 4
		
		V_infinity_halo = 205
		R_H = 20
		c_halo = 10
		function V_bulge(r, V0, a)
		    return V0 * (r ./ (r.^2 .+ a^2).^(3/4))
		end
end

# ╔═╡ cb94e83a-2e31-45d9-b4e0-b352d4595e99
function filter_galactic_plane(df1; threshold=2.5)
	    # Convert RA and Dec to radians
	    ra_rad = deg2rad.(df1.ra)
	    dec_rad = deg2rad.(df1.dec)
	
	    # Galactic north pole coordinates in radians
	    ra_gp = deg2rad(192.8595)
	    dec_gp = deg2rad(27.1284)
	
	    # Compute sin(b)
	    sin_b = sin.(dec_rad) .* sin(dec_gp) .+ cos.(dec_rad) .* cos(dec_gp) .* cos.(ra_rad .- ra_gp)
	
	    # Compute Galactic latitude b in degrees
	    b = rad2deg.(asin.(sin_b))
	
	    # Filter stars near the Galactic plane (|b| < threshold)
	    return df1[abs.(b) .< threshold, :]
	end

# ╔═╡ 2a36ce70-5dc5-401d-864e-9b5e479b7913
begin

	
	# Apply filtering: only stars with non-missing parallax, positive parallax, and proper motions
	filtered_df = filter_galactic_plane(df1)
	
	N = filter(row -> !ismissing(row.parallax) && row.parallax > 0 && !ismissing(row.pmra) && !ismissing(row.pmdec), filtered_df)
end;


# ╔═╡ cc962df0-b6a0-481c-a95b-9265930d38d7
begin
	# Further filtering for stars with non-missing parallax, pmra, and pmdec

		N[!, :distance_kpc_error] = (1.0 ./ N.parallax_error .+ N.parallax) - (1.0 ./ N.parallax .- N.parallax_error )
		N[!, :distance_kpc] = 1.0 ./ N.parallax
		N[!, :proper_motion] = sqrt.(N.pmra.^2 .+ N.pmdec.^2)
		N[!, :proper_motion_error] = sqrt.(N.pmra_error.^2 .+ N.pmdec_error.^2)
		N[!, :tangential_velocity_error] = 4.74 .* N.proper_motion_error .* 		N.distance_kpc_error
		N[!, :tangential_velocity] = 4.74 .* N.proper_motion .* N.distance_kpc
end

# ╔═╡ 6f370423-0a7e-4a20-9c69-81fe352020a0
begin
    # Compute true velocity (km/s), filtering out missing radial velocities
    N[!, :true_velocity_error] = sqrt.(N.tangential_velocity_error.^2 .+ coalesce.(N.radial_velocity_error, 0.0).^2)
    N[!, :tangential_velocity_corrected] = N.tangential_velocity
    N[!, :true_velocity] = sqrt.(N.tangential_velocity_corrected.^2 .+ coalesce.(N.radial_velocity, 0.0).^2)
end

# ╔═╡ 0ccfae55-6162-4c97-a04d-2dd3c5d9370e
begin
	function compute_galactocentric_kinematics(ra, dec, l, b, d_kpc, pmra_masyr, pmdec_masyr, v_r_kms)
	    # Constants
	    distance_sun_to_galactic_center = 8.5  # kpc (Distance from Galactic Center to Sun)
	    Θ_0 = 220.0  # km/s (Galactic rotation at Sun)
	    U_sun, V_sun, W_sun = 11.1, 12.24, 7.25  # km/s
	
	    # 1. Convert star position to Galactocentric Cartesian coordinates
	    X = d_kpc * cosd(b) * cosd(l)
	    Y = d_kpc * cosd(b) * sind(l)
	    Z = d_kpc * sind(b)
	
	    X_GC = distance_sun_to_galactic_center - X
	    Y_GC = -Y  # Sign flip for Galactic rotation direction
	    Z_GC = Z
	
	    R_GC = sqrt(X_GC^2 + Y_GC^2 + Z_GC^2)  # Galactocentric radius
	
	    # 2. Compute velocities in Galactic coordinates
	    μ_l = pmra_masyr / cosd(b)  # Corrected proper motion in l
	    μ_b = pmdec_masyr

		
		# v_tangential -  km/s= 4.74(conversion factor) * distance - kpc * proper motion in milliarcseconds/year
		
	    v_l = 4.74 * d_kpc * μ_l  # km/s
	    v_b = 4.74 * d_kpc * μ_b  # km/s
	
	    # 3D velocity components (U, V, W)
	    U = v_r_kms * cosd(b) * cosd(l) - v_l * sind(l) - v_b * cosd(l) * sind(b)
	    V = v_r_kms * cosd(b) * sind(l) + v_l * cosd(l) - v_b * sind(l) * sind(b)
	    W = v_r_kms * sind(b) + v_b * cosd(b)
	
	    # Correct for Sun's motion and Galactic rotation
	    U_corr = U + U_sun
	    V_corr = V + V_sun + Θ_0
	    W_corr = W + W_sun
	
	    # 3. Compute orbital (azimuthal) velocity
	    V_phi = (X_GC * V_corr - Y_GC * U_corr) / sqrt(X_GC^2 + Y_GC^2)
	
	    return R_GC, U_corr, V_corr, W_corr, V_phi
	end
	
	# Apply to your dataframe
	results = map(row -> compute_galactocentric_kinematics(row.ra, row.dec, row.l, row.b, row.distance_kpc, row.pmra, row.pmdec, row.radial_velocity), eachrow(N))
	
	# Add columns to dataframe
	N[!, :Galactocentric_Radius] = [r[1] for r in results]
	N[!, :V_X] = [r[2] for r in results]
	N[!, :V_Y] = [r[3] for r in results]
	N[!, :V_Z] = [r[4] for r in results]
	N[!, :V_phi] = [r[5] for r in results]
	
	# Filter (if needed)
	N_New = filter(row -> row.Galactocentric_Radius > 0 && row.V_phi > 0,  N)
end

# ╔═╡ efda5197-b8ff-467a-8d32-6d9b7992a7ec
begin

    ra_rad = deg2rad.(N_New.ra)
    dec_rad = deg2rad.(N_New.dec)
    
    x = N_New.distance_kpc .* cos.(dec_rad) .* cos.(ra_rad)
    y = N_New.distance_kpc .* cos.(dec_rad) .* sin.(ra_rad)
    z = N_New.distance_kpc .* sin.(dec_rad)

    # Creating scatter plot
    p = scatter(x, y, z, markersize = 2, title = "A Map of our DataSet (Over the Disc of the Milky Way", xlabel = "kpc_x", ylabel = "kpc_y", zlabel = "kpc_z", legend=false)
    
    # Setting manual limits
    xlims!(p, -20, 20)
    ylims!(p, -20, 20)
    zlims!(p, -20, 20)
    
    # Create a galactic disc (radius 13 kpc)
    θ = range(0, 2π, length=100)  # Angle for circle
    rr = range(0, 13, length=20)   # Radius steps
    
    # Creating disc points (initially in XY plane)
    disc_x = [ri * cos(θi) for ri in rr, θi in θ]
    disc_y = [ri * sin(θi) for ri in rr, θi in θ] .- 7.7
    disc_z = zeros(size(disc_x))
    
    # Known galactic plane tilt (approximate)
    tilt_angle = deg2rad(28)  
    
    # Rotating disc to match galactic plane tilt
    disc_y_rot = disc_y .* cos(tilt_angle) .- disc_z .* sin(tilt_angle)
    disc_z_rot = disc_y .* sin(tilt_angle) .+ disc_z .* cos(tilt_angle)
    
    # Plotting the disc
    surface!(p, disc_x, disc_y_rot, disc_z_rot, alpha=0.5, color=:blue, legend = false )
    
    p  # Displaying the plot
end

# ╔═╡ fb84175b-122e-4121-8498-085810627024
begin

	scatter(N_New.Galactocentric_Radius, N_New.V_phi, marker=:circle, xlabel="Galactocentric Radius (kpc)", ylabel="Orbital Velocity (km/s)", title="Orbital Velocity and the Keplerian", legend=true, label="Observed Object", alpha = 0.4, xlims=(0, xmax), ylims=(0, 500))


	G = 4.302e-6  # kpc * (km/s)^2 / Msun
	
	keplerian_velocity(r, M) = sqrt(G * M / r)
	
	
	### Galactic rotation curve
	r_values = 0.1:0.1:xmax # Range of distances (kpc)
	v_kepler = keplerian_velocity.(r_values, mass_slider)  # Keplerian velocity based on current mass
	
	### Generating scatter plot and overlaying Keplerian curve
	
	plot!(r_values, v_kepler, label="Keplerian Curve", legend = true, linewidth=2, color=:red)
end

# ╔═╡ c93a3123-22d0-4433-a6d3-34cc1395016a
begin
	
	# 1. Prepare the data
	# Replace these with your actual data vectors
	r_data = N_New.Galactocentric_Radius  # Example Galactocentric Radius data
	v_data = N_New.V_phi  # Example Orbital Velocity data
	
	# Remove any invalid points
	valid = .!(isnan.(r_data) .| isnan.(v_data))
	r_data = r_data[valid]
	v_data = v_data[valid]
	
	# Create a DataFrame
	df = DataFrame(Radius = r_data, Velocity = v_data)
	
	# 2. Fit a linear regression model
	modelLR = lm(@formula(Velocity ~ Radius), df)
	
	# 3. Generate predictions
	r_plot = range(minimum(r_data), maximum(r_data), length=200)
	df_plot = DataFrame(Radius = r_plot)
	predictions = predict(modelLR, df_plot)
	
	# 4. Create plot
	scatter(r_data, v_data, label="Data", color=:blue, alpha=0.6,
	        xlabel="Galactocentric Radius (kpc)", ylabel="Orbital Velocity (km/s)",
	        title="Linear Regression Fit", legend=:topleft)
	
	plot!(r_plot, predictions, label="Best-Fit Line", color=:red, lw=2)
	
end

# ╔═╡ 1f1cab74-5f8d-4ab4-a12e-5e53a89d9080
begin
	# 5. Print model summary
	println("Linear Regression Model Summary:")
	println(coeftable(modelLR))
end

# ╔═╡ d62f9aab-db3d-49f3-9844-ed9b091903ef
begin
	
	
	# Scatter plot with filtered data (Assuming `N_filtered` is already defined)
	scatter(N_New.Galactocentric_Radius, N_New.V_phi, marker=:circle, label="Data Points", xlabel="Galactocentric Radius (kpc)", ylabel="Orbital Velocity (km/s)", title="Orbital Velocity vs. Galactocentric Radius", legend=false, alpha=0.3 , xlims=(0, x_max4), ylims=(0, 500))
	

	# Customize the plot further
	xlabel!("Galactocentric Radius (kpc)")
	ylabel!("Orbital Velocity (km/s)")
	title!("(BigManBlastoise)")
	xlims!(0, x_max4)
	ylims!(0, 450)

	model1(x) = scaler1
	
	# Overlay the model curve
	plot!(model1, 0:0.1:25, color=:red, lw=2, label="Model", xlims= (0,x_max4))

end


# ╔═╡ bb8551a1-8e35-49a4-8aa8-1fc89167ae7e
begin
	
	
	# Scatter plot with filtered data (Assuming `N_filtered` is already defined)
	scatter(N_New.Galactocentric_Radius, N_New.V_phi, marker=:circle, label="Data Points", xlabel="Galactocentric Radius (kpc)", ylabel="Orbital Velocity (km/s)", title="Orbital Velocity vs. Galactocentric Radius", legend=false, alpha=0.3 , xlims=(0, x_max), ylims=(0, 500))
	

	# Customize the plot further
	xlabel!("Galactocentric Radius (kpc)")
	ylabel!("Orbital Velocity (km/s)")
	title!("Simple Model")
	xlims!(0, x_max)
	ylims!(0, 450)

	model(x) = scaler * (1 - exp(-0.6 * x)) - 0.15
	
	# Overlay the model curve
	plot!(model, 0:0.1:25, color=:red, lw=2, label="Model", xlims= (0,x_max))
end


# ╔═╡ ddba2999-edab-4899-86d7-44c822e33ac1
begin
	
	# Defining the range of r (distance in kpc)
	r = 0:0.1:xmax_2  # 0 to user selected value for kpc
	
	
	# Calculating each component
	V_bulge_curve = V_bulge(r, V0_bulge, a_bulge)
	V_disk_curve = V_disk(r, V0_disk, R_d)
	V_HI_curve = V_HI(r, V0_HI, b_HI)
	V_H2_curve = V_H2(r, V0_H2, R_H2)
	V_halo_curve = V_halo(r, V_infinity_halo, c_halo)
	
	# Total rotation curve (quadrature sum)
	V_total_curve = sqrt.(V_bulge_curve.^2 + V_disk_curve.^2 + V_HI_curve.^2 + V_H2_curve.^2 + V_halo_curve.^2)
	
	# Scatter plot with filtered data (Assuming `N_filtered` is already defined)
	scatter(N_New.Galactocentric_Radius, N_New.V_phi, marker=:circle, label="Data Points", xlabel="Galactocentric Radius (kpc)", ylabel="Orbital Velocity (km/s)", title="Orbital Velocity vs. Galactocentric Radius", legend=false, alpha=0.1 , xlims=(0, xmax_2), ylims=(0, 500))
	
	# Overlay the rotation curve components and the total curve
	plot!(r, V_bulge_curve, label="Bulge", color=:blue, linewidth=2)
	plot!(r, V_disk_curve, label="Stellar Disk", color=:green, linewidth=2)
	plot!(r, V_HI_curve, label="HI Layer", color=:yellow, linewidth=2)
	plot!(r, V_H2_curve, label="H2 Layer", color=:pink, linewidth=2)
	plot!(r, V_halo_curve, label="Dark Halo", color=:gray, linewidth=2)
	plot!(r, V_total_curve, label="Total", color=:black, linewidth=3)
	
	# Customize the plot further
	xlabel!("Galactocentric Radius (kpc)")
	ylabel!("Orbital Velocity (km/s)")
	title!("Orbital Velocity vs. Galactocentric Radius")
	xlims!(0, xmax_2)
	ylims!(0, 450)
end


# ╔═╡ 6a9d29ad-ecb7-4f9a-a86f-7e67e658a7de
plot_residualsFlat(N_New, r, model1)

# ╔═╡ 1ad1dc9b-f12d-4e9a-8ea8-c3c85ab738c3
plot_residualsSimple(N_New, r, model)

# ╔═╡ 5e980df1-7333-4210-97c1-7aa3f0de9a67
plot_residuals(N_New, r, V_total_curve)

# ╔═╡ 50dde6f7-8c8b-439d-9dbe-498ef41aa5b3
begin

    # Real data from Gaia analysis
    r_data1 = N_New.Galactocentric_Radius
    V_obs1 = N_New.V_phi

    # Define component models

 # Initial parameter guess (can be tuned)
    initial_guess = [
        200, 1.5,  # bulge
        120, 3.0,  # disk
        90, 4.0,   # HI
        70, 3.5,   # H2
        190, 15.0  # halo
    ]
	
	function rotation_model(p, r)
	    V0_bulge, a_bulge, V0_disk, R_d, V0_HI, b_HI, V0_H2, R_H2, V_halo_inf, c_halo = p
	    sqrt.(
	        V_bulge(r, V0_bulge, a_bulge).^2 .+  
	        V_disk(r, V0_disk, R_d).^2 .+
	        V_HI(r, V0_HI, b_HI).^2 .+
	        V_H2(r, V0_H2, R_H2).^2 .+
	        V_halo(r, V_halo_inf, c_halo).^2
	    )
	end
	
    # Fit the model to real Gaia-based data
    fit = curve_fit(rotation_model, r_data1, V_obs1, initial_guess)
    fitted_params = fit.param

    println("\n🔍 Recovered Parameters from Gaia Data:")
    param_labels = [
        "V0_bulge", "a_bulge", "V0_disk", "R_d", "V0_HI", "b_HI",
        "V0_H2", "R_H2", "V_infinity_halo", "c_halo"
    ]
    for (label, val) in zip(param_labels, fitted_params)
        println("$label: ", round(val, digits=3))
    end

    # Plot: Observed vs Fitted
    plot(r_data1, V_obs1, label="Observed V_phi", seriestype=:scatter, alpha=0.3, color=:gray)
    plot!(r_data1, rotation_model(fitted_params, r_data1), label="Fitted Curve", lw=3, color=:blue)
    xlabel!("Galactocentric Radius (kpc)")
    ylabel!("Orbital Velocity (km/s)")
    title!("Nonlinear Fit to Gaia-based Orbital Velocity Data")
end


# ╔═╡ 868ae38e-5355-42cf-b645-1ae6d6b443c8
begin
	function f_test_rss(RSS1, RSS2, k1, k2, n)
	    numerator = (RSS1 - RSS2) / (k2 - k1)
	    denominator = RSS2 / (n - k2)
	    F_stat = numerator / denominator
	    return F_stat
	end
	
	RSS1=RSSSimple(N_New, r, model)
	RSS2=RSS(N_New, r, V_total_curve)
	k1 = 3
	k2 = 9
	n = 5000
	F = f_test_rss(RSS1, RSS2,k1, k2, n)
	println("F stat value of, ", F)
end
		

# ╔═╡ 6e0c44d6-af43-4dbd-ab2f-d36053e2559a
begin
	w = k2 - k1
	w1 = n - k2
	p_value = 1 - cdf(FDist(w, w1), F)
	println("p-value = ", p_value)
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
ColorSchemes = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
GLM = "38e38edf-8417-5370-95a0-9cbb8c7f171a"
Interact = "c601a237-2ae4-5e1e-952c-7a85b0c7eef1"
Interpolations = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
Loess = "4345ca2d-374a-55d4-8d30-97f9976e7612"
LsqFit = "2fda8390-95c7-5789-9bda-21331edee243"
Optim = "429524aa-4258-5aef-a3af-852621145aeb"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Polynomials = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
StatsPlots = "f3b207a7-027a-5e70-b257-86293d7955fd"

[compat]
CSV = "~0.10.15"
ColorSchemes = "~3.29.0"
DataFrames = "~1.7.0"
Distributions = "~0.25.117"
GLM = "~1.9.0"
Interact = "~0.10.5"
Interpolations = "~0.15.1"
LaTeXStrings = "~1.4.0"
Loess = "~0.6.4"
LsqFit = "~0.15.1"
Optim = "~1.10.0"
Plots = "~1.40.9"
PlutoUI = "~0.7.61"
Polynomials = "~4.0.19"
Statistics = "~1.11.1"
StatsPlots = "~0.15.7"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.2"
manifest_format = "2.0"
project_hash = "40e85fdc456b948f0126c19b5a007899e15b9d8d"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "cd8b948862abee8f3d3e9b73a102a9ca924debb0"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.2.0"
weakdeps = ["SparseArrays", "StaticArrays"]

    [deps.Adapt.extensions]
    AdaptSparseArraysExt = "SparseArrays"
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Arpack]]
deps = ["Arpack_jll", "Libdl", "LinearAlgebra", "Logging"]
git-tree-sha1 = "9b9b347613394885fd1c8c7729bfc60528faa436"
uuid = "7d9fca2a-8960-54d3-9f78-7d1dccf2cb97"
version = "0.5.4"

[[deps.Arpack_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "OpenBLAS_jll", "Pkg"]
git-tree-sha1 = "5ba6c757e8feccf03a1554dfaf3e26b3cfc7fd5e"
uuid = "68821587-b530-5797-8361-c406ea357684"
version = "3.5.1+1"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra"]
git-tree-sha1 = "017fcb757f8e921fb44ee063a7aafe5f89b86dd1"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.18.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceCUDSSExt = "CUDSS"
    ArrayInterfaceChainRulesCoreExt = "ChainRulesCore"
    ArrayInterfaceChainRulesExt = "ChainRules"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceReverseDiffExt = "ReverseDiff"
    ArrayInterfaceSparseArraysExt = "SparseArrays"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    CUDSS = "45b445bb-4962-46a0-9369-b4df9d0f772e"
    ChainRules = "082447d4-558c-5d27-93f4-14fc19e9eca2"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    ReverseDiff = "37e2e3b7-166d-5795-8a7a-e32c996b4267"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.AssetRegistry]]
deps = ["Distributed", "JSON", "Pidfile", "SHA", "Test"]
git-tree-sha1 = "b25e88db7944f98789130d7b503276bc34bc098e"
uuid = "bf4720bc-e11a-5d0c-854e-bdca1663c893"
version = "0.1.0"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1b96ea4a01afe0ea4090c5c8039690672dd13f2e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.9+0"

[[deps.CSSUtil]]
deps = ["Colors", "JSON", "Markdown", "Measures", "WebIO"]
git-tree-sha1 = "b9fb4b464ec10e860abe251b91d4d049934f7399"
uuid = "70588ee8-6100-5070-97c1-3cb50ed05fe8"
version = "0.1.1"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "PrecompileTools", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings", "WorkerUtilities"]
git-tree-sha1 = "deddd8725e5e1cc49ee205a1964256043720a6c3"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.15"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "2ac646d71d0d24b44f3f8c84da8c9f4d70fb67df"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.4+0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "1713c74e00545bfe14605d2a2be1712de8fbcb58"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.25.1"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.Clustering]]
deps = ["Distances", "LinearAlgebra", "NearestNeighbors", "Printf", "Random", "SparseArrays", "Statistics", "StatsBase"]
git-tree-sha1 = "3e22db924e2945282e70c33b75d4dde8bfa44c94"
uuid = "aaaa29a8-35af-508c-8bc3-b662a17a0fe5"
version = "0.15.8"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "962834c22b66e32aa10f7611c08c8ca4e20749a9"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.8"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "403f2d8e209681fcbd9468a8514efff3ea08452e"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.29.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "362a287c3aa50601b0bc359053d5c2468f0e7ce0"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.11"

[[deps.CommonSubexpressions]]
deps = ["MacroTools"]
git-tree-sha1 = "cda2cfaebb4be89c9084adaca7dd7333369715c5"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.1"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "8ae8d32e09f0dcf42a36b90d4e17f5dd2e4c4215"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.16.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "d9d26935a0bcffc87d2613ce14c527c99fc543fd"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.5.0"

[[deps.ConstructionBase]]
git-tree-sha1 = "76219f1ed5771adbb096743bff43fb5fdd4c1157"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.8"

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseLinearAlgebraExt = "LinearAlgebra"
    ConstructionBaseStaticArraysExt = "StaticArrays"

    [deps.ConstructionBase.weakdeps]
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "DataStructures", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrecompileTools", "PrettyTables", "Printf", "Random", "Reexport", "SentinelArrays", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "fb61b4812c49343d7ef0b533ba982c46021938a6"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.7.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "4e1fe97fdaed23e9dc21d4d664bea76b65fc50a0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.22"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"
version = "1.11.0"

[[deps.Dbus_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fc173b380865f70627d7dd1190dc2fce6cc105af"
uuid = "ee1fde0b-3d02-5ea6-8484-8dfef6360eab"
version = "1.14.10+0"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "c7e3a542b999843086e2f29dac96a618c105be1d"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.12"
weakdeps = ["ChainRulesCore", "SparseArrays"]

    [deps.Distances.extensions]
    DistancesChainRulesCoreExt = "ChainRulesCore"
    DistancesSparseArraysExt = "SparseArrays"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "03aa5d44647eaec98e1920635cdfed5d5560a8b9"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.117"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
git-tree-sha1 = "e7b7e6f178525d17c720ab9c081e4ef04429f860"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.4"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a4be429317c42cfae6a7fc03c31bad1970c310d"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+1"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "d36f682e590a83d63d1c7dbd287573764682d12a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.11"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d55dffd9ae73ff72f1c0482454dcf2ec6c6c4a63"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.6.5+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "53ebe7511fa11d33bec688a9178fac4e49eeee00"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.2"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "466d45dc38e15794ec7d5d63ec03d776a9aff36e"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.4+1"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "7de7c78d681078f027389e067864a8d53bd7c3c9"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.8.1"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4d81ed14783ec49ce9f2e168208a12ce1815aa25"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+3"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates"]
git-tree-sha1 = "3bab2c5aa25e7840a4b065805c0cdfc01f3068d2"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.24"
weakdeps = ["Mmap", "Test"]

    [deps.FilePathsBase.extensions]
    FilePathsBaseMmapExt = "Mmap"
    FilePathsBaseTestExt = "Test"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "6a70198746448456524cb442b8af316927ff3e1a"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.13.0"
weakdeps = ["PDMats", "SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Setfield"]
git-tree-sha1 = "f089ab1f834470c525562030c8cfde4025d5e915"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.27.0"

    [deps.FiniteDiff.extensions]
    FiniteDiffBandedMatricesExt = "BandedMatrices"
    FiniteDiffBlockBandedMatricesExt = "BlockBandedMatrices"
    FiniteDiffSparseArraysExt = "SparseArrays"
    FiniteDiffStaticArraysExt = "StaticArrays"

    [deps.FiniteDiff.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "301b5d5d731a0654825f1f2e906990f7141a106b"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.16.0+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "a2df1b776752e3f344e5116c06d75a10436ab853"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.38"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "2c5512e11c791d1baed2049c5652441b28fc6a31"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "846f7026a9decf3679419122b49f8a1fdb48d2d5"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.16+0"

[[deps.FunctionalCollections]]
deps = ["Test"]
git-tree-sha1 = "04cb9cfaa6ba5311973994fe3496ddec19b6292a"
uuid = "de31a74c-ac4f-5751-b3fd-e18cd04993ca"
version = "0.5.0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"
version = "1.11.0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "fcb0584ff34e25155876418979d4c8971243bb89"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.0+2"

[[deps.GLM]]
deps = ["Distributions", "LinearAlgebra", "Printf", "Reexport", "SparseArrays", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns", "StatsModels"]
git-tree-sha1 = "273bd1cd30768a2fddfa3fd63bbc746ed7249e5f"
uuid = "38e38edf-8417-5370-95a0-9cbb8c7f171a"
version = "1.9.0"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Qt6Wayland_jll", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "0ff136326605f8e06e9bcf085a356ab312eef18a"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.13"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "9cb62849057df859575fc1dda1e91b82f8609709"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.13+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "b0036b392358c80d2d2124746c2bf3d48d457938"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.82.4+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a6dbda1fd736d60cc477d99f2e7a042acfa46e8"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.15+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "PrecompileTools", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "c67b33b085f6e2faf8bf79a61962e7339a81129c"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.15"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll"]
git-tree-sha1 = "55c53be97790242c29031e5cd45e8ac296dadda3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "8.5.0+0"

[[deps.HypergeometricFunctions]]
deps = ["LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "2bd56245074fab4015b9174f24ceba8293209053"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.27"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.InlineStrings]]
git-tree-sha1 = "6a9fde685a7ac1eb3495f8e812c5a7c3711c2d5e"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.3"

    [deps.InlineStrings.extensions]
    ArrowTypesExt = "ArrowTypes"
    ParsersExt = "Parsers"

    [deps.InlineStrings.weakdeps]
    ArrowTypes = "31f734f8-188a-4ce0-8406-c8a06bd891cd"
    Parsers = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "0f14a5456bdc6b9731a5682f439a672750a09e48"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2025.0.4+0"

[[deps.Interact]]
deps = ["CSSUtil", "InteractBase", "JSON", "Knockout", "Observables", "OrderedCollections", "Reexport", "WebIO", "Widgets"]
git-tree-sha1 = "c5091992248c7134af7c90554305c600d5d9012b"
uuid = "c601a237-2ae4-5e1e-952c-7a85b0c7eef1"
version = "0.10.5"

[[deps.InteractBase]]
deps = ["Base64", "CSSUtil", "Colors", "Dates", "JSExpr", "JSON", "Knockout", "Observables", "OrderedCollections", "Random", "WebIO", "Widgets"]
git-tree-sha1 = "aa5daeff326db0a9126a225b58ca04ae12f57259"
uuid = "d3863d7c-f0c8-5437-a7b4-3ae773c01009"
version = "0.10.10"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"
version = "1.11.0"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "88a101217d7cb38a7b481ccd50d21876e1d1b0e0"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.15.1"
weakdeps = ["Unitful"]

    [deps.Interpolations.extensions]
    InterpolationsUnitfulExt = "Unitful"

[[deps.InvertedIndices]]
git-tree-sha1 = "6da3c4316095de0f5ee2ebd875df8721e7e0bdbe"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.1"

[[deps.IrrationalConstants]]
git-tree-sha1 = "e2222959fbc6c19554dc15174c81bf7bf3aa691c"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.4"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLFzf]]
deps = ["REPL", "Random", "fzf_jll"]
git-tree-sha1 = "1d4015b1eb6dc3be7e6c400fbd8042fe825a6bac"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.10"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "a007feb38b422fbdab534406aeca1b86823cb4d6"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.7.0"

[[deps.JSExpr]]
deps = ["JSON", "MacroTools", "Observables", "WebIO"]
git-tree-sha1 = "b413a73785b98474d8af24fd4c8a975e31df3658"
uuid = "97c1335a-c9c5-57fe-bc5d-ec35cebe8660"
version = "0.5.4"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "eac1206917768cb54957c65a615460d87b455fc1"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.1.1+0"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "7d703202e65efa1369de1279c162b915e245eed1"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.9"

[[deps.Knockout]]
deps = ["JSExpr", "JSON", "Observables", "Test", "WebIO"]
git-tree-sha1 = "91835de56d816864f1c38fb5e3fad6eb1e741271"
uuid = "bcebb21b-c2e3-54f8-a781-646b90f6d2cc"
version = "0.2.6"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "170b660facf5df5de098d866564877e119141cbd"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.2+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "aaafe88dccbd957a8d82f7d05be9b69172e0cee3"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "4.0.1+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "78211fb6cbc872f77cad3fc0b6cf647d923f4929"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "18.1.7+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c602b1127f4751facb671441ca72715cc95938a"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.3+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "dda21b8cbd6a6c40d9d02a73230f9d70fed6918c"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.4.0"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "cd714447457c660382fe634710fb56eb255ee42e"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.6"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SparseArraysExt = "SparseArrays"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"
version = "1.11.0"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "27ecae93dd25ee0909666e6835051dd684cc035e"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+2"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll"]
git-tree-sha1 = "d77592fa54ad343c5043b6f38a03f1a3c3959ffe"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.11.1+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "ff3b4b9d35de638936a525ecd36e86a8bb919d11"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.7.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "df37206100d39f79b3376afb6b9cee4970041c61"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.51.1+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "be484f5c92fad0bd8acfef35fe017900b0b73809"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.18.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a31572773ac1b745e0343fe5e2c8ddda7a37e997"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.41.0+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "4ab7581296671007fc33f07a721631b8855f4b1d"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.7.1+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "321ccef73a96ba828cd51f2ab5b9f917fa73945a"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.41.0+0"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "e4c3be53733db1051cc15ecf573b1042b3a712a1"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.3.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.Loess]]
deps = ["Distances", "LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "f749e7351f120b3566e5923fefdf8e52ba5ec7f9"
uuid = "4345ca2d-374a-55d4-8d30-97f9976e7612"
version = "0.6.4"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "13ca9e2586b89836fd20cccf56e57e2b9ae7f38f"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.29"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "f02b56007b064fbfddb4c9cd60161b6dd0f40df3"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.1.0"

[[deps.LsqFit]]
deps = ["Distributions", "ForwardDiff", "LinearAlgebra", "NLSolversBase", "Printf", "StatsAPI"]
git-tree-sha1 = "f386224fa41af0c27f45e2f9a8f323e538143b43"
uuid = "2fda8390-95c7-5789-9bda-21331edee243"
version = "0.15.1"

[[deps.MIMEs]]
git-tree-sha1 = "c64d943587f7187e751162b3b84445bbbd79f691"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "1.1.0"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "5de60bc6cb3899cd318d80d627560fae2e2d99ae"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2025.0.1+1"

[[deps.MacroTools]]
git-tree-sha1 = "72aebe0b5051e5143a079a4685a46da330a40472"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.15"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.MultivariateStats]]
deps = ["Arpack", "Distributions", "LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI", "StatsBase"]
git-tree-sha1 = "816620e3aac93e5b5359e4fdaf23ca4525b00ddf"
uuid = "6f286f6a-111f-5878-ab1e-185364afe411"
version = "0.10.3"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "a0b464d183da839699f4c79e7606d9d186ec172c"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.3"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "cc0a5deefdb12ab3a096f00a6d42133af4560d71"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.1.2"

[[deps.NearestNeighbors]]
deps = ["Distances", "StaticArrays"]
git-tree-sha1 = "8a3271d8309285f4db73b4f662b1b290c715e85e"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.21"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Observables]]
git-tree-sha1 = "7438a59546cf62428fc9d1bc94729146d37a7225"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.5"

[[deps.OffsetArrays]]
git-tree-sha1 = "5e1897147d1ff8d98883cda2be2187dcf57d8f0c"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.15.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "38cb508d080d21dc1128f7fb04f20387ed4c0af4"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.3"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a9697f1d06cc3eb3fb3ad49cc67f2cfabaac31ea"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.16+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1346c9208249809840c91b26703912dff463d335"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.6+0"

[[deps.Optim]]
deps = ["Compat", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "ab7edad78cdef22099f43c54ef77ac63c2c9cc64"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.10.0"

    [deps.Optim.extensions]
    OptimMOIExt = "MathOptInterface"

    [deps.Optim.weakdeps]
    MathOptInterface = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6703a85cb3781bd5909d48730a67205f3f31a575"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.3+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "cc4054e898b852042d7b503313f7ad03de99c3dd"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.8.0"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "966b85253e959ea89c53a9abebbf2e964fbf593b"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.32"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "3b31172c032a1def20c98dae3f2cdc9d10e3b561"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.56.1+0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pidfile]]
deps = ["FileWatching", "Test"]
git-tree-sha1 = "2d8aaf8ee10df53d0dfb9b8ee44ae7c04ced2b03"
uuid = "fa939f87-e72e-5be4-a000-7fc836dbe307"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "db76b1ecd5e9715f3d043cec13b2ec93ce015d53"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.44.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"
weakdeps = ["REPL"]

    [deps.Pkg.extensions]
    REPLExt = "REPL"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "41031ef3a1be6f5bbbf3e8073f210556daeae5ca"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.3.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "StableRNGs", "Statistics"]
git-tree-sha1 = "3ca9a356cd2e113c420f2c13bea19f8d3fb1cb18"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.3"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "TOML", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "dae01f8c2e069a683d3a6e17bbae5070ab94786f"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.40.9"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "7e71a55b87222942f0f9337be62e26b1f103d3e4"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.61"

[[deps.Polynomials]]
deps = ["LinearAlgebra", "OrderedCollections", "RecipesBase", "Requires", "Setfield", "SparseArrays"]
git-tree-sha1 = "555c272d20fc80a2658587fb9bbda60067b93b7c"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "4.0.19"

    [deps.Polynomials.extensions]
    PolynomialsChainRulesCoreExt = "ChainRulesCore"
    PolynomialsFFTWExt = "FFTW"
    PolynomialsMakieCoreExt = "MakieCore"
    PolynomialsMutableArithmeticsExt = "MutableArithmetics"

    [deps.Polynomials.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    FFTW = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
    MakieCore = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
    MutableArithmetics = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "36d8b4b899628fb92c2749eb488d884a926614d3"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.3"

[[deps.PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "1101cd475833706e4d0e7b122218257178f48f34"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.4.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"
version = "1.11.0"

[[deps.PtrArrays]]
git-tree-sha1 = "1d36ef11a9aaf1e8b74dacc6a731dd1de8fd493d"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.3.0"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "492601870742dcd38f233b23c3ec629628c1d724"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.7.1+1"

[[deps.Qt6Declarative_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6ShaderTools_jll"]
git-tree-sha1 = "e5dd466bf2569fe08c91a2cc29c1003f4797ac3b"
uuid = "629bc702-f1f5-5709-abd5-49b8460ea067"
version = "6.7.1+2"

[[deps.Qt6ShaderTools_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "1a180aeced866700d4bebc3120ea1451201f16bc"
uuid = "ce943373-25bb-56aa-8eca-768745ed7b5a"
version = "6.7.1+1"

[[deps.Qt6Wayland_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6Declarative_jll"]
git-tree-sha1 = "729927532d48cf79f49070341e1d918a65aba6b0"
uuid = "e99dba38-086e-5de3-a5b1-6e4c66e897c3"
version = "6.7.1+1"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9da16da70037ba9d701192e27befedefb91ec284"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.11.2"

    [deps.QuadGK.extensions]
    QuadGKEnzymeExt = "Enzyme"

    [deps.QuadGK.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "StyledStrings", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"
version = "1.11.0"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "62389eeff14780bfe55195b7204c0d8738436d64"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.1"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "852bd0f55565a9e973fcfee83a84413270224dc4"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.8.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "58cdd8fb2201a6267e1db87ff148dd6c1dbd8ad8"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.5.1+0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "712fb0231ee6f9120e005ccd56297abbc053e7e0"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.8"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "c5391c6ace3bc430ca630251d02ea9687169ca68"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.2"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"
version = "1.11.0"

[[deps.ShiftedArrays]]
git-tree-sha1 = "503688b59397b3307443af35cd953a13e8005c16"
uuid = "1277b4bf-5013-50f5-be3d-901d8477a67a"
version = "2.0.0"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "f305871d2f381d21527c770d4788c06c097c9bc1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.2.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.11.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "64cca0c26b4f31ba18f13f6c12af7c85f478cfde"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.5.0"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StableRNGs]]
deps = ["Random"]
git-tree-sha1 = "83e6cce8324d49dfaf9ef059227f91ed4441a8e5"
uuid = "860ef19b-820b-49d6-a774-d7a799459cd3"
version = "1.0.2"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "0feb6b9031bd5c51f9072393eb5ab3efd31bf9e4"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.13"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"
weakdeps = ["SparseArrays"]

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["AliasTables", "DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "29321314c920c26684834965ec2ce0dacc9cf8e5"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.4"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "b423576adc27097764a90e163157bcfc9acf0f46"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.2"

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

    [deps.StatsFuns.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.StatsModels]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "Printf", "REPL", "ShiftedArrays", "SparseArrays", "StatsAPI", "StatsBase", "StatsFuns", "Tables"]
git-tree-sha1 = "9022bcaa2fc1d484f1326eaa4db8db543ca8c66d"
uuid = "3eaba693-59b7-5ba5-a881-562e759f1c8d"
version = "0.7.4"

[[deps.StatsPlots]]
deps = ["AbstractFFTs", "Clustering", "DataStructures", "Distributions", "Interpolations", "KernelDensity", "LinearAlgebra", "MultivariateStats", "NaNMath", "Observables", "Plots", "RecipesBase", "RecipesPipeline", "Reexport", "StatsBase", "TableOperations", "Tables", "Widgets"]
git-tree-sha1 = "3b1dcbf62e469a67f6733ae493401e53d92ff543"
uuid = "f3b207a7-027a-5e70-b257-86293d7955fd"
version = "0.15.7"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "725421ae8e530ec29bcbdddbe91ff8053421d023"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.4.1"

[[deps.StyledStrings]]
uuid = "f489334b-da3d-4c2e-b8f0-e476e12c162b"
version = "1.11.0"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.7.0+0"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableOperations]]
deps = ["SentinelArrays", "Tables", "Test"]
git-tree-sha1 = "e383c87cf2a1dc41fa30c093b2a19877c83e1bc1"
uuid = "ab02a1b2-a7df-11e8-156e-fb1833f50b87"
version = "1.2.0"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "598cd7c1f68d1e205689b1c2fe65a9f85846f297"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"
version = "1.11.0"

[[deps.TranscodingStreams]]
git-tree-sha1 = "0c45878dcfdcfa8480052b6ab162cdd138781742"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.3"

[[deps.Tricks]]
git-tree-sha1 = "6cae795a5a9313bbb4f60683f7263318fc7d1505"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.10"

[[deps.URIs]]
git-tree-sha1 = "cbbebadbcc76c5ca1cc4b4f3b0614b3e603b5000"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.2"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"
version = "1.11.0"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "c0667a8e676c53d390a09dc6870b3d8d6650e2bf"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.22.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "975c354fcd5f7e1ddcc1f1a23e6e091d99e99bc8"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.4"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "85c7811eddec9e7f22615371c3cc81a504c508ee"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+2"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "5db3e9d307d32baba7067b13fc7b5aa6edd4a19a"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.36.0+0"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.WebIO]]
deps = ["AssetRegistry", "Base64", "Distributed", "FunctionalCollections", "JSON", "Logging", "Observables", "Pkg", "Random", "Requires", "Sockets", "UUIDs", "WebSockets", "Widgets"]
git-tree-sha1 = "0eef0765186f7452e52236fa42ca8c9b3c11c6e3"
uuid = "0f1e0344-ec1d-5b48-a673-e5cf874b6c29"
version = "0.8.21"

[[deps.WebSockets]]
deps = ["Base64", "Dates", "HTTP", "Logging", "Sockets"]
git-tree-sha1 = "4162e95e05e79922e44b9952ccbc262832e4ad07"
uuid = "104b5d7c-a370-577a-8038-80a2059c5097"
version = "1.6.0"

[[deps.Widgets]]
deps = ["Colors", "Dates", "Observables", "OrderedCollections"]
git-tree-sha1 = "e9aeb174f95385de31e70bd15fa066a505ea82b9"
uuid = "cc8bc4a8-27d6-5769-a93b-9d913e69aa62"
version = "0.6.7"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c1a7aa6219628fcd757dede0ca95e245c5cd9511"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.0.0"

[[deps.WorkerUtilities]]
git-tree-sha1 = "cd1659ba0d57b71a464a29e64dbc67cfe83d54e7"
uuid = "76eceee3-57b5-4d4a-8e66-0e911cebbf60"
version = "1.6.1"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "b8b243e47228b4a3877f1dd6aee0c5d56db7fcf4"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.6+1"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "7d1671acbe47ac88e981868a078bd6b4e27c5191"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.42+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "56c6604ec8b2d82cc4cfe01aa03b00426aac7e1f"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.6.4+1"

[[deps.Xorg_libICE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "326b4fea307b0b39892b3e85fa451692eda8d46c"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.1.1+0"

[[deps.Xorg_libSM_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libICE_jll"]
git-tree-sha1 = "3796722887072218eabafb494a13c963209754ce"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.4+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "9dafcee1d24c4f024e7edc92603cedba72118283"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+3"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e9216fdcd8514b7072b43653874fd688e4c6c003"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.12+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "807c226eaf3651e7b2c468f687ac788291f9a89b"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.3+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "89799ae67c17caa5b3b5a19b8469eeee474377db"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.5+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "d7155fea91a4123ef59f42c4afb5ab3b4ca95058"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.6+3"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "6fcc21d5aea1a0b7cce6cab3e62246abd1949b86"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "6.0.0+0"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "984b313b049c89739075b8e2a94407076de17449"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.8.2+0"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll"]
git-tree-sha1 = "a1a7eaf6c3b5b05cb903e35e8372049b107ac729"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.5+0"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "b6f664b7b2f6a39689d822a6300b14df4668f0f4"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.4+0"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "a490c6212a0e90d2d55111ac956f7c4fa9c277a6"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.11+1"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c57201109a9e4c0585b208bb408bc41d205ac4e9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.2+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "1a74296303b6524a0472a8cb12d3d87a78eb3612"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.0+3"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "dbc53e4cf7701c6c7047c51e17d6e64df55dca94"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.2+1"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "04341cb870f29dcd5e39055f895c39d016e18ccd"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.4+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "ab2221d309eda71020cdda67a973aa582aa85d69"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.6+1"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "691634e5453ad362044e2ad653e79f3ee3bb98c3"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.39.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a63799ff68005991f9d9491b6e95bd3478d783cb"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.6.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "446b23e73536f84e8037f5dce465e92275f6a308"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.7+1"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "gperf_jll"]
git-tree-sha1 = "431b678a28ebb559d224c0b6b6d01afce87c51ba"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.9+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6e50f145003024df4f5cb96c7fce79466741d601"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.56.3+0"

[[deps.gperf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0ba42241cb6809f1a278d0bcb976e0483c3f1f2d"
uuid = "1a1c6b14-54f6-533d-8383-74cd7377aa70"
version = "3.1.1+1"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "522c1df09d05a71785765d19c9524661234738e9"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.11.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "e17c115d55c5fbb7e52ebedb427a0dca79d4484e"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.2+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.libdecor_jll]]
deps = ["Artifacts", "Dbus_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pango_jll", "Wayland_jll", "xkbcommon_jll"]
git-tree-sha1 = "9bf7903af251d2050b467f76bdbe57ce541f7f4f"
uuid = "1183f4f0-6f2a-5f1a-908b-139f9cdfea6f"
version = "0.2.2+0"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "141fe65dc3efabb0b1d5ba74e91f6ad26f84cc22"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.11.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8a22cf860a7d27e4f3498a0fe0811a7957badb38"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.3+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "ad50e5b90f222cfe78aa3d5183a20a12de1322ce"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.18.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "068dfe202b0a05b8332f1e8e6b4080684b9c7700"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.47+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "490376214c4721cdaca654041f635213c6165cb3"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+2"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "814e154bdb7be91d78b6802843f76b6ece642f11"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.6+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d5a767a3bb77135a99e433afe0eb14cd7f6914c3"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2022.0.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "63406453ed9b33a0df95d570816d5366c92b7809"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+2"
"""

# ╔═╡ Cell order:
# ╟─f269ff32-3e80-4548-97c6-d337827db7aa
# ╟─25699170-0fcd-11f0-2202-772042e9fd47
# ╟─80ec79c4-9971-47a2-a884-5417254a31fd
# ╟─ee84b436-12a5-4f87-8187-e064574632bb
# ╟─d1af481b-5a8e-469f-982a-401fbffd4590
# ╟─efda5197-b8ff-467a-8d32-6d9b7992a7ec
# ╟─2a36ce70-5dc5-401d-864e-9b5e479b7913
# ╟─72bfcd0b-47a4-4535-8a9d-2e059a6e86a1
# ╟─cc962df0-b6a0-481c-a95b-9265930d38d7
# ╟─30c74cb8-8da9-4521-8f4a-48319390f5e6
# ╟─6f370423-0a7e-4a20-9c69-81fe352020a0
# ╟─b07f485a-7168-4d87-acc2-026c1e2b4ad2
# ╟─34394d4f-5c13-4fa5-977e-1f0aee78dc2c
# ╟─8c1e2c97-381d-4c3b-af76-e269e46feb8c
# ╟─0ccfae55-6162-4c97-a04d-2dd3c5d9370e
# ╟─95432f3f-e4e5-48e4-abc2-522a1af4d402
# ╟─7d85f4cd-e4e8-4076-afc5-33daf90c08fe
# ╟─fb84175b-122e-4121-8498-085810627024
# ╟─3ae662a1-2c66-4810-9770-514cc93897ff
# ╟─077c4e59-79e9-4108-85c8-cf6e06458394
# ╟─c93a3123-22d0-4433-a6d3-34cc1395016a
# ╠═1f1cab74-5f8d-4ab4-a12e-5e53a89d9080
# ╟─63e51762-c84e-4839-9160-fb2dcf82e8ac
# ╟─d62f9aab-db3d-49f3-9844-ed9b091903ef
# ╟─bab951b8-227f-4c89-a4a3-08f65a16fc35
# ╟─627e2f9c-65d9-4645-a1a4-48c612876bd5
# ╠═6a9d29ad-ecb7-4f9a-a86f-7e67e658a7de
# ╟─711fdcf5-333a-44b9-a821-30cea343d89c
# ╟─bb8551a1-8e35-49a4-8aa8-1fc89167ae7e
# ╟─57fa7cfd-8521-46ca-b2fd-70dc506e2593
# ╟─9a619a6a-0814-4c21-a3be-50041c3e198f
# ╟─1ad1dc9b-f12d-4e9a-8ea8-c3c85ab738c3
# ╟─4d4b096d-a87a-482f-a19c-7918e7af47fa
# ╠═ddba2999-edab-4899-86d7-44c822e33ac1
# ╟─ab4b6e35-bbd6-4c5d-ae6c-cf3cb2092a85
# ╠═5e980df1-7333-4210-97c1-7aa3f0de9a67
# ╠═50dde6f7-8c8b-439d-9dbe-498ef41aa5b3
# ╟─d7213a5c-f663-420d-99fd-1ea15758c4bc
# ╟─868ae38e-5355-42cf-b645-1ae6d6b443c8
# ╟─6e0c44d6-af43-4dbd-ab2f-d36053e2559a
# ╠═4d5dbdaa-7859-45fe-ba01-5d783400543c
# ╟─93763c84-79f2-4cca-9a97-a871361c92a3
# ╟─f7c435ef-deb8-49b5-8406-ec0b36116aac
# ╟─bc5edd76-7aa4-494c-b78e-2949a032b0a9
# ╟─f0807783-6915-49e5-a49f-eb770695b0ce
# ╟─7f1f1fc5-1077-4f56-94de-6a56521ad0e1
# ╟─f5902359-7f0b-4244-8ebd-191358e3a227
# ╟─de2c68cb-640c-4fbd-83f8-eb8ef48104c0
# ╟─efa3f064-bc0f-491b-b06a-d8cd6d42a174
# ╠═54ce63de-9640-4c00-b685-f8163f5f4bd8
# ╟─cb94e83a-2e31-45d9-b4e0-b352d4595e99
# ╠═dab5f085-ac90-4342-b2eb-60b097320de9
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002

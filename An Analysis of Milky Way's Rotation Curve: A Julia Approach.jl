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

# ╔═╡ e25a61ff-03e8-4cd4-bec2-52a00834b011
begin
	import Pkg
	
	# Install required packages
	Pkg.add([
	    "CSV",
		"LinearAlgebra",
	    "DataFrames", 
	    "Plots"])
end

# ╔═╡ 0b0fbb27-6b38-4ab7-8b14-18a3023698b6
begin
	using CSV, DataFrames
	
	file_path = "DataSets/TwentyKTWO.csv"
	
	# Load the Gaia data
	df = CSV.read(file_path, DataFrame)
	# Drop rows with missing RA or Dec values
	df1 = dropmissing(df, [:ra, :dec])
	
	# Define the function to filter stars based on their proximity to the Galactic plane
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
	
	# Apply filtering: only stars with non-missing parallax, positive parallax, and proper motions
	filtered_df = filter_galactic_plane(df1)
	
	# Further filtering for stars with non-missing parallax, pmra, and pmdec
		N = filter(row -> !ismissing(row.parallax) && row.parallax > 0 && !ismissing(row.pmra) && !ismissing(row.pmdec), filtered_df)

end


# ╔═╡ c6d7fed2-e29b-4396-9d19-c07e96a02bd9
begin
	using LinearAlgebra
	# Constants
	
	# Function to convert RA, Dec, d to Galactocentric radius and velocities
	function galactocentric_transform(ra, dec, d, v_ra, v_dec, v_r)
	    # Convert to radians
	    α = deg2rad(ra)
	    δ = deg2rad(dec)
	    
	    # Convert Equatorial to Cartesian coordinates relative to Sun
	    x = d * cos(δ) * cos(α)
	    y = d * cos(δ) * sin(α)
	    z = d * sin(δ)
	    
	    # Rotation matrix from Equatorial to Galactic coordinates
	    R = [-0.05487556 -0.87343709 -0.48383502;
	         +0.49410943 -0.44482963 +0.74698224;
	         -0.86766615 -0.19807637 +0.45598378]
	    
	    # Galactic Cartesian coordinates relative to Sun
	    galactic_coords = R * [x, y, z]
	    X, Y, Z = galactic_coords
	    
	    # Compute Galactocentric radius (ensure it's in kpc)
	    R_G = sqrt(X^2 + Y^2 + Z^2)  # kpc
	    
	    # Velocity transformation
	    v_xyz = R * [v_ra, v_dec, v_r]
	    v_X, v_Y, v_Z = v_xyz
	    
	    # Compute orbital velocity (V_phi in the Galactic plane)
	    V_phi = (X * v_Y - Y * v_X) / sqrt(X^2 + Y^2)
	    
	    return R_G, v_X, v_Y, v_Z, V_phi
	end
	
	# Apply transformation and filter
	results = [galactocentric_transform(row.ra, row.dec, row.distance_kpc, row.pmra, row.pmdec, row.radial_velocity) for row in eachrow(N)]

	N[!, :Galactocentric_Radius] = [r[1] for r in results]
	N[!, :Vx] = [r[2] for r in results]
	N[!, :Vy] = [r[3] for r in results]
	N[!, :Vz] = [r[4] for r in results]
	N[!, :V_phi] = [r[5] for r in results]  # Orbital velocity
	
	N_New = filter(row -> row.Galactocentric_Radius > 0, N)  # Keeping only positive distances
end

# ╔═╡ 103b6c4c-a50c-4c61-8499-d78df265fae1
begin
	using PlutoUI
	md"""
	**Adjust X-Axis Scale (View Window)  
	Max R (kpc):** $(@bind xmax Slider(0:1:50; default=25, show_value=true))
	
	Galactic Enclosed Mass / Solar Masses:** $(@bind mass_slider Slider(1e10:1e10:1e12, show_value=true, default=1e11))
		"""
end

# ╔═╡ 059cb579-672a-4ff9-9477-c3decc2c785e
begin
	using Plots
	
	scatter(N_New.Galactocentric_Radius, N_New.V_phi.+220, marker=:circle, xlabel="Galactocentric Radius (kpc)", ylabel="Orbital Velocity (km/s)", title="Orbital Velocity vs. Galactocentric Radius", legend=true, label="Observed Object", xlims=(0, xmax), ylims=(0, 500))

	G = 4.302e-6  # kpc * (km/s)^2 / Msun
	
	keplerian_velocity(r, M) = sqrt(G * M / r)
	
	
	### Galactic rotation curve
	r_values = 0.1:0.1:xmax # Range of distances (kpc)
	v_kepler = keplerian_velocity.(r_values, mass_slider)  # Keplerian velocity based on current mass
	
	### Generating scatter plot and overlaying Keplerian curve
	
	plot!(r_values, v_kepler, label="Keplerian Curve", legend = true, linewidth=2, color=:red)
end

# ╔═╡ d7c0a6ae-9a3b-4955-bda4-7325d15f08d9
md"""
#### Parag Chettri ~ Santiago Berumen ~ Isaac Whitson
# An Analysis of Milky Way's Rotation Curve: A Julia Approach
## Data from Gaia DR 2
"""

# ╔═╡ 8d3e1dc6-7d84-42e3-b228-9cf73313fc2b
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
	
	We will choose a value of $b = 5^\circ$.
	"""
end

# ╔═╡ d3c12dbb-0fc1-4d31-bf5e-5812c5e51fb8
begin
	md"""
	### Distance Calculation from Parallax 
	"""
end

# ╔═╡ bceadb12-7d88-4a8b-aae0-a696366627ce
begin
	md""" 
	
	The distance $D$ in kiloparsecs $(kpc)$ is calculated using the parallax $p$ in milliarcseconds $(mas)$:  
	
	$D = \frac{1}{p}$
	
	where: $D$ is in kiloparsecs $(kpc)$  
	
	where: $p$ is the parallax in milliarcseconds $(mas)$ 

	The array of distances in $kpc$ is shown below.
	"""
end

# ╔═╡ 0f93be8e-c045-4cba-8755-44bf01faad1e
begin
	N.distance_kpc_error = 1.0 ./ (N.parallax_error)
	N.distance_kpc = 1.0 ./ N.parallax
end

# ╔═╡ 26c44b6e-0e5a-47d1-a2d8-52cdc8a8b3f9
begin
	md"""
	### Finding the stars' Tangential Velocity
	"""
end

# ╔═╡ 37675bef-ab4b-4b7a-a946-102148c6dce6
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

# ╔═╡ ae23b057-cb23-408f-801d-db86b29ce617
begin
	# Compute total proper motion in mas/yr
    N.proper_motion = sqrt.(N.pmra.^2 .+ N.pmdec.^2)
	N.proper_motion_error = sqrt.(N.pmra_error.^2 .+ N.pmdec_error.^2)
    # Compute tangential velocity (km/s)
	N.tangential_velocity_error = 4.74 .* N.proper_motion_error .* N.distance_kpc_error
	N.tangential_velocity = 4.74 .* N.proper_motion .* N.distance_kpc

end

# ╔═╡ 5580ffda-b8fa-4504-830e-588c91bcdbda
begin
	md"""
	### True Velocity from Radial and Tangential Velocities
	"""
end

# ╔═╡ 4cb64ec9-4c1c-47ec-857d-45e764df2e56
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

# ╔═╡ db27aa1b-23da-4ace-bfda-9a34ecf4554a
begin
 # Compute true velocity (km/s), filtering out missing radial velocities
	N.true_velocity_error = sqrt.(N.tangential_velocity_error.^2 .+ coalesce.(N.radial_velocity_error, 0.0).^2)
	N.tangential_velocity_corrected = N.tangential_velocity
	N.true_velocity = sqrt.(N.tangential_velocity_corrected.^2 .+ coalesce.(N.radial_velocity, 0.0).^2)
end

# ╔═╡ 947e2007-d7e9-4049-b324-0b41d13c80b3
begin
	md"""
	Now there is something still wrong with our data. 
	
	Gaia, our source, does not correct for the orbital velocity of our sun and the solar system. This means, in our frame of reference, all the stars within our neighborhood with similar velocities seem to move with near-zero velocities. 

	We must now find a way to correct our data, particularly by shifting the rest frame to the Galactic Center and have our independant value of Distance be the Galactocentric Radius instead!
	
	"""
end

# ╔═╡ 58c42bec-6967-4023-8d52-eea5376ec0b8
begin
	md"""
	# Translating our data to a Galactocentric Rest Frame
	"""
end

# ╔═╡ 8ee908d3-f28e-4146-8a10-811fde6938cc
begin
	md"""
	## **Transformation Function**
	The function `galactocentric_transform` converts equatorial coordinates (RA, Dec, distance) and velocities to Galactocentric values.
	
	### **1. Coordinate Conversion**
	Given:
	- Right Ascension ($\alpha$) in degrees
	- Declination ($\delta$) in degrees
	- Distance ($d$) in kpc
	- Proper motions ($v_{ra}$, $v_{dec}$) and radial velocity ($v_r$)
	
	We convert RA and Dec to radians:
	$\alpha = \deg2rad(ra), \quad \delta = \deg2rad(dec)$
	
	Using spherical-to-Cartesian transformation, the position $(x, y, z)$ in the equatorial frame is:
		
	$\begin{aligned}
		x &= d \cos(\delta) \cos(\alpha) \\
		y &= d \cos(\delta) \sin(\alpha) \\
		z &= d \sin(\delta)
	\end{aligned}$
	
	### **2. Rotation to Galactic Frame**
	The transformation matrix $R$ converts equatorial coordinates to Galactic coordinates:
	
	$R = \begin{bmatrix}
	    -0.05487556 & -0.87343709 & -0.48383502 \\
	    +0.49410943 & -0.44482963 & +0.74698224 \\
	    -0.86766615 & -0.19807637 & +0.45598378
	\end{bmatrix}$

	
	Multiplying $R$ by the position vector:
	
	$\begin{bmatrix} X \\ Y \\ Z \end{bmatrix} = R \cdot \begin{bmatrix} x \\ y \\ z \end{bmatrix}$
	
	### **3. Computing Galactocentric Radius**
	The Galactocentric radius is given by:
	
	$R_G = \sqrt{X^2 + Y^2 + Z^2} \text{ kpc}$
	
	### **4. Velocity Transformation**
	Applying the same rotation matrix $R$ to the velocity components:
	
	$\begin{bmatrix} v_X \\ v_Y \\ v_Z \end{bmatrix} = R \cdot \begin{bmatrix} v_{ra} \\ v_{dec} \\ v_r \end{bmatrix}$
	
	### **5. Computing Orbital Velocity**
	The azimuthal (orbital) velocity $V_{\phi}$ in the Galactic plane is derived using:
	
	$V_{\phi} = \frac{X v_Y - Y v_X}{\sqrt{X^2 + Y^2}}$
	
	## **Applying the Transformation**
	The function is applied to each row in dataset `N`, extracting transformed values:
	
	$N[!, :Galactocentric\_Radius] = R_G$
	
	$N[!, :Vx] = v_X, \quad N[!, :Vy] = v_Y, \quad N[!, :Vz] = v_Z, \quad N[!, :V_\phi] = V_\phi$
	
	Finally, a filtering step ensures only positive Galactocentric distances are retained:
	
	$N\_New = \{ R_G > 0 \}$
	
	This completes the transformation of equatorial coordinates and velocities into a Galactocentric reference frame.
	"""


end

# ╔═╡ 52be580b-0add-4c32-aba1-719b5898c5fa
begin
	md"""
	This calculation of $V_{\phi}$ is jus the given deviation from the Galactic mean orbital velocity which we can assign as a normalizing value on our graph for disk stars as +220 km/s [(Stellar Kinematics)](https://en.wikipedia.org/wiki/Stellar_kinematics).
	"""
end

# ╔═╡ ddd1d5c6-cbfd-48e4-89b4-bfb3767d7c14
begin
	md"""
	### The Keplerian Orbit
	"""
end

# ╔═╡ fcbbfa5d-bbf3-480b-b1b0-c4b835525113
begin
	md"""

	As we know, Keplerian orbits are a fundamental part of Newtonian Physics, which we can overlay to check againt in this case. 

	$V_{\text{Keplerian}}(R) = \sqrt{\frac{GM_{\text{enclosed}}}{R}}$

	Where: $V_{\text{Keplerian}}$ is the orbital velocity in $km/s$

	Where: $M_{\text{enclosed}}$ is the mass within the galactic center that other galactic objects orbit.

	Where: $R$ is the distance from that mass.

	"""
end
	

# ╔═╡ 2b7a8b7e-a6d1-4360-b155-54c1bc771ee1
begin
	md"""
	Needless to say, the rotation of the galaxy is clearly not bound just to Keplarian parameters.
	"""
end

# ╔═╡ 81f8893e-4c4b-4318-8e3d-70259cf4e044
begin
	md"""
	# Our Prior - Milky Way Rotation Curve Components
	
	The total rotation curve of the Milky Way is the quadrature sum of multiple contributions from the bulge, stellar disk, HI layer, H2 layer, and dark halo:
	
	
	$V_{\text{total}}(r) = \sqrt{V_{\text{bulge}}^2 + V_{\text{disk}}^2 + V_{\text{HI}}^2 + V_{\text{H2}}^2 + V_{\text{halo}}^2}$
	
	#### **1. Bulge (Blue Line)**
	The bulge follows a steep Keplerian-like rise and then quickly levels off:
	
	$V_{\text{bulge}}(r) = V_0 \frac{r}{(r^2 + a^2)^{3/4}}$
	
	where: $V_0$ is a scaling factor $(~300 km/s)$,
	where: $a$ is a scale radius $(~0.5–1 kpc)$.
	
	#### **2. Stellar Disk (Green Line)**
	The stellar disk follows an exponential profile:
	
	$V_{\text{disk}}(r) = V_0 \left(\frac{r}{R_d}\right) e^{-r/(2 R_d)}$
	
	where: $R_d$ is the disk scale length $(~3 kpc)$,
	where: $V_0$ is a normalization constant $(~200 km/s)$.
	
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
	
	where: $V_{\infty}$ is the asymptotic velocity $(~220 km/s)$,
	where: $R_H$ is the halo scale radius $(~20 kpc)$,
	where: $c$ is a core radius $(~10 kpc)$.

	#### **4. Final Estimation (Black Line)**
	
	$V_{\text{total}}(r)$

	This is our final estimation as a weighted sum of the given distributions.

	Source: [(Ueshima, 2010)](https://www-sk.icrr.u-tokyo.ac.jp/xmass/publist/ueshima_PhD.pdf)
	"""
end


# ╔═╡ 9f4404c3-ca98-4be7-9c8f-187ff582250e
begin
	md"""
	### Prior over the Velocity of Sources
	"""
end

# ╔═╡ 3e5f61b8-edbc-41b2-a047-0656bf419bdb
begin
	md"""
	**Adjust X-Axis Scale (View Window)  
	Max R (kpc):** $(@bind xmax_2 Slider(0:1:50; default=25, show_value=true))
	"""
end

# ╔═╡ 97e6ab22-be14-4b2a-9032-f5432839ac23
begin
	
	# Define the rotation curve components
	function V_bulge(r, V0, a)
	    return V0 * (r ./ (r.^2 .+ a^2).^(3/4))
	end
	
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
	
	# Define the range of r (distance in kpc)
	r = 0:0.1:xmax_2  # 0 to user selected value for kpc
	
	# Parameters for the components (from the description)
	V0_bulge = 300
	a_bulge = 1
	
	V0_disk = 200
	R_d = 3
	
	V0_HI = 100
	b_HI = 5
	
	V0_H2 = 80
	R_H2 = 4
	
	V_infinity_halo = 220
	R_H = 20
	c_halo = 10
	
	# Calculate each component
	V_bulge_curve = V_bulge(r, V0_bulge, a_bulge)
	V_disk_curve = V_disk(r, V0_disk, R_d)
	V_HI_curve = V_HI(r, V0_HI, b_HI)
	V_H2_curve = V_H2(r, V0_H2, R_H2)
	V_halo_curve = V_halo(r, V_infinity_halo, c_halo)
	
	# Total rotation curve (quadrature sum)
	V_total_curve = sqrt.(V_bulge_curve.^2 + V_disk_curve.^2 + V_HI_curve.^2 + V_H2_curve.^2 + V_halo_curve.^2)
	
	# Scatter plot with filtered data (Assuming `N_filtered` is already defined)
	scatter(N_New.Galactocentric_Radius, N_New.V_phi .+ 220, marker=:circle, label="Data Points", xlabel="Galactocentric Radius (kpc)", ylabel="Orbital Velocity (km/s)", title="Orbital Velocity vs. Galactocentric Radius", legend=false, alpha=0.1 , xlims=(0, xmax_2), ylims=(0, 500))
	
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


# ╔═╡ 057a643b-128e-461e-9697-8b8344ce6331
begin
	md"""
	**Why are there so many data points between 0 and 7 kpc?**
	
	The answer: The galactic center bulge is dense with stars and we are about 8 kpc from the center of the galaxy.
	"""
end

# ╔═╡ 6b1f254e-9926-40af-b1b9-4a30603ff595
begin
	md"""
	# NEXT STEPS

	Parag - Fitted the Data with Swag

	Santi - TwentyK[NUMBER] files. Make this user-selectable.

	Isaac - List all statistical things from the labs which we can use here.

	
	
	"""
end

# ╔═╡ 729d9763-c8ca-4c2d-b208-33373bf66bf5
begin
	md"""
	
	### **Citations**


	
	[Ueshima, K. (2010). Study of pulse shape discrimination and low background techniques for liquid xenon dark matter detectors. Department of Physics, School of Science, University of Tokyo.](https://www-sk.icrr.u-tokyo.ac.jp/xmass/publist/ueshima_PhD.pdf)

	"""
end

# ╔═╡ Cell order:
# ╟─d7c0a6ae-9a3b-4955-bda4-7325d15f08d9
# ╟─0b0fbb27-6b38-4ab7-8b14-18a3023698b6
# ╟─8d3e1dc6-7d84-42e3-b228-9cf73313fc2b
# ╟─d3c12dbb-0fc1-4d31-bf5e-5812c5e51fb8
# ╟─bceadb12-7d88-4a8b-aae0-a696366627ce
# ╟─0f93be8e-c045-4cba-8755-44bf01faad1e
# ╟─26c44b6e-0e5a-47d1-a2d8-52cdc8a8b3f9
# ╟─37675bef-ab4b-4b7a-a946-102148c6dce6
# ╟─ae23b057-cb23-408f-801d-db86b29ce617
# ╟─5580ffda-b8fa-4504-830e-588c91bcdbda
# ╟─4cb64ec9-4c1c-47ec-857d-45e764df2e56
# ╟─db27aa1b-23da-4ace-bfda-9a34ecf4554a
# ╟─947e2007-d7e9-4049-b324-0b41d13c80b3
# ╟─58c42bec-6967-4023-8d52-eea5376ec0b8
# ╟─8ee908d3-f28e-4146-8a10-811fde6938cc
# ╟─c6d7fed2-e29b-4396-9d19-c07e96a02bd9
# ╟─52be580b-0add-4c32-aba1-719b5898c5fa
# ╟─ddd1d5c6-cbfd-48e4-89b4-bfb3767d7c14
# ╟─fcbbfa5d-bbf3-480b-b1b0-c4b835525113
# ╟─059cb579-672a-4ff9-9477-c3decc2c785e
# ╟─103b6c4c-a50c-4c61-8499-d78df265fae1
# ╟─2b7a8b7e-a6d1-4360-b155-54c1bc771ee1
# ╟─81f8893e-4c4b-4318-8e3d-70259cf4e044
# ╟─9f4404c3-ca98-4be7-9c8f-187ff582250e
# ╟─97e6ab22-be14-4b2a-9032-f5432839ac23
# ╟─3e5f61b8-edbc-41b2-a047-0656bf419bdb
# ╟─057a643b-128e-461e-9697-8b8344ce6331
# ╟─6b1f254e-9926-40af-b1b9-4a30603ff595
# ╟─729d9763-c8ca-4c2d-b208-33373bf66bf5
# ╟─e25a61ff-03e8-4cd4-bec2-52a00834b011

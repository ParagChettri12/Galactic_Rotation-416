println("Making sure Pluto is installed")
import Pkg; 
Pkg.add(name="Pluto", version="0.20.4");
import Pluto; 

println("Installing dashboard.jl")
Pluto.activate_notebook_environment("../dashboard.jl"); 
Pkg.instantiate(); 

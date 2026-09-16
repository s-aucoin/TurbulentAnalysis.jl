module TurbulentAnalysis

using ExtraStats       # for extra stats and least squares fitting functions
using DataMethods      # Extra data processing and analysis functions
using StatsBase        # for Statistics
using HypothesisTests  # For even more stats
using LinearAlgebra    # For linear algebra #
using Measurements     # For dealing with uncertainty
using Dates            # For dealing with time
import Base.Threads.@threads # For parallel computing


include("velocity_structure_functions.jl") # functions for calculating structure functions
include("quality_control.jl") # functions for applying quality control to the turbulence calculations
include("misc.jl") # miscellaneous functions for turbulence analysis

end
